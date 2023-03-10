
'''
# Name: MGI_DiseaseModel.py
#
# List of human genes with one-to-one relationship with mouse ortholog
# 
# Rules:
# 1: if human gene has one-to-one relationship with the mouse ortholog
# 2: those that have a DO annotation from the Alliance gene-disease relationship (_annottype_key = 1022)
# 3: report 1 only: exclude: NOT annotations (exclude qualifier = 1614157)
#
# Report #1 columns:
# 1:  Human gene symbol
# 2:  Human gene name
# 3:  HGNC id for human gene
# 4:  DO term name associated with human gene (if >1, show in separate row, duplicate other columns 1-3, and 5)
# 5:  DO term ID associated with human gene
# 6:  Mouse genotype IDs
# 7:  Mouse gene (related to human gene in column 1)
# 8:  Mouse MGI ID
# 9:  Facilities
# 10: Repository ID
#
# Report #2 columns:
# 1:  DO term name associated with human gene (if >1, show in separate row, duplicate other columns 1-3, and 5)
# 2:  DO term ID associated with human gene
# 3:  NOT model : NOT or blank
# 4:  Allele Pairs 
# 5:  Strain Background 
# 6:  Allele symbol from column 4 (exclude wild type allele and allele with attribute = recombinase or transactivator)
# 7:  Allele MGI ID 
# 8:  Total number of allele references for allele in column 6/7
# 9:  Repository ID for allele
# 10: Mouse genotype RR id
# 11: Mouse gene symbol from allele in column 6/7
# 12: Mouse gene ID
# 13: Repository ID for gene 
#
# History:
#
# lec   01/23/2023
#       wts2-1095/fl2-149/Mouse models of human disease report (1)
#       wts2-1120/fl2-204/Mouse models of human disease report (2)
#
'''

import sys
import os
import csv
import db
import reportlib

db.setTrace()

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

IMSR_CSV = os.environ['IMSR_STRAINS_CSV']
genotypeidLookup = {}
rridLookup = {}
facilityLookup = {}
repositoryAlleleLookup = {}
repositoryGeneLookup = {}

def initialize():
        """
        build lookups
        start building the mouse/human one-to-one
        """

        global genotypeidLookup
        global rridLookup
        global facilityLookup
        global repositoryAlleleLookup
        global repositoryGeneLookup

        # cache DO/Genotype (1020)
        results = db.sql('''
                select a.accID as doID, aa.accID, aa._LogicalDB_key
                from VOC_Annot vg, ACC_Accession a, ACC_Accession aa
                where a._LogicalDB_key = 191
                and a._Object_key = vg._Term_key
                and vg._AnnotType_key = 1020
                and vg._Object_key = aa._Object_key
                and aa._MGIType_key = 12
                and aa._LogicalDB_key in (1, 179)
                order by doID
        ''', 'auto')
        for r in results:
                key = r['doID']
                value = r['accID']
                if r['_LogicalDB_key'] == 1:
                        if key not in genotypeidLookup:
                                genotypeidLookup[key] = []
                        genotypeidLookup[key].append(value)
                else:
                        if key not in rridLookup:
                                rridLookup[key] = []
                        rridLookup[key].append(value)

        #print(genotypeidLookup)
        #print(rridLookup)

        # cache marker & provider
        csvfile = open(IMSR_CSV, 'r')
        reader = csv.reader(csvfile)
        for r in reader:
                allele_ids = r[0]
                marker_ids = r[1]
                provider = r[3]
                strain_id = r[5].replace(':EuMMCR','EuMMCR')
                for id in allele_ids.split(','):
                    repositoryAlleleLookup.setdefault(id,[]).append(strain_id)
                for id in marker_ids.split(','):
                    facilityLookup.setdefault(id,[]).append(provider)
                    repositoryGeneLookup.setdefault(id,[]).append(strain_id)

        # unique and sort the provider list
        for id, facilities in list(facilityLookup.items()):
                facilityLookup[id] = list(set(facilities))
                facilityLookup[id].sort()
        #print(facilityLookup)
        #print(repositoryAlleleLookup['MGI:1855953'])
        #print(repositoryGeneLookup)

        # mouse/human orthologs
        db.sql('''
        select distinct
               m1._Marker_key as m_Marker_key,
               m1.symbol as msymbol,
               m1.name as mname,
               a.accID as markerID,
               m2._Marker_key as h_Marker_key,
               m2.symbol as hsymbol,
               m2.name as hname
        into temporary table homology
        from MRK_Cluster mc,
             MRK_ClusterMember h1,
             MRK_ClusterMember h2,
             MRK_Marker m1, ACC_Accession a,
             MRK_Marker m2, MRK_Status ms
        where mc._ClusterSource_key = 75885739
        and mc._Cluster_key = h1._Cluster_key
        and h1._Marker_key = m1._Marker_key
        and m1._Organism_key = 1
        and m1._Marker_Status_key = ms._Marker_Status_key
        and h1._Cluster_key = h2._Cluster_key
        and h2._Marker_key = m2._Marker_key
        and m2._Organism_key = 2
        and m1._Marker_key = a._Object_key
        and a._MGIType_key = 2
        and a._LogicalDB_key = 1
        and a.preferred = 1
        ''', None)
        db.sql('create index homology_idx1 on homology(m_Marker_key)', None)
        db.sql('create index homology_idx2 on homology(h_Marker_key)', None)

        # one-to-one
        db.sql('''
        select h.m_Marker_key
        into temporary table homology_mouse
        from homology h
        group by h.m_Marker_key having count(*) = 1
        ''', None)
        db.sql('create index mouse_idx on homology_mouse(m_Marker_key)', None)

        db.sql('''
        select h.h_Marker_key
        into temporary table homology_human
        from homology h
        group by h.h_Marker_key having count(*) = 1
        ''', None)
        db.sql('create index human_idx on homology_human(h_Marker_key)', None)

def processReport1():
        """
        final mouse/human one-to-one with DO/Human annotations (1022)
        excluding "NOT" qualifiers (_vocab_key = 53, _qualifier_key = 1614157 = NOT)
        i.e. only include _qualifier_key = 1614158
        """

        # report 1
        fp1 = reportlib.init('MGI_DiseaseGeneModel', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
        fp1.write('Human gene symbol' + TAB)
        fp1.write('Human gene name' + TAB)
        fp1.write('HGNC id for human gene' + TAB)
        fp1.write('DO term name associated with human gene' + TAB)
        fp1.write('DO term ID associated with human gene' + TAB)
        fp1.write('Mouse genotype IDs' + TAB)
        fp1.write('Mouse gene' + TAB)
        fp1.write('Mouse MGI ID' + TAB)
        fp1.write('Facilities for mouse gene' + TAB)
        fp1.write('Repository ID' + CRT)

        # final orthology results with DO associations and qualifier != NOT
        db.sql('''
        select c.*, t.term as doTerm, a.accID as doID
        into temporary table results1
        from homology_mouse m, homology_human h, homology c, VOC_Annot va, VOC_Term t, ACC_Accession a
        where c.m_Marker_key = m.m_Marker_key
        and c.h_Marker_key = h.h_Marker_key
        and c.h_Marker_key = va._Object_key
        and va._AnnotType_key = 1022
        and va._Qualifier_key = 1614158
        and va._Term_key = t._Term_key
        and t._Term_key = a._Object_key
        and a._MGIType_key = 13
        and a.preferred = 1
        and a._LogicalDB_key = 191
        ''', 'auto')
        db.sql('create index results1_idx1 on results1(m_Marker_key)',  None)
        db.sql('create index results1_idx2 on results1(h_Marker_key)',  None)

        # cache HGNC acc ids
        hgncLookup = {}
        results = db.sql('''
        select distinct r.h_Marker_key, a.accID
        from results1 r, ACC_Accession a
        where r.h_Marker_key = a._Object_key
        and a._MGIType_key = 2
        and a._LogicalDB_key = 64
        and a.preferred = 1
        order by r.h_Marker_key, a.accID
        ''', 'auto')
        for r in results:
                key = r['h_Marker_key']
                value = r['accID']
                if key not in hgncLookup:
                        hgncLookup[key] = []
                hgncLookup[key].append(value)
        #print(hgncLookup)

        results = db.sql('select * from results1 order by hsymbol, doTerm', 'auto')
        for r in results:

                # 1:  Human gene symbol
                # 2:  Human gene name
                fp1.write(r['hsymbol'] + TAB)
                fp1.write(r['hname'] + TAB)

                # 3:  HGNC id for human gene
                key = r['h_Marker_key']
                if key in hgncLookup:
                        fp1.write('|'.join(hgncLookup[key]))
                fp1.write(TAB)

                # 4:  DO term associated with human gene
                # 5:  DO term ID associated with human gene
                doID = r['doID']
                fp1.write(r['doTerm'] + TAB)
                fp1.write(r['doID'] + TAB)

                # 6:  Mouse genotype IDs
                if doID in genotypeidLookup:
                        fp1.write('|'.join(genotypeidLookup[doID]))
                fp1.write(TAB)

                # 7:  Mouse gene (related to human gene in column 1)
                # 8:  Mouse MGI ID
                markerID = r['markerID']
                fp1.write(r['msymbol'] + TAB)
                fp1.write(markerID + TAB)

                # 9:  Facilities with strains carrying mutations in this gene
                if markerID in facilityLookup:
                        fp1.write('|'.join(facilityLookup[markerID]))
                fp1.write(TAB)

                # 10: ID from repository-stock number
                if markerID in repositoryGeneLookup:
                        fp1.write('|'.join(repositoryGeneLookup[markerID]))
                fp1.write(CRT)
                
        reportlib.finish_nonps(fp1)

def processReport2():
        """
        final mouse/human one-to-one with DO/Human annotations (1022)
        """

        # 1:  DO term name associated with human gene (if >1, show in separate row, duplicate other columns 1-3, and 5)
        # 2:  DO term ID associated with human gene
        # 3:  NOT model : NOT or blank
        # 4:  Allele Pairs 
        # 5:  Strain Background 
        # 6:  Allele symbol (exclude wild type allele and allele with attribute = recombinase or transactivator)
        # 7:  Allele MGI ID 
        # 8:  Total number of allele references for allele
        # 9:  Repository ID for allele
        # 10: Mouse genotype RR id
        # 11: Mouse gene symbol from allele
        # 12: Mouse gene ID
        # 13: Repository ID for gene

        # report 2
        fp2 = reportlib.init('MGI_DiseaseMouseModel', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
        fp2.write('DO term name associated with human gene' + TAB)
        fp2.write('DO term ID associated with human gene' + TAB)
        fp2.write('NOT model' + TAB)
        fp2.write('Allele Pairs' + TAB)
        fp2.write('Strain Background' + TAB)
        fp2.write('Allele symbol from col 4' + TAB)
        fp2.write('Allele MGI ID' + TAB)
        fp2.write('Total number of allele references' + TAB)
        fp2.write('Repository ID from allele' + TAB)
        fp2.write('Mouse genotype RR id' + TAB)
        fp2.write('Marker symbol from col 4' + TAB)
        fp2.write('Marker MGI ID' + TAB)
        fp2.write('Repository ID from gene' + CRT)

        # final orthology results where DO associations exist (qualifier != NOT)
        # where human marker -> DO/Human Marker (1022) exists
        db.sql('''
        select distinct c.m_marker_key, c.msymbol, c.markerID, t1._term_key, t1.term as doTerm, t2.term as qualifierTerm, a.accID as doID
        into temporary table results2
        from homology_mouse m, homology_human h, homology c, VOC_Annot va, VOC_Term t1, VOC_Term t2, ACC_Accession a
        where c.m_Marker_key = m.m_Marker_key
        and c.h_Marker_key = h.h_Marker_key
        and c.h_Marker_key = va._Object_key
        and va._AnnotType_key = 1022
        and va._Term_key = t1._Term_key
        and va._Qualifier_key = t2._Term_key
        and va._Term_key = a._Object_key
        and a._MGIType_key = 13
        and a.preferred = 1
        and a._LogicalDB_key = 191
        order by t1.term
        ''', 'auto')
        db.sql('create index results2_idx1 on results2(m_Marker_key)',  None)

        # cache DO/Genotype (1020), accids
        doAlleleLookup = {}
        results = db.sql('''
        select distinct r.m_Marker_key, g._genotype_key,
                aa._allele_key, aa.symbol as alleleSymbol, aaa.accID as alleleID, n.note as allelePairs, aa.isWildType, s.strain,
                case when exists (select 1 from VOC_Annot va
                        where aa._Allele_key = va._Object_key
                        and va._AnnotType_key = 1014
                        and va._Term_key in (11025588,13289567)
                        ) 
                        then 1 else 0 end as skipAllele
        from results2 r, VOC_Annot va, GXD_AlleleGenotype ag, ALL_Allele aa, GXD_Genotype g, ACC_Accession aaa, PRB_Strain s, MGI_Note n
        where r.m_Marker_key = ag._Marker_key
        and ag._Allele_key = aa._Allele_key
        and ag._Genotype_key = g._Genotype_key
        and g._Strain_key = s._Strain_key
        and ag._Genotype_key = n._Object_key
        and n._MGIType_key = 12
        and n._NoteType_key = 1016
        and ag._Genotype_key = va._Object_key
        and va._AnnotType_key = 1020
        and aa._Allele_key = aaa._Object_key
        and aaa._MGIType_key = 11
        and aaa._LogicalDB_key = 1
        ''', 'auto')
        for r in results:
                key = r['m_Marker_key']
                allelePairs = r['allelePairs'].replace('\n', '')
                value = (r['_allele_key'], r['alleleSymbol'], r['alleleID'], allelePairs, r['isWildType'], r['skipAllele'], r['strain'])
                if key not in doAlleleLookup:
                        doAlleleLookup[key] = []
                doAlleleLookup[key].append(value)
        #print(doAlleleLookup)

        # number of allele references for isWildType = 0 and skipAllele = 0
        alleleRefLookup = {}
        results = db.sql('''select a._Object_key, count(distinct a._Refs_key) as counter
                from results2 r, GXD_AlleleGenotype ag, ALL_Allele aa, MGI_Reference_Assoc a
                where r.m_Marker_key = ag._Marker_key
                and ag._Allele_key = aa._Allele_key
                and aa.isWildType = 0
                and ag._Allele_key = a._Object_key
                and a._MGIType_key = 11
                and not exists (select 1 from VOC_Annot va
                        where aa._Allele_key = va._Object_key
                        and va._AnnotType_key = 1014
                        and va._Term_key in (11025588,13289567)
                        )
                group by a._Object_key
                ''', 'auto')
        for r in results:
                key = r['_Object_key']
                value = r['counter']
                if key not in alleleRefLookup:
                        alleleRefLookup[key] = []
                alleleRefLookup[key].append(value)
        #print(alleleRefLookup)

        results = db.sql('select * from results2 order by doterm, msymbol', 'auto')
        for r in results:

                key = r['m_Marker_key']
                doID = r['doID']

                # 1:  DO term name associated with human gene
                # 2:  DO term ID associated with human gene
                fp2.write(r['doTerm'] + TAB)
                fp2.write(doID + TAB)

                # 3: NOT models
                if r['qualifierTerm'] != None:
                        fp2.write(r['qualifierTerm'])
                fp2.write(TAB)

                if key in doAlleleLookup:
                        d = doAlleleLookup[key][0]
                        alleleKey = d[0]
                        alleleSymbol = d[1]
                        alleleID = d[2]
                        allelePairs = d[3]
                        isWildType = d[4]
                        skipAllele = d[5]
                        strain = d[6]

                        # 4: Allele Pairs
                        fp2.write(allelePairs + TAB)
        
                        # 5: Strain Background
                        fp2.write(strain + TAB)

                        # 6:  Allele symbol (exclude wild type allele and allele with attribute = recombinase or rransactivator)
                        # 7:  Allele MGI ID 
                        # 8:  Total number of allele references for allele in column 6/7
                        if isWildType == 0 and skipAllele == 0:
                                fp2.write(alleleSymbol + TAB)
                                fp2.write(alleleID + TAB)
                                fp2.write(str(alleleRefLookup[alleleKey][0]) + TAB)
                        else:
                                fp2.write(TAB*3)

                        # 9: Repository ID for allele
                        if alleleID in repositoryAlleleLookup:
                                fp2.write('|'.join(repositoryAlleleLookup[alleleID]))
                        fp2.write(TAB)
                else:
                        fp2.write(TAB*6)

                # 10: Mouse genotype RR id
                if doID in rridLookup:
                        fp2.write('|'.join(rridLookup[doID]))
                fp2.write(TAB)

                # 11: Mouse gene symbol from allele 
                # 12: Mouse gene ID
                # 13: Repository ID for gene (if available, can be null)
                markerID = r['markerID']
                fp2.write(r['msymbol'] + TAB)
                fp2.write(markerID + TAB)
                if markerID in repositoryGeneLookup:
                        fp2.write('|'.join(repositoryGeneLookup[markerID]))
                fp2.write(CRT)

        reportlib.finish_nonps(fp2)

#
# main
#

initialize()
processReport1()
processReport2()

