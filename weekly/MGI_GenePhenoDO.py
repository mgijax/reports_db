#!/usr/local/bin/python

'''
#
# MGI_GenePheno.py
#
# Report:
#
#	See MGI_PhenoGenoMP.py
#
#       Phenotype/Genotype MP Annotations
#
#	Tab-delimited report where a row in the report
#	represents a Genotype annotation to a MP term
#	and all of the references (J:) used to annotate that term.
#
#	MGI_GenePheno
#
#	field 1: Allelic Composition
#	field 2: pipe-delimited list of Allele symbols
#	field 3: pipe-delimited list of Allele ids
#	field 4: Genetic Background
#	field 5: Mammalian Phenotype ID
#	field 6: PubMed ID
#	field 7: MGI Marker Accession ID (comma-delimited)
#	field 8: MGI Genotype Accession ID (comma-delimited)
#
#	MGI_Gene_DiseaseDO
#	field 8: DO ID(s) (comma-delimited)
#	where DO annotations do not include NOT
#
#	MGI_Gene_NotDiseaseDO
#	field 9: DO ID(s) (comma-delimited)
#	where DO annotations includes NOT
#
# Usage:
#       MGI_GenePheno.py
#
# History:
#
# lec	11/03/016
#	- TR12427/Disease Ontology (DO)
#	copied from MGI_GenePheno.py
#
'''
 
import sys 
import os
import string
import reportlib
import db

db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE(False)

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

fp1 = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp2 = reportlib.init('MGI_Geno_DiseaseDO', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp3 = reportlib.init('MGI_Geno_NotDiseaseDO', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# select all MP annotations
#

db.sql('''select distinct a._Object_key, a._Term_key, e._Refs_key 
	into temporary table mp 
	from VOC_Annot a, VOC_Evidence e, VOC_Term t 
	where a._AnnotType_key = 1002 
	and a._Qualifier_key = t._Term_key 
	and t.term is null 
	and a._Annot_key = e._Annot_key
	''', None)
db.sql('create index idx1 on mp(_Object_key)', None)
db.sql('create index idx2 on mp(_Term_key)', None)
db.sql('create index idx3 on mp(_Refs_key)', None)

#
# resolve MP ids
#
results = db.sql('''select distinct m._Term_key, a.accID 
	from mp m, ACC_Accession a 
	where m._Term_key = a._Object_key 
	and a._MGIType_key = 13 
	and a.preferred = 1
	''', 'auto')
mpID = {}
for r in results:
    key = r['_Term_key']
    value = r['accID']
    mpID[key] = value
    
#
# resolve References
#
results = db.sql('''select distinct m._Object_key, m._Term_key, a.accID 
	from mp m, ACC_Accession a 
	where m._Refs_key = a._Object_key 
	and a._LogicalDB_key = 29
	''', 'auto')
mpRef = {}
for r in results:
    key = `r['_Object_key']` + `r['_Term_key']`
    value = r['accID']
    if not mpRef.has_key(key):
	mpRef[key] = []
    mpRef[key].append(value)

#
# resolve Genotype Strains
# only include conditional = 0 (i.e. exclude conditional)
#
results = db.sql('''select distinct m._Object_key, s.strain 
	from mp m, GXD_Genotype g, PRB_Strain s 
	where m._Object_key = g._Genotype_key 
	and g.isConditional = 0 
	and g._Strain_key = s._Strain_key
	''', 'auto')
mpStrain = {}
for r in results:
    key = r['_Object_key']
    value = r['strain']
    mpStrain[key] = value

#
# resolve Allele Combination display
#
results = db.sql('''select distinct m._Object_key, nc.note 
	from mp m, MGI_Note n, MGI_NoteChunk nc 
	where m._Object_key = n._Object_key 
	and n._NoteType_key = 1016 
	and n._Note_key = nc._Note_key
	''', 'auto')
mpDisplay = {}
for r in results:
    key = r['_Object_key']
    value = string.replace(string.strip(r['note']), '\n', ',')
    mpDisplay[key] = value

#
# resolve Alleles
# only include markers of type 'gene' (1) or complex/cluster/region (10)
# only include single allele pairs
#
results = db.sql('''select distinct m._Object_key, a.symbol, aa.accID, a.isWildType
	from mp m, GXD_AlleleGenotype ag, ALL_Allele a, MRK_Marker mm, ACC_Accession aa
	where m._Object_key = ag._Genotype_key 
	and ag._Allele_key = a._Allele_key 
	and ag._Marker_key = mm._Marker_key 
	and mm._Marker_Type_key in (1,10)
	and ag._Allele_key = aa._Object_key
	and aa._MGIType_key = 11
	and aa._LogicalDB_key = 1
	and aa.prefixPart = 'MGI:'
	and aa.preferred = 1
        and not exists (select 1 from GXD_AllelePair a2
	        where ag._Genotype_key = a2._Genotype_key
		and a2.sequenceNum > 1)
	''', 'auto')
mpAlleles = {}
mpAlleleIDs = {}
for r in results:
    key = r['_Object_key']
    value = r['symbol']
    value2 = r['accID']

    if not mpAlleles.has_key(key):
	mpAlleles[key] = []
    mpAlleles[key].append(value)

    if r['isWildType'] == 0:
        if not mpAlleleIDs.has_key(key):
	    mpAlleleIDs[key] = []
        mpAlleleIDs[key].append(value2)

#
# resolve Marker ID
#
results = db.sql('''select distinct m._Object_key, a.accID
	from mp m, GXD_AlleleGenotype g, ACC_Accession a
	where m._Object_key = g._Genotype_key
	and g._Marker_key = a._Object_key
	and a._MGIType_key = 2
	and a._LogicalDB_key = 1
	and a.prefixPart = 'MGI:'
	and a.preferred = 1
	''', 'auto')
mpMarker = {}
for r in results:
    key = r['_Object_key']
    value = r['accID']
    if not mpMarker.has_key(key):
	mpMarker[key] = []
    mpMarker[key].append(value)

#
# resolve Genotype ID
#
results = db.sql('''select distinct m._Object_key, a.accID
	from mp m, GXD_AlleleGenotype g, ACC_Accession a
	where m._Object_key = g._Genotype_key
	and g._Genotype_key = a._Object_key
	and a._MGIType_key = 12
	and a._LogicalDB_key = 1
	and a.prefixPart = 'MGI:'
	and a.preferred = 1
	''', 'auto')
mpGenotype= {}
for r in results:
    key = r['_Object_key']
    value = r['accID']
    if not mpGenotype.has_key(key):
	mpGenotype[key] = []
    mpGenotype[key].append(value)

#
# DO annotations that do not have "NOT" qualifier
#
results = db.sql('''select distinct m._Object_key, a.accID
           from mp m, VOC_Annot va, ACC_Accession a
           where m._Object_key = va._Object_key 
           and va._AnnotType_key in (1020) 
	   and va._Qualifier_key = 1614158
	   and va._Term_key = a._Object_key
	   and a._LogicalDB_key = 191
	   and a.preferred = 1
	   ''', 'auto')
mpDO1 = {}
for r in results:
    key = r['_Object_key']
    value = r['accID']
    if not mpDO1.has_key(key):
	mpDO1[key] = []
    mpDO1[key].append(value)

#
# DO annotations that have "NOT" qualifier
#
results = db.sql('''select distinct m._Object_key, a.accID
           from mp m, VOC_Annot va, ACC_Accession a
           where m._Object_key = va._Object_key 
           and va._AnnotType_key in (1020) 
	   and va._Qualifier_key = 1614157
	   and va._Term_key = a._Object_key
	   and a._LogicalDB_key = 191
	   and a.preferred = 1
	   ''', 'auto')
mpDO2 = {}
for r in results:
    key = r['_Object_key']
    value = r['accID']
    if not mpDO2.has_key(key):
	mpDO2[key] = []
    mpDO2[key].append(value)

#
# process results
#
results = db.sql('select distinct _Object_key, _Term_key ' + \
	'from mp order by _Object_key, _Term_key', 'auto')

for r in results:

    genotype = r['_Object_key']
    term = r['_Term_key']
    refKey = `genotype` + `term`

    # we only want to list the Genotypes that have Allele Pairs
    if not mpDisplay.has_key(genotype):
       continue

    if not mpStrain.has_key(genotype):
       continue

    if not mpAlleles.has_key(genotype):
       continue

    fp1.write(mpDisplay[genotype] + TAB)
    fp1.write(string.join(mpAlleles[genotype], '|') + TAB)

    if mpAlleleIDs.has_key(genotype):
        fp1.write(string.join(mpAlleleIDs[genotype], '|'))
    fp1.write(TAB)

    fp1.write(mpStrain[genotype] + TAB)
    fp1.write(mpID[term] + TAB)
    if mpRef.has_key(refKey):
        fp1.write(string.join(mpRef[refKey], ','))
    fp1.write(TAB + string.join(mpMarker[genotype], ',') + TAB)
    fp1.write(TAB + string.join(mpGenotype[genotype], ',') + CRT)

    #
    # DO report 1
    #

    if mpDO1.has_key(genotype):

        fp2.write(mpDisplay[genotype] + TAB)
        fp2.write(string.join(mpAlleles[genotype], '|') + TAB)
        if mpAlleleIDs.has_key(genotype):
            fp2.write(string.join(mpAlleleIDs[genotype], '|'))
        fp2.write(TAB)
        fp2.write(mpStrain[genotype] + TAB)
        fp2.write(mpID[term] + TAB)
        if mpRef.has_key(refKey):
            fp2.write(string.join(mpRef[refKey], ','))
        fp2.write(TAB + string.join(mpMarker[genotype], ',') + TAB)
        fp2.write(string.join(mpDO1[genotype], ',') + CRT)

    #
    # DO report 2
    #

    if mpDO2.has_key(genotype):

        fp3.write(mpDisplay[genotype] + TAB)
        fp3.write(string.join(mpAlleles[genotype], '|') + TAB)
        if mpAlleleIDs.has_key(genotype):
            fp3.write(string.join(mpAlleleIDs[genotype], '|'))
        fp3.write(TAB)
        fp3.write(mpStrain[genotype] + TAB)
        fp3.write(mpID[term] + TAB)
        if mpRef.has_key(refKey):
            fp3.write(string.join(mpRef[refKey], ','))
        fp3.write(TAB + string.join(mpMarker[genotype], ',') + TAB)
        fp3.write(string.join(mpDO2[genotype], ',') + CRT)

reportlib.finish_nonps(fp1)	# non-postscript file
reportlib.finish_nonps(fp2)	# non-postscript file
reportlib.finish_nonps(fp3)	# non-postscript file

