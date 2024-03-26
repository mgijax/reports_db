
'''
# ALL_CellLine_Targeted.py
#
# Report:
#       Tab-delimited file of all mutant cell lines related to targeted alleles
#
# Usage:
#       ALL_CellLine_Targeted.py
#
# Used by:
#       NCBI Deanna Church
#
# Output Format:
#
#   1. mclID
#   2. vector
#   3. mclCreator
#   4. mclLibrary - derivationName
#   5. parentCellLine
#   6. parentCellLineStrain
#   7. alleleID
#   8. alleleSymbol
#   9. alleleName
#   10. alleleType
#   11. markerID
#   12. markerSymbol
#   13. esIMSR - list of IMSR providers holding es cell lines
#       associated to the Gene/Allele
#   14. strainIMSR - list of IMSR providers holding strain stock
#       associated to the Gene/Allele
#
# kstone 12/18/2014
#	- Removed IMSR query, replaced with Solr csv file
#
# lec   10/27/2014
#       - TR11750/postres complient
#
# 12/31/2013	lec
#	- TR11515/changes to _Allele_Type_key
#
# 04/09/2012    lec
#       - convert to using one db.sql() per sql command in preparation for postgres
#
# 12/28/2011	lec
#	- changed non-ansi-standard query to left outer join
#
# 10/06/2011 sc
#       - added to public reports
#
# 05/17/2010 marka
#       - TR9864 created
#
'''

import csv
import sys
import os
import reportlib
import db
import db

db.setTrace()

CRT = reportlib.CRT
TAB = reportlib.TAB

IMSR_CSV = os.environ['IMSR_STRAINS_CSV']

# functions

def createImsrEsDict():
        """
        Reads in the IMSR allStrains.csv file and generates
        a map of allele/marker ID to facility abbreviations
        for es cells only
        """
        facilityMap = {}
        csvfile = open(IMSR_CSV, 'r')
        reader = csv.reader(csvfile)
        for row in reader:
                allele_ids = row[0]
                marker_ids = row[1]
                provider = row[3]
                strain_states = row[4]
                if strain_states:
                        for strain_state in strain_states.split(','):
                                if strain_state.lower() != 'es cell':
                                        continue
                                if allele_ids:
                                        for id in allele_ids.split(','):
                                                facilityMap.setdefault(id, []).append(provider)
                                if marker_ids:
                                        for id in marker_ids.split(','):
                                                facilityMap.setdefault(id, []).append(provider)

        # unique and sort the provider list
        for id, facilities in list(facilityMap.items()):
                facilityMap[id] = list(set(facilities))
                facilityMap[id].sort()

        return facilityMap

def createImsrStrainDict():
        """
        Reads in the IMSR allStrains.csv file and generates
        a map of allele/marker ID to facility abbreviations
        for any type except es cells
        """
        facilityMap = {}
        csvfile = open(IMSR_CSV, 'r')
        reader = csv.reader(csvfile)
        for row in reader:
                allele_ids = row[0]
                marker_ids = row[1]
                provider = row[3]
                strain_states = row[4]
                if strain_states:
                        for strain_state in strain_states.split(','):
                                if strain_state.lower() == 'es cell':
                                        continue
                                if allele_ids:
                                        for id in allele_ids.split(','):
                                                facilityMap.setdefault(id, []).append(provider)
                                if marker_ids:
                                        for id in marker_ids.split(','):
                                                facilityMap.setdefault(id, []).append(provider)

        # unique and sort the provider list
        for id, facilities in list(facilityMap.items()):
                facilityMap[id] = list(set(facilities))
                facilityMap[id].sort()

        return facilityMap

#create report
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('#mclID' + TAB)
fp.write('vector' + TAB)
fp.write('mclCreator' + TAB)
fp.write('mclLibrary' + TAB)
fp.write('parentCellLine' + TAB)
fp.write('parentCellLineStrain' + TAB)
fp.write('alleleID' + TAB)
fp.write('alleleSymbol' + TAB)
fp.write('alleleName' + TAB)
fp.write('alleleType' + TAB)
fp.write('markerID' + TAB)
fp.write('markerSymbol' + TAB)
fp.write('esIMSR' + TAB)
fp.write('strainIMSR' + CRT)
                
#
# _Allele_Type_key = 847116 => 'Targeted'
#

db.sql('''
        select distinct c.cellLine as mclID, c.parentCellLine, 
            c.vector, c.creator as mclCreator, c.derivationName as library,  
            c.parentCellLineStrain, ac._Allele_key, aa.accID as alleleID, 
            avv._Marker_key, ma.accID as markerID 
        into temporary table gt_seqs 
        from ALL_CellLine_View c, ALL_Allele_CellLine ac, 
                ACC_Accession cc, ACC_Accession aa, ACC_Accession ma,
                ALL_Allele avv
        where c.cellLine = cc.accID 
            and cc._MGIType_key = 28 
            and cc._LogicalDB_key in (108, 109, 137) 
            and c._CellLine_key = ac._MutantCellLine_key 
            and ac._Allele_key = avv._Allele_key 
            and avv._Allele_Status_key in (847114, 3983021) 
            and avv._Allele_Type_key = 847116
            and ac._Allele_key = aa._Object_key
            and aa._LogicalDB_key = 1 
            and aa._MGIType_key = 11 
            and aa.private = 0  
            and aa.preferred = 1  
            and avv._Marker_key = ma._Object_key
            and ma._LogicalDB_key = 1 
            and ma._MGIType_key = 2 
            and ma.private = 0  
            and ma.preferred = 1
            ''', None)
        
db.sql('create index idx_mrk on gt_seqs (markerID)', None)
db.sql('create index idx_allp on gt_seqs (alleleID)', None)
db.sql('create index idx_gtallid on gt_seqs (_Allele_key)', None)
db.sql('create index idx_gtmrkid on gt_seqs (_Marker_key)', None)

# get imsr allele & marker ID => providers for ES Cells, and for Mice
es = createImsrEsDict()
strain = createImsrStrainDict()

results = db.sql('''
        select g.mclID, g.vector, g.mclCreator, g.library, 
        g.parentCellLine, g.parentCellLineStrain,  
            g.alleleID, al.symbol as alleleSymbol, 
            al.name as alleleName, t.term as alleleType,  
            mrk.symbol as markerSymbol, g.markerID 
        from gt_seqs g, ALL_Allele_View al, MRK_Marker mrk, VOC_Term t
        where g._Allele_key = al._Allele_key 
            and al._Allele_Type_key = t._Term_key 
            and g._Marker_key = mrk._Marker_key
            ''', 'auto')
for row in results:
        
        esProviders = set()
        strainProviders = set()
        
        if row['alleleID'] in es:
            esProviders = esProviders.union(set(es[row['alleleID']]))
        if row['markerID'] in es:
            esProviders = esProviders.union(set(es[row['markerID']]))
        if row['alleleID'] in strain:
            strainProviders = strainProviders.union(set(strain[row['alleleID']]))
        if row['markerID'] in strain:
            strainProviders = strainProviders.union(set(strain[row['markerID']]))
        
        esString =  ", ".join(esProviders)
        strString =  ", ".join(strainProviders)
        
        fp.write(row['mclID'] + TAB)
        fp.write(row['vector'] + TAB)
        fp.write(row['mclCreator'] + TAB)
        fp.write(row['library'] + TAB)
        fp.write(row['parentCellLine'] + TAB)
        fp.write(row['parentCellLineStrain'] + TAB)
        fp.write(row['alleleID'] + TAB)
        fp.write(row['alleleSymbol'] + TAB)
        fp.write(row['alleleName'] + TAB)
        fp.write(row['alleleType'] + TAB)
        fp.write(str(row['markerID']) + TAB)
        fp.write(str(row['markerSymbol']) + TAB)
        fp.write(str(esString) + TAB)
        fp.write(str(strString) + CRT)

reportlib.finish_nonps(fp)
