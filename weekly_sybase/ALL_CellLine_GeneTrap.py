#!/usr/local/bin/python

'''
# ALL_CellLine_GeneTrap.py
# 
# Report:
# 	Tab-delimited file of all mutant cell lines related to gene trap alleles
#
# Usage:
#	ALL_CellLine_GeneTrap.py
#
# Used by:
#	NCBI Deanna Church
#
#
# Gene Trap Alleles, 
# 	. markers
#	. mutant cell line data (sequence)
#	. cell line information (vector, strain, etc.)
#	. imsr counts of alleles by es cell line
#	. imsr counts of markers by es cellline
#	. imsr counts of alleles by strain
#	. imsr counts of markers by strain
#
# Output Format:
#
#   1.  mclID
#   2.  vector
#   3.  mclCreator
#   4.  mclLibrary - derivationName
#   5.  parentCellLine
#   6.  parentCellLineStrain
#   7.  alleleID
#   8.  alleleSymbol
#   9.  alleleName
#   10. alleleType
#   11. markerID
#   12. markersymbol
#   13. esIMSR - list of IMSR providers holding es cell lines 
#	associated to the Gene/Allele
#   14. strainIMSR - list of IMSR providers holding strain stock 
#	associated to the Gene/Allele
#   15. id.version
#   16. tag = Tag Method
#
# kstone 12/18/2014
#	- Removed IMSR db query, replaced with Sorl csv file
#
# lec   10/27/2014
#       - TR11750/postres complient
#
# 11/18/2013	lec
#	- TR11530/bug : column 10 is printing allele-status not allele-type.
#
# 04/26/2012	lec
#	- TR11035/re-organized this report to begin from the perspective of the Gene Trap Alleles
#	- the previous version took 25-30 minutes to run 8 minutes
#	
# 04/09/2012	lec
#	- convert to using one db.sql() per sql command in preparation for postgres
#
# 12/28/2011	lec
#       - changed non-ansi-standard query to left outer join
#
# 10/06/2011 sc
#	- added to public reports
#
# 05/17/2010 marka
#	- TR9864 created
#
'''

import csv
import sys
import os
import reportlib

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db

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
        for id, facilities in facilityMap.items():
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
        for id, facilities in facilityMap.items():
                facilityMap[id] = list(set(facilities))
                facilityMap[id].sort()
        
        return facilityMap


db.useOneConnection(1)

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
fp.write('strainIMSR' + TAB)
fp.write('id.version' + TAB)
fp.write('tag' + CRT)

#
# gene trap alleles
#
db.sql('''
	select a._Allele_key, 
	       a.symbol, 
	       a.name,
	       t1.term as alleleType,
	       ac.accID as alleleID
	into #genetrap
	from ALL_Allele a, VOC_Term t1, ACC_Accession ac
	where a._Allele_Status_key in (847114, 3983021)
	and a._Allele_Type_key = t1._Term_key
	and a._Allele_Type_key = 847121
        and a._Allele_key = ac._Object_key
        and ac._LogicalDB_key = 1 
        and ac._MGIType_key = 11 
        and ac.preferred = 1
	''', None)
db.sql('create index genetrap_idx1 on #genetrap(_Allele_key)', None)
db.sql('create index genetrap_idx2 on #genetrap(alleleID)', None)

#
# gene trap/markers
#
db.sql('''
	select a._Allele_key, m.symbol as markersymbol, ma.accID as markerID
	into #markers
	from #genetrap a, ALL_Marker_Assoc ama, MRK_Marker m, ACC_Accession ma
	where a._Allele_key = ama._Allele_key
	and ama._Marker_key = m._Marker_key
        and m._Marker_key = ma._Object_key
        and ma._LogicalDB_key = 1
        and ma._MGIType_key = 2
        and ma.preferred = 1
	''', 'auto')
db.sql('create index markers_idx1 on #markers(_Allele_key)', None)
db.sql('create index markers_idx2 on #markers(markerID)', None)

#
# gene trap sequences
#
db.sql('''
	select g._Allele_key, s.accID, ss.version, t.term as tag
	into #sequences 
	from #genetrap g, 
	     SEQ_Summary_View s, 
	     SEQ_Sequence ss, 
	     SEQ_GeneTrap sg, 
	     SEQ_Allele_Assoc sa, 
	     VOC_Term t
	where g._Allele_key = sa._Allele_key
	and sa._Sequence_key = s._Object_key
	and s._LogicalDB_key = 9  
	and s._MGIType_key = 19
	and s._Object_key = ss._Sequence_key 
	and s._Object_key = sg._Sequence_key 
	and sg._TagMethod_key = t._Term_key 
	''', None)
db.sql('create index sequences_idx1 on #sequences (_Allele_key)', None)

#
# cell line info
#
db.sql('''
	select g._Allele_key, 
	       cv.cellLine as mclID,
	       cv.cellLineStrain as parentCellLineStrain,
	       cv.creator as mclCreator,
	       cv.vector, 
	       cv.parentCellLine, 
	       cv.derivationName as mclLibrary
	into #celllines
	from #genetrap g, ALL_Allele_CellLine_View cv
	where g._Allele_key = cv._Allele_key and cv.isMutant = 1  
	''', None)
db.sql('create index celllines_idx1 on #celllines (_Allele_key)', None)

# get imsr allele & marker ID => providers for ES Cells, and for Mice
es = createImsrEsDict()
strain = createImsrStrainDict()

#
# ready to print
#
results = db.sql('''
	select 
	    g._Allele_key, 
	    g.symbol, 
	    g.name,
	    g.alleleType,
	    g.alleleID,
	    m.markerID, 
	    m.markersymbol,
	    s.accID, 
	    s.version, 
	    s.tag, 
	    cv.mclID,
	    cv.parentCellLineStrain,
	    cv.mclCreator,
	    cv.vector, 
	    cv.parentCellLine, 
	    cv.mclLibrary
	from #genetrap g 
		LEFT OUTER JOIN #markers m on (g._Allele_key = m._Allele_key),
	     #sequences s, #celllines cv
	where g._Allele_key = s._Allele_key
	      and g._Allele_key = cv._Allele_key
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
        fp.write(row['mclLibrary'] + TAB)  
        fp.write(row['parentCellLine'] + TAB)
        fp.write(row['parentCellLineStrain'] + TAB)       
        fp.write(row['alleleID'] + TAB)
        fp.write(row['symbol'] + TAB)
        fp.write(row['name'] + TAB)
        fp.write(row['alleleType'] + TAB)
        fp.write(str(row['markerID']) + TAB)
        fp.write(str(row['markersymbol']) + TAB)
        fp.write(str(esString) + TAB)
        fp.write(str(strString) + TAB)
        fp.write(row['accID'] + '.' + row['version'] + TAB)
        fp.write(row['tag'] + CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
