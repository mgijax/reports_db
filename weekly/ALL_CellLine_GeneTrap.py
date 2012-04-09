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
#	associated to the Gene/Allele
#   14. strainIMSR - list of IMSR providers holding strain stock 
#	associated to the Gene/Allele
#   15. id.version
#   16. tag = Tag Method
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

import sys
import os
import db
import reportlib

CRT = reportlib.CRT
TAB = reportlib.TAB

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

db.sql('''select distinct s._Object_key, s.accID, ss.version, 
	    ac.accID as alleleID, 
	    a.symbol, a.name, m.accID as markerID, 
	    a.markerSymbol as markerSymbol,
	    allt.term as alleleType, a._Allele_key, t.term 
	into #gt_seqs 
	from SEQ_Summary_View s, SEQ_Sequence ss, SEQ_GeneTrap sg, VOC_Term t, 
	    SEQ_Allele_Assoc sa
                LEFT OUTER JOIN ALL_Allele_View a on (
                            sa._Allele_key = a._Allele_key
                            and a._Allele_Status_key in (847114, 3983021))
                LEFT OUTER JOIN ACC_Accession ac on (
                            a._Allele_key = ac._Object_key
	                    and ac._LogicalDB_key = 1 
	                    and ac._MGIType_key = 11 
	                    and ac.private = 0 
	                    and ac.preferred = 1
                            )
                LEFT OUTER JOIN VOC_Term allt on (a._Allele_Type_key = allt._Term_key)
                LEFT OUTER JOIN ACC_Accession m on (
                            a._Marker_key = m._Object_key
                            and m._LogicalDB_key = 1
                            and m._MGIType_key = 2
                            and m.private = 0
                            and m.preferred = 1
                            )
	where s._LogicalDB_key = 9  
	    and s._MGIType_key = 19
	    and s._Object_key = ss._Sequence_key 
	    and s._Object_key = sg._Sequence_key 
	    and sg._TagMethod_key = t._Term_key 
	    and s._Object_key = sa._Sequence_key
	    ''', None)

db.sql("create index idx_gtallelekey on #gt_seqs (_Allele_key)", None)
db.sql("create index idx_gtallele on #gt_seqs (alleleID)", None)
db.sql("create index idx_gtmarker on #gt_seqs (markerID)", None)

db.sql('''select g._Allele_key, cv.vector, cv.parentCellLine, s2.strain 
	into #cells
	from #gt_seqs g, ALL_Allele_CellLine_View cv, ALL_CellLine p, 
	    PRB_Strain s2 
	where g._Allele_key = cv._Allele_key  
	    and cv.isMutant = 1  
	    and cv.parentCellLine_key = p._CellLine_key 
	    and p._Strain_key = s2._Strain_key''', None)
db.sql('''create index idx_alleleKey on #cells (_Allele_key)''', None)
db.sql('''create table #imsrCounts(accID varchar(30), abbrevName varchar(30), cType int)''', None)
	
db.sql('''insert into #imsrCounts select distinct ac.accID, f.abbrevName, 1
	from #gt_seqs g, imsr..StrainFacilityAssoc sfa, imsr..SGAAssoc sga, 
	    imsr..Accession ac, imsr..Facility f 
	where ac.accID = g.alleleID 
	and ac._IMSRType_key = 3 
	and ac._Object_key = sga._Allele_key 
	and sga._Strain_key = sfa._Strain_key
	and sfa._StrainState_key = 2 
	and sfa._Facility_key = f._Facility_key''', None)

db.sql('''insert into #imsrCounts select distinct ac.accID, f.abbrevName, 1 
	from #gt_seqs g, imsr..StrainFacilityAssoc sfa, imsr..SGAAssoc sga, 
	    imsr..Accession ac, imsr..Facility f 
	where ac.accID = g.markerID 
	and ac._IMSRType_key = 2 
	and ac._Object_key = sga._Gene_key 
	and sga._Strain_key = sfa._Strain_key
	and sfa._StrainState_key = 2
	and sfa._Facility_key = f._Facility_key''', None)

db.sql('''insert into #imsrCounts select distinct ac.accID, f.abbrevName, 2
	from #gt_seqs g, imsr..StrainFacilityAssoc sfa, imsr..SGAAssoc sga,  
	    imsr..Accession ac, imsr..Facility f 
	where ac.accID = g.alleleID 
	and ac._IMSRType_key = 3 
	and ac._Object_key = sga._Allele_key 
	and sga._Strain_key = sfa._Strain_key
	and sfa._StrainState_key <> 2 
	and sfa._Facility_key = f._Facility_key''', None)

db.sql('''insert into #imsrCounts select distinct ac.accID, f.abbrevName, 2
	from #gt_seqs g, imsr..StrainFacilityAssoc sfa, imsr..SGAAssoc sga,
	    imsr..Accession ac, imsr..Facility f 
	where ac.accID = g.markerID 
	    and ac._IMSRType_key = 2 
	    and ac._Object_key = sga._Gene_key 
	    and sga._Strain_key = sfa._Strain_key 
	    and sfa._StrainState_key <> 2 
	    and sfa._Facility_key = f._Facility_key''', None)
	
results = db.sql('''select distinct accID, abbrevName 
	from #imsrCounts 
	where cType = 1''', 'auto')
es = {}
for c in results:
	accID = c['accID']
	if accID not in es:
		es[accID] = [c['abbrevName']]
	else:
		es[accID].append(c['abbrevName'])
	
results = db.sql('''select distinct accID, abbrevName 
	from #imsrCounts 
	where cType = 2''', 'auto')
strain = {}
for c in results:
	accID = c['accID']
	if accID not in strain:
		strain[accID] = [c['abbrevName']]
	else:
		strain[accID].append(c['abbrevName'])
	
results = db.sql('''select distinct g.accID, g.version, g.term as tag, 
	    cv.cellLine as mclID, cv.derivationName as mclLibrary, 
	    c.parentCellLine, c.strain as parentCellLineStrain, 
	    c.vector, cv.creator as mclCreator, g.alleleID, g.symbol, 
	    g.name, g.alleleType, g.markerID, g.markerSymbol 
	from #gt_seqs g, #cells c, ALL_Allele_CellLine ac, 
	    ALL_Allele_CellLine_View cv 
	where g._Allele_key = ac._Allele_key 
	    and g._Allele_key = c._Allele_key 
	    and ac._MutantCellLine_key = cv._MutantCellLine_key''', 'auto')

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
        fp.write(str(row['markerSymbol']) + TAB)
        fp.write(str(esString) + TAB)
        fp.write(str(strString) + TAB)
        fp.write(row['accID'] + '.' + row['version'] + TAB)
        fp.write(row['tag'] + CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
