#!/usr/local/bin/python

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
# 10/06/2011 sc
#       - added to public reports
#
# 05/17/2010 marka
#       - TR9864 created
#
'''

import string
import sys
import db
import reportlib
import os

CRT = reportlib.CRT
TAB = reportlib.TAB

db.useOneConnection(1)

cmds=[]
sel=[]


db.sql('''select distinct c.cellLine as mclID, c.parentCellLine, 
	    c.vector, c.creator as mclCreator, c.derivationName as library,  
	    c.parentCellLineStrain, ac._Allele_key, aa.accID as alleleID, 
	    avv._Marker_key, ma.accID as markerID 
	into #gt_seqs 
	from ALL_CellLine_View c, ALL_Allele_CellLine ac, ALL_Allele_View avv, 
	ACC_Accession cc, ACC_Accession aa, ACC_Accession ma
	where c.cellLine = cc.accID 
	    and cc._MGIType_key = 28 
	    and cc._LogicalDB_key in (108, 109, 137) 
	    and c._CellLine_key = ac._MutantCellLine_key 
	    and ac._Allele_key = aa._Object_key
	    and aa._LogicalDB_key = 1 
	    and aa._MGIType_key = 11 
	    and aa.private = 0  
	    and aa.preferred = 1  
	    and ac._Allele_key = avv._Allele_key 
	    and avv._Allele_Status_key in (847114, 3983021) 
	    and avv._Allele_Type_key 
		in (847116,847117,847118,847119,847120)
	    and avv._Marker_key *= ma._Object_key
	    and ma._LogicalDB_key = 1 
	    and ma._MGIType_key = 2 
	    and ma.private = 0  
	    and ma.preferred = 1''', None)
	
	
db.sql("create index idx_mrk on #gt_seqs (markerID)", None)

db.sql("create index idx_allp on #gt_seqs (alleleID)", None)

db.sql("create index idx_gtallid on #gt_seqs (_Allele_key)", None)

db.sql("create index idx_gtmrkid on #gt_seqs (_Marker_key)", None)

cmds.append("create table #imsrCounts(accID varchar(30), abbrevName varchar(30), cType int)")

cmds.append('''insert into #imsrCounts select distinct ac.accID, f.abbrevName, 1
	from #gt_seqs g, imsr..StrainFacilityAssoc sfa, imsr..SGAAssoc sga, 
	    imsr..Accession ac, imsr..Facility f 
	where ac.accID = g.alleleID 
	and ac._IMSRType_key = 3 
	and ac._Object_key = sga._Allele_key 
	and sga._Strain_key = sfa._Strain_key
	and sfa._StrainState_key = 2 
	and sfa._Facility_key = f._Facility_key''')

cmds.append('''insert into #imsrCounts select distinct ac.accID, f.abbrevName, 1 
	from #gt_seqs g, imsr..StrainFacilityAssoc sfa, imsr..SGAAssoc sga, 
	    imsr..Accession ac, imsr..Facility f 
	where ac.accID = g.markerID 
	and ac._IMSRType_key = 2 
	and ac._Object_key = sga._Gene_key 
	and sga._Strain_key = sfa._Strain_key
	and sfa._StrainState_key = 2
	and sfa._Facility_key = f._Facility_key''')

cmds.append('''insert into #imsrCounts select distinct ac.accID, f.abbrevName, 2
	from #gt_seqs g, imsr..StrainFacilityAssoc sfa, imsr..SGAAssoc sga,  
	    imsr..Accession ac, imsr..Facility f 
	where ac.accID = g.alleleID 
	and ac._IMSRType_key = 3 
	and ac._Object_key = sga._Allele_key 
	and sga._Strain_key = sfa._Strain_key
	and sfa._StrainState_key <> 2 
	and sfa._Facility_key = f._Facility_key''')

cmds.append('''insert into #imsrCounts select distinct ac.accID, f.abbrevName, 2
	from #gt_seqs g, imsr..StrainFacilityAssoc sfa, imsr..SGAAssoc sga,
	    imsr..Accession ac, imsr..Facility f 
	where ac.accID = g.markerID 
	    and ac._IMSRType_key = 2 
	    and ac._Object_key = sga._Gene_key 
	    and sga._Strain_key = sfa._Strain_key 
	    and sfa._StrainState_key <> 2 
	    and sfa._Facility_key = f._Facility_key''')

sel.append("select distinct accID, abbrevName from #imsrCounts where cType = 1")

sel.append("select distinct accID, abbrevName from #imsrCounts where cType = 2")
	
sel.append('''select g.mclID, g.vector, g.mclCreator, g.library, 
	g.parentCellLine, g.parentCellLineStrain,  
	    g.alleleID, al.symbol as alleleSymbol, 
	    al.name as alleleName, t.term as alleleType,  
	    mrk.symbol as markerSymbol, g.markerID 
	from #gt_seqs g, ALL_Allele_View al, MRK_Marker mrk, VOC_Term t
	where g._Allele_key = al._Allele_key 
	    and al._Allele_Type_key = t._Term_key 
	    and g._Marker_key = mrk._Marker_key''')

db.sql(cmds, 'auto')

# get results
results = db.sql(sel, 'auto')

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
		
es = {}
for c in results [-3]:
	accID = c['accID']
	if accID not in es:
		es[accID] = [c['abbrevName']]
	else:
		es[accID].append(c['abbrevName'])
	
strain = {}
for c in results [-2]:
	accID = c['accID']
	if accID not in strain:
		strain[accID] = [c['abbrevName']]
	else:
		strain[accID].append(c['abbrevName'])

for row in results [-1]:
	
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
db.useOneConnection(0)
