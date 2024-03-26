
'''
#
# NorCOMM_Allele.py
#
# Reports:
#
# TR 10549
#
# Generate a tab delimited public report containing the list of 
# official nomenclature for all NorCOMM alleles in MGI db 
# (this includes ES Cells and mice being generated by NorCOMM).
#
# The NorCOMM alleles can be identified by their mutant cell lines having logical 
# database NorCOMM
#
# Format of the report should be one line per NorCOMM allele w/ fields:
#
#   mutant cell line ID (from one of the logical dbs above)
#   logical db name from above (to differentiate CSD from Regeneron)
#   MGI id of the allele
#   allele symbol
#   allele name
#   MGI id of the marker for the allele
#   marker symbol
#
# Usage:
#       NorCOMM_Allele.py
#
# History:
#
# sc	02/01/2011
#	- created, TR 10549
#
'''
 
import sys 
import os
import reportlib
import db

db.setTrace()

CRT = reportlib.CRT
TAB = reportlib.TAB

#
# Main
#

fp = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])
fp.write('#\n')
fp.write('# This report lists the official nomenclature for all NorCOMM alleles in MGI\n')
fp.write('#\n')
fp.write('#  column 1: project ID\n')
fp.write('#  column 2: logical db name\n')
fp.write('#  column 3: MGI id of the allele\n')
fp.write('#  column 4: allele symbol\n')
fp.write('#  column 5: allele name\n')
fp.write('#  column 6: MGI id of the marker\n')
fp.write('#  column 7: marker symbol\n')
fp.write('#  column 8: mutant cell line IDs\n')
fp.write('#\n\n')

#
# select all NorCOMM alleles
#

db.sql('''select a._Allele_key, a._Marker_key, a.symbol, a.name,
        aa.accID as projectID, ldb.name as ldb, m.symbol as markerSymbol
        into temporary table norcomm
        from ALL_Allele a, ACC_Accession aa, ACC_LogicalDB ldb, MRK_Marker m 
        where a._Allele_key = aa._Object_key 
        and aa._MGIType_key = 11 
        and aa._LogicalDB_key = 143
        and aa._LogicalDB_key = ldb._LogicalDB_key 
        and a._Marker_key = m._Marker_key''', None)

db.sql('create index idx1 on norcomm(_Allele_key)', None)
db.sql('create index idx2 on norcomm(_Marker_key)', None)

#
# allele ids
#

results = db.sql('''select k._Allele_key, a.accID 
      from norcomm k, ACC_Accession a 
      where k._Allele_key = a._Object_key 
      and a._MGIType_key = 11 
      and a._LogicalDB_key = 1 
      and a.prefixPart = 'MGI:' 
      and a.preferred = 1''', 'auto')
alleleID = {}
for r in results:
    key = r['_Allele_key']
    value = r['accID']
    alleleID[key] = value

# marker ids

results = db.sql('''select distinct k._Marker_key, a.accID
      from norcomm k, ACC_Accession a 
      where k._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 1 
      and a.prefixPart = 'MGI:' 
      and a.preferred = 1''', 'auto')
markerID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    markerID[key] = value

#
# mutant cell line ids
#

results = db.sql('''select distinct k._Allele_key, a.accID
        from norcomm k, ALL_Allele_CellLine c, ACC_Accession a
        where k._Allele_key = c._Allele_key
        and c._MutantCellLine_key = a._Object_key
        and a._MGIType_key = 28
        and a._LogicalDB_key = 142''', 'auto')
mutantID = {}
for r in results:
    key = r['_Allele_key']
    value = r['accID']
    if key not in mutantID:
        mutantID[key] = []
    mutantID[key].append(value)

#
# process results
#
results = db.sql('select * from norcomm order by projectID', 'auto')

for r in results:

    fp.write(r['projectID'] + TAB + \
             r['ldb'] + TAB + \
             alleleID[r['_Allele_key']] + TAB + \
             r['symbol'] + TAB + \
             r['name'] + TAB + \
             markerID[r['_Marker_key']] + TAB + \
             r['markerSymbol'] + TAB)

    if r['_Allele_key'] in mutantID:
        fp.write(','.join(mutantID[r['_Allele_key']]))
    fp.write(CRT)

reportlib.finish_nonps(fp)	# non-postscript file
