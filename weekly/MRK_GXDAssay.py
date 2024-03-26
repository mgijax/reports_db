
'''
#
# MRK_GXDAssay.py 08/28/2002
#
# Report:
#       Tab-delimited file of MGI Genes with associated GXD Assay records
#
# Usage:
#       MRK_GXDAssay.py
#
# Used by:
#	Deanna Church of NCBI
#
# History:
#
# lec	03/27/2009
#	- TR 9568; do not include recombinase assay (11)
#
# lec	04/23/2008
#	- TR 8775; do not include new assay types
#
# lec	08/28/2002
#	- TR 4027
#
'''
 
import sys
import os
import reportlib
import db

db.setTrace()

TAB = reportlib.TAB
CRT = reportlib.CRT

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

results = db.sql('''
        select distinct a.accID, g._Marker_key 
        from GXD_Assay g, ACC_Accession a 
        where g._AssayType_key not in (10,11) 
        and g._Assay_key = a._Object_key 
        and a._MGIType_key = 8 
        and a.prefixPart = 'MGI:' 
        and a._LogicalDB_key = 1 
        and a.preferred = 1
        ''', 'auto')
ids = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if key not in ids:
        ids[key] = []
    ids[key].append(value)

results = db.sql('''
        select distinct a.accID, m._Marker_key, m.symbol 
        from MRK_Marker m, GXD_Assay g, ACC_Accession a 
        where m._Marker_key = g._Marker_key 
        and g._AssayType_key not in (10,11) 
        and m._Marker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a.prefixPart = 'MGI:' 
        and a._LogicalDB_key = 1 
        and a.preferred = 1
        ''', 'auto')

for r in results:
        fp.write(r['accID'] + TAB + r['symbol'] + TAB + ','.join(ids[r['_Marker_key']]) + CRT)

reportlib.finish_nonps(fp)
