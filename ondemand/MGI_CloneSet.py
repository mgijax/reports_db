'''
#
# MGI_CloneSet.py
#
# Report:
#
#       MGI Clone ID
#       Clone Name
#       MGI Marker ID
#       MGI Marker Symbol
#       Clone Set
#
# Usage:
#       MGI_CloneSet.py [list of logical DB names]
#
#       1.  separate logical db names by a single comma, no spaces
#       2.  if one logical DB name is the superset, then list that one last
#
#       ex. "Image"
#       ex. "NIA 15K,NIA 7.4K,NIA"
#       ex. "RIKEN (FANTOM),RIKEN"
#
# History:
#
# sc    03/21/20 python 3 upgrade
#
# lec   11/01/2004
#       - new
#
'''
 
import sys
import os
import string
import mgi_utils
import reportlib
import db

db.setTrace()

TAB = reportlib.TAB
CRT = reportlib.CRT

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

db.sql('''
        select pa._Object_key, pa._LogicalDB_key, ldb.name as db
        into temporary table clone1 
        from ACC_Accession pa, ACC_LogicalDB ldb
        where pa._MGIType_key = 3 
        and pa._LogicalDB_key = ldb._LogicalDB_key
        and ldb.name in ('IMAGE', 'NIA 15K', 'NIA 7.4K','NIA', 'RIKEN (FANTOM)', 'RIKEN', 'RPCI-23', 'RPCI-24')
        ''', None)
db.sql('create index clone1_idx1 on clone1(_Object_key)', None)
db.sql('create index clone1_idx2 on clone1(_LogicalDB_key)', None)

# grab all associated markers

db.sql('''
        select n._Object_key, n.db, pm._Marker_key 
        into temporary table clone2 
        from clone1 n, PRB_Marker pm 
        where n._Object_key = pm._Probe_key
        ''', None)
db.sql('create index clone2_idx1 on clone2(_Marker_key)', None)

# grab marker symbol, accID

results = db.sql('''
        select n._Object_key, m.symbol, ma.accID as markerID
        from clone2 n, MRK_Marker m, ACC_Accession ma 
        where n._Marker_key = m._Marker_key  
        and n._Marker_key = ma._Object_key 
        and ma._MGIType_key = 2 
        and ma._LogicalDB_key = 1 
        and ma.prefixPart = 'MGI:' 
        and ma.preferred = 1
        ''', 'auto')
markers = {}
for r in results:
    key = r['_Object_key']
    value = r
    if key not in markers:
        markers[key] = []
    markers[key].append(value) 

# grab segment name, accID

results = db.sql('''
        select n._Object_key, n.db, p.name, pa.accID as segmentID
        from clone2 n, PRB_Probe p, ACC_Accession pa 
        where n._Object_key = p._Probe_key  
        and n._Object_key = pa._Object_key  
        and pa._MGIType_key = 3 
        and pa._LogicalDB_key = 1 
        and pa.prefixPart = 'MGI:' 
        and pa.preferred = 1 
        order by n.db, name
        ''', 'auto')

for r in results:
    key = r['_Object_key']

    if key in markers:
        allMarkers = markers[key]
        for m in allMarkers:
            fp.write(r['segmentID'] + TAB + \
                     r['name'] + TAB + \
                     m['markerID'] + TAB + \
                     m['symbol'] + TAB + \
                     r['db'] + CRT)
    else:
        fp.write(r['segmentID'] + TAB + \
                 r['name'] + TAB + \
                 TAB + \
                 TAB + \
                 r['db'] + CRT)

reportlib.finish_nonps(fp)

