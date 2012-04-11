#!/usr/local/bin/python

'''
#
# MGI_CloneSet.py
#
# Report:
#
#	MGI Clone ID
#	Clone Name
#	MGI Marker ID
#	MGI Marker Symbol
#	Clone Set
#
# Usage:
#       MGI_CloneSet.py [list of logical DB names]
#
#	1.  separate logical db names by a single comma, no spaces
#	2.  if one logical DB name is the superset, then list that one last
#
#	ex. "Image"
#	ex. "NIA 15K,NIA 7.4K,NIA"
#	ex. "RIKEN (FANTOM),RIKEN"
#
# History:
#
# lec	11/01/2004
#	- new
#
'''
 
import sys
import os
import string
import mgi_utils
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


TAB = reportlib.TAB
CRT = reportlib.CRT

#
# Main
#

l1 = string.split(sys.argv[1], ',')
lKeys = []
for l in l1:
    results = db.sql('select _LogicalDB_key from ACC_LogicalDB where name = "%s"' % (l), 'auto')
    for r in results:
        lKeys.append(r['_LogicalDB_key'])

reportName = string.split(l1[0], ' ')
fp = reportlib.init('MGI_CloneSet_' + reportName[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# get all clones

db.sql('''
	select pa._Object_key, pa._LogicalDB_key 
	into #clone1 
	from ACC_Accession pa 
	where pa._MGIType_key = 3 
	and pa._LogicalDB_key = %s
	''' % (lKeys[0]), None)
db.sql('create index idx1 on #clone1(_Object_key)', None)

for l in lKeys[1:]:

    db.sql(''''
	insert into #clone1 
	select pa._Object_key, pa._LogicalDB_key 
	from ACC_Accession pa 
	where pa._MGIType_key = 3 
	and pa._LogicalDB_key = %s
	and not exists (select 1 from #clone1 n where pa._Object_key = n._Object_key)
	''' % (l), None)

db.sql('create index clone1_idx2 on #clone1(_LogicalDB_key)', None)

# grab logical DB name

db.sql('''
	select n._Object_key, db = ldb.name 
	into #clone2 
	from #clone1 n, ACC_LogicalDB ldb 
	where n._LogicalDB_key = ldb._LogicalDB_key
	''', None)
db.sql('create index clone2_idx1 on #clone2(_Object_key)', None)

# grab all associated markers

db.sql('''
	select n._Object_key, pm._Marker_key 
	into #clone3 
	from #clone2 n, PRB_Marker pm 
	where n._Object_key = pm._Probe_key
	''', None)
db.sql('create index clone3_idx1 on #clone3(_Marker_key)', None)

# grab marker symbol, accID

results = db.sql('''
	select n._Object_key, m.symbol, ma.accID as markerID
	from #clone3 n, MRK_Marker m, ACC_Accession ma 
	where n._Marker_key = m._Marker_key  
	and n._Marker_key = ma._Object_key 
	and ma._MGIType_key = 2 
	and ma._LogicalDB_key = 1 
	and ma.prefixPart = "MGI:" 
	and ma.preferred = 1
	''', 'auto')
markers = {}
for r in results:
    key = r['_Object_key']
    value = r
    if not markers.has_key(key):
	markers[key] = []
    markers[key].append(value) 

# grab segment name, accID

results = db.sql('''
	select n._Object_key, n.db, p.name, pa.accID as segmentID
	from #clone2 n, PRB_Probe p, ACC_Accession pa 
	where n._Object_key = p._Probe_key  
	and n._Object_key = pa._Object_key  
	and pa._MGIType_key = 3 
	and pa._LogicalDB_key = 1 
	and pa.prefixPart = "MGI:" 
	and pa.preferred = 1 
	order by n.db, name
	''', 'auto')
for r in results:
    key = r['_Object_key']

    if markers.has_key(key):
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

