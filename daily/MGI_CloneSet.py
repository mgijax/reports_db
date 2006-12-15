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
import db
import mgi_utils
import reportlib

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

print 'query 1 begin...%s' % (mgi_utils.date())
db.sql('select pa._Object_key, pa._LogicalDB_key ' + \
	'into #clone1 ' + \
	'from ACC_Accession pa ' + \
	'where pa._MGIType_key = 3 ' + \
	'and pa._LogicalDB_key = %s' % (lKeys[0]), None)
db.sql('create index idx1 on #clone1(_Object_key)', None)

for l in lKeys[1:]:

    print 'query 1(%s) begin...%s' % (l, mgi_utils.date())
    cmd = 'insert into #clone1 ' + \
	    'select pa._Object_key, pa._LogicalDB_key ' + \
	    'from ACC_Accession pa ' + \
	    'where pa._MGIType_key = 3 ' + \
	    'and pa._LogicalDB_key = %s ' % (l) + \
	    'and not exists (select 1 from #clone1 n where pa._Object_key = n._Object_key)'
    db.sql(cmd, None)
    print 'query 1(%s) end...%s' % (l, mgi_utils.date())

db.sql('create index idx2 on #clone1(_LogicalDB_key)', None)
print 'query 1 end...%s' % (mgi_utils.date())

# grab logical DB name

print 'query 2 begin...%s' % (mgi_utils.date())
db.sql('select n._Object_key, db = ldb.name ' + \
	'into #clone2 ' + \
	'from #clone1 n, ACC_LogicalDB ldb ' + \
	'where n._LogicalDB_key = ldb._LogicalDB_key', None)
db.sql('create index idx1 on #clone2(_Object_key)', None)
print 'query 2 end...%s' % (mgi_utils.date())

# grab all associated markers

print 'query 3 begin...%s' % (mgi_utils.date())
db.sql('select n._Object_key, pm._Marker_key ' + \
	'into #clone3 ' + \
	'from #clone2 n, PRB_Marker pm ' + \
	'where n._Object_key = pm._Probe_key', None)
db.sql('create index idx1 on #clone3(_Marker_key)', None)
print 'query 3 end...%s' % (mgi_utils.date())

# grab marker symbol, accID

print 'query 4 begin...%s' % (mgi_utils.date())
cmd = 'select n._Object_key, markerSymbol = m.symbol, markerID = ma.accID ' + \
	'from #clone3 n, MRK_Marker m, ACC_Accession ma ' + \
	'where n._Marker_key = m._Marker_key  ' + \
	'and n._Marker_key = ma._Object_key ' + \
	'and ma._MGIType_key = 2 ' + \
	'and ma._LogicalDB_key = 1 ' + \
	'and ma.prefixPart = "MGI:" ' + \
	'and ma.preferred = 1'
results = db.sql(cmd, 'auto')
print 'query 4 end...%s' % (mgi_utils.date())
markers = {}
for r in results:
    key = r['_Object_key']
    value = r
    if not markers.has_key(key):
	markers[key] = []
    markers[key].append(value) 

# grab segment name, accID

print 'query 5 begin...%s' % (mgi_utils.date())
cmd = 'select n._Object_key, n.db, segmentName = p.name, segmentID = pa.accID ' + \
	'from #clone2 n, PRB_Probe p, ACC_Accession pa ' + \
	'where n._Object_key = p._Probe_key  ' + \
	'and n._Object_key = pa._Object_key  ' + \
	'and pa._MGIType_key = 3 ' + \
	'and pa._LogicalDB_key = 1 ' + \
	'and pa.prefixPart = "MGI:" ' + \
	'and pa.preferred = 1 ' + \
	'order by n.db, segmentName'
results = db.sql(cmd, 'auto')
print 'query 5 end...%s' % (mgi_utils.date())

for r in results:
    key = r['_Object_key']

    if markers.has_key(key):
	allMarkers = markers[key]
	for m in allMarkers:
            fp.write(r['segmentID'] + TAB + \
	             r['segmentName'] + TAB + \
	             m['markerID'] + TAB + \
	             m['markerSymbol'] + TAB + \
		     r['db'] + CRT)
    else:
        fp.write(r['segmentID'] + TAB + \
	         r['segmentName'] + TAB + \
	         TAB + \
	         TAB + \
		 r['db'] + CRT)

reportlib.finish_nonps(fp)

