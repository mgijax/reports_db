#!/usr/local/bin/python

'''
#
# PRB_ESTMapped.py 06/29/2000
# TR 1734
#
# Report:
#       Tab-delimited file
#       Mapped ESTs
#
# Usage:
#       PRB_ESTMapped.py
#
# Used by:
#       Incyte
#
# Notes:
#
# History:
#
# lec	11/04/2004
#	- include NIA
#
# lec	06/29/2000
#	- created
#
'''
 
import sys
import os
import db
import reportlib
import mgi_utils

#
# Main
#

db.useOneConnection(1)

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

# Image

cmds = []
cmds.append('select _Object_key into #ests from ACC_Accession where _MGIType_key = 3 and _LogicalDB_key = 17')
cmds.append('create index idx1 on #ests(_Object_key)')
db.sql(cmds, None)
print 'query 1 end...%s' % (mgi_utils.date())

# NIA

#if we need to add more....for l in [49, ...]:
db.sql('insert into #ests ' + \
    'select a._Object_key from ACC_Accession a ' + \
    'where a._MGIType_key = 3  ' + \
    'and a._LogicalDB_key = 49 ' + \
    'and not exists (select 1 from #ests e where e._Object_key = a._Object_key)', None)
print 'query 2 end...%s' % (mgi_utils.date())

# mapped markers 

cmds = []
cmds.append('select m._Marker_key, m.symbol, m.name, m.chromosome, o.offset ' + \
	'into #mappedmarkers ' + \
	'from MRK_Marker m, MRK_Offset o ' + \
	'where m._Organism_key = 1 ' + \
	'and m._Marker_Status_key in (1, 3) ' + \
	'and m.chromosome != "UN" ' + \
	'and m._Marker_key = o._Marker_key ' + \
	'and o.source = 0 ' + \
	'order by m.symbol')
cmds.append('create index idx1 on #mappedmarkers(_Marker_key)')
db.sql(cmds, None)
print 'query 3 end...%s' % (mgi_utils.date())

# mapped markers for ESTs

cmds = []
cmds.append('select e._Object_key, m.* ' + \
	'into #ests2 ' + \
	'from  #ests e, PRB_Marker pm, #mappedmarkers m ' + \
	'where e._Object_key = pm._Probe_key ' + \
	'and pm._Marker_key = m._Marker_key ')
cmds.append('create index idx1 on #ests2(_Object_key)')
cmds.append('create index idx2 on #ests2(symbol)')
db.sql(cmds, None)
print 'query 4 end...%s' % (mgi_utils.date())

# MGI Acc ID for ESTs

results = db.sql('select e._Object_key, a.accID ' + \
      'from #ests2 e, ACC_Accession a ' + \
      'where e._Object_key = a._Object_key ' + \
      'and a._MGIType_key = 3 ' + \
      'and a.prefixPart = "MGI:"' + \
      'and a._LogicalDB_key = 1 ' + \
      'and a.preferred = 1', 'auto')
print 'query 5 end...%s' % (mgi_utils.date())

mgiIDs = {}
for r in results:
    key = r['_Object_key']
    value = r['accID']
    mgiIDs[key] = value

# GenBank ID for ESTs

results = db.sql('select e._Object_key, a.accID ' + \
      'from #ests2 e, ACC_Accession a ' + \
      'where e._Object_key = a._Object_key ' + \
      'and a._MGIType_key = 3 ' + \
      'and a._LogicalDB_key = 9', 'auto')
print 'query 6 end...%s' % (mgi_utils.date())

genbankIDs = {}
for r in results:
    key = r['_Object_key']
    value = r['accID']
    if not genbankIDs.has_key(key):
        genbankIDs[key] = []
    genbankIDs[key].append(value)

results = db.sql('select _Object_key, symbol, name, chromosome, offset from #ests2 order by symbol', 'auto')
print 'query 7 end...%s' % (mgi_utils.date())

for r in results:

	if r['offset'] == -1.0:
	    offset = 'syntenic'
	elif r['offset'] == -999.0:
	    offset = 'N/A'
	else:
	    offset = str(r['offset'])

	if genbankIDs.has_key(r['_Object_key']):
	    for g in genbankIDs[r['_Object_key']]:
	        fp.write(mgiIDs[r['_Object_key']] + reportlib.TAB)
	        fp.write(g + reportlib.TAB)
	        fp.write(r['chromosome'] + reportlib.TAB)
	        fp.write(offset + reportlib.TAB)
	        fp.write(r['symbol'] + reportlib.TAB)
	        fp.write(r['name'] + reportlib.CRT)
        else:
	    fp.write(mgiIDs[r['_Object_key']] + reportlib.TAB)
	    fp.write(reportlib.TAB)
	    fp.write(r['chromosome'] + reportlib.TAB)
	    fp.write(offset + reportlib.TAB)
	    fp.write(r['symbol'] + reportlib.TAB)
	    fp.write(r['name'] + reportlib.CRT)

reportlib.finish_nonps(fp)

db.useOneConnection(0)
