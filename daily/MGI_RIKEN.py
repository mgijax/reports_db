#!/usr/local/bin/python

'''
#
# MGI_RIKEN.py 10/27/2003
#
# Report:
#       Tab-delimited file for RIKEN
#
#	. RIKEN Clone ID
#	. MGI Accession ID for Molecular Segment
#	. "Problem" Sequence
#	. MGI Accession ID for Marker
#	. Gene Symbol
#	. Gene Name
#
# Usage:
#       MGI_RIKEN.py
#
# Used by:
#
#	RIKEN
#
# History:
#
# lec	10/27/2003
#	- created per TR 5252
#
'''
 
import sys
import os
import string
import db
import mgi_utils
import reportlib

CRT = reportlib.CRT
TAB = reportlib.TAB

problemNote = 'MGI curatorial staff have found evidence of artifact in the sequence of this molecular segment.'

#
# Main
#

fp = reportlib.init(sys.argv[0], printHeading = 0, outputdir = os.environ['REPORTOUTPUTDIR'])

cmds = []

cmds.append('select a1._Object_key, cloneID = a1.accID, mgiID = a2.accID ' + \
	'into #riken ' + \
	'from PRB_Acc_View a1, PRB_Acc_View a2 ' + \
	'where a1._LogicalDB_key = 26 ' + \
	'and a1._Object_key = a2._Object_key ' + \
	'and a2._LogicalDB_key = 1 ' + \
	'and a2.prefixPart = "MGI:" ' + \
	'and a2.preferred = 1')

cmds.append('create nonclustered index idx_key on #riken(_Object_key)')

cmds.append('select r._Object_key, a.accID ' + \
    'from #riken r, PRB_Acc_View a ' + \
    'where r._Object_key = a._Object_key ' + \
    'and a._LogicalDB_key = 9')

# problem sequences
cmds.append('select distinct r._Object_key ' + \
	'from #riken r, PRB_Notes n ' + \
	'where r._Object_key = n._Probe_key ' + \
	'and n.note like "%curatorial staff have found evidence of artifact in the sequence of this molecular%"')

cmds.append('select r.*, markerID = ma.accID, m.symbol, m.name ' + \
	'from #riken r, PRB_Marker pm, MRK_Acc_View ma, MRK_Marker m ' + \
	'where r._Object_key = pm._Probe_key ' + \
	'and pm._Marker_key = ma._Object_key ' + \
	'and ma._LogicalDB_key = 1 ' + \
	'and ma.prefixPart = "MGI:" ' + \
	'and ma.preferred = 1 ' + \
	'and pm._Marker_key = m._Marker_key ' + \
	'union ' + \
        'select r.*, markerID = null, symbol = null, name = null ' + \
	'from #riken r ' + \
	'where not exists (select 1 from PRB_Marker pm ' + \
	'where r._Object_key = pm._Probe_key) ' + \
	'order by r.cloneID')

results = db.sql(cmds, 'auto')

seqIDs = {}
for r in results[-3]:
    key = r['_Object_key']
    value = r['accID']
    seqIDs[key] = value

problemClones = []
for r in results[-2]:
    problemClones.append(r['_Object_key'])

for r in results[-1]:

	if r['_Object_key'] in problemClones:
	    isProblem = 'Problem'
	else:
	    isProblem = ''

	fp.write(r['cloneID'] + TAB + \
	    r['mgiID'] + TAB)

        if seqIDs.has_key(r['_Object_key']):
	    fp.write(seqIDs[r['_Object_key']])

        fp.write(TAB + isProblem + TAB + \
	    mgi_utils.prvalue(r['markerID']) + TAB + \
	    mgi_utils.prvalue(r['symbol']) + TAB + \
	    mgi_utils.prvalue(r['name']) + CRT)

reportlib.finish_nonps(fp)
