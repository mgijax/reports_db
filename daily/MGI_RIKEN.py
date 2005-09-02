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
#	. Sequence ID
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
# lec	08/12/2005
#	- just select HTC sequences for RIKEN clones
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

fp.write('#\n')
fp.write('#  tab-delimited file of all MGI RIKEN clone/marker associations\n\n')
fp.write('#  column 1: RIKEN Clone ID\n')
fp.write('#  column 2: MGI Clone ID\n')
fp.write('#  column 3: Sequence ID\n')
fp.write('#  column 4: "Problem" Clone\n')
fp.write('#  column 5: MGI Marker ID\n')
fp.write('#  column 6: Gene Symbol\n')
fp.write('#  column 7: Gene Name\n')
fp.write('#\n')

#
# select RIKEN clones
#

db.sql('select a1._Object_key, cloneID = a1.accID, mgiID = a2.accID ' + \
	'into #riken ' + \
	'from ACC_Accession a1, ACC_Accession a2 ' + \
	'where a1._MGIType_key = 3 ' + \
	'and a1._LogicalDB_key = 26 ' + \
	'and a1._Object_key = a2._Object_key ' + \
	'and a2._MGIType_key = 3 ' + \
	'and a2._LogicalDB_key = 1 ' + \
	'and a2.prefixPart = "MGI:" ' + \
	'and a2.preferred = 1', None)

db.sql('create nonclustered index idx_key on #riken(_Object_key)', None)

#
# select GenBank/EMBL/DDBJ:HTC Riken Sequence
# or GenBank/EMBL/DDBJ (dummy Sequence)
# there should be at most one per RIKEN clone
#

results = db.sql('select r._Object_key, a.accID ' + \
    'from #riken r, SEQ_Probe_Cache p, SEQ_Sequence s, ACC_Accession a ' + \
    'where r._Object_key = p._Probe_key ' + \
    'and p._Sequence_key = s._Sequence_key ' + \
    'and s._SequenceProvider_key in (316375, 316380) ' + \
    'and p._Sequence_key = a._Object_key ' + \
    'and a._MGIType_key = 19 ' + \
    'and a._LogicalDB_key = 9', 'auto')

seqIDs = {}
for r in results:
    key = r['_Object_key']
    value = r['accID']
    seqIDs[key] = value

#
# problem sequences
#

results = db.sql('select distinct r._Object_key ' + \
	'from #riken r, PRB_Notes n ' + \
	'where r._Object_key = n._Probe_key ' + \
	'and n.note like "%curatorial staff have found evidence of artifact in the sequence of this molecular%"', 'auto')

problemClones = []
for r in results:
    problemClones.append(r['_Object_key'])

#
# final results
#

results = db.sql('select r.*, markerID = ma.accID, m.symbol, m.name ' + \
	'from #riken r, PRB_Marker pm, ACC_Accession ma, MRK_Marker m ' + \
	'where r._Object_key = pm._Probe_key ' + \
	'and pm._Marker_key = ma._Object_key ' + \
	'and ma._MGIType_key = 2 ' + \
	'and ma._LogicalDB_key = 1 ' + \
	'and ma.prefixPart = "MGI:" ' + \
	'and ma.preferred = 1 ' + \
	'and pm._Marker_key = m._Marker_key ' + \
	'union ' + \
        'select r.*, markerID = null, symbol = null, name = null ' + \
	'from #riken r ' + \
	'where not exists (select 1 from PRB_Marker pm ' + \
	'where r._Object_key = pm._Probe_key) ' + \
	'order by r.cloneID', 'auto')

for r in results:

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

