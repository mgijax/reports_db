#!/usr/local/bin/python

'''
#
# PRB_PrimerSeq.py 11/05/99
#
# Report:
#       TR 1047 Primers and Sequences
#
# Usage:
#       PRB_PrimerSeq.py
#
# Notes:
#	- all reports use mgireport directory for output file
#	- all reports use db default of public login
#	- all reports use server/database default of environment
#	- use lowercase for all SQL commands (i.e. select not SELECT)
#	- all public SQL reports require the header and footer
#	- all private SQL reports require the header
#
# History:
#
# lec	11/05/1999
#	- created
#
# dbm	02/11/2003
#	- added marker AccID, chromosome and offset for TR 4506
#
# dbm	02/12/2003
#	- added marker name
#
'''
 
import sys 
import os
import db
import reportlib
import string
import mgi_utils

CRT = reportlib.CRT
TAB = reportlib.TAB

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

db.sql('select p._Probe_key, pname = p.name, p.primer1sequence, p.primer2sequence, p.productSize, pm._Marker_key, ' + \
      'm.symbol, mname = m.name, m.chromosome ' + \
      'into #primers ' + \
      'from PRB_Probe p, PRB_Marker pm, MRK_Marker m ' + \
      'where p._SegmentType_key = 63473 and ' + \
      'p._Probe_key = pm._Probe_key and ' + \
      'pm._Marker_key = m._Marker_key', None)
db.sql('create index idx1 on #primers(_Probe_key)', None)
db.sql('create index idx2 on #primers(_Marker_key)', None)

cmd = 'select p.*, probeID = a1.accID, markerID = a2.accID, o.offset ' + \
      'from #primers p, ACC_Accession a1, ACC_Accession a2, MRK_Offset o ' + \
      'where p._Probe_key = a1._Object_key and ' + \
	    'a1._MGIType_key = 3 and ' + \
            'a1._LogicalDB_key = 1 and ' + \
            'a1.prefixPart = "MGI:" and ' + \
            'a1.preferred = 1 and ' + \
            'p._Marker_key = a2._Object_key and ' + \
	    'a2._MGIType_key = 2 and ' + \
            'a2._LogicalDB_key = 1 and ' + \
            'a2.prefixPart = "MGI:" and ' + \
            'a2.preferred = 1 and ' + \
            'p._Marker_key = o._Marker_key and ' + \
            'o.source = 0 ' + \
      'order by p.symbol'

results = db.sql(cmd, 'auto')

for r in results:
    mname = r['mname']
    pname = r['pname']
    p1seq = r['primer1sequence']
    p2seq = r['primer2sequence']
    prodSize = r['productSize']

    fp.write(r['symbol'] + TAB +
	mname + TAB +
	pname + TAB +
        r['markerID'] + TAB + 
	r['probeID'] + TAB +
        mgi_utils.prvalue(p1seq) + TAB + 
	mgi_utils.prvalue(p2seq) + TAB + 
	mgi_utils.prvalue(prodSize) + TAB +
        r['chromosome'] + TAB + 
	str(r['offset']) + CRT)

reportlib.finish_nonps(fp)
