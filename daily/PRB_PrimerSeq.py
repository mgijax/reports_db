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
'''
 
import sys 
import os
import db
import reportlib

CRT = reportlib.CRT
TAB = reportlib.TAB

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmd = 'select p.name, p.primer1sequence, p.primer2sequence, p.productSize, m.symbol, p.mgiID ' + \
	'from PRB_Primer_View p, PRB_Marker_View m ' + \
	'where p._Probe_key = m._Probe_key ' + \
	'order by m.symbol'

results = db.sql(cmd, 'auto')

fp.write('%-40s' % ("name") + TAB)
fp.write('%-25s' % ("symbol") + TAB)
fp.write('%-80s' % ("primer1sequence") + TAB)
fp.write('%-80s' % ("primer2sequence") + TAB)
fp.write('%-40s' % ("productSize") + TAB)
fp.write('%-30s' % ("mgiID") + CRT)

for r in results:
	fp.write('%-40s' % (r['name']) + TAB)
	fp.write('%-25s' % (r['symbol']) + TAB)
	fp.write('%-80s' % (str(r['primer1sequence'])) + TAB)
	fp.write('%-80s' % (str(r['primer2sequence'])) + TAB)
	fp.write('%-40s' % (str(r['productSize'])) + TAB)
	fp.write('%-30s' % (r['mgiID']) + CRT)

reportlib.trailer(fp)
reportlib.finish_nonps(fp)	# non-postscript file

