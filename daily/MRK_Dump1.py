#!/usr/local/bin/python

'''
#
# MRK_Dump1.py 11/16/98
#
# Report:
#       Tab-delimited file of MGI Mouse Markers
#	excluding Withdrawns Symbols
#
# Usage:
#       MRK_Dump1.py
#
# Used by:
#	Those who want to create WWW links to MGI Marker details.
#
# Notes:
#
# History:
#
# lec	01/13/98
#	- added comments
#
'''

import sys
import os
import db
import reportlib

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []
cmds.append('select _Marker_key, symbol ' + \
	'into #markers ' + \
	'from MRK_Marker ' + \
	'where _Organism_key = 1 ' + \
	'and _Marker_Status_key in (1,3)')
cmds.append('create index idx1 on #markers(_Marker_key)')
cmds.append('create index idx2 on #markers(symbol)')
db.sql(cmds, None)

results = db.sql('select m.symbol, a.accID ' + \
	'from #markers m, ACC_Accession a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1 ' + \
	'order by m.symbol', 'auto')

for r in results:
	fp.write(r['accID'] + reportlib.TAB + \
	         r['symbol'] + reportlib.CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)

