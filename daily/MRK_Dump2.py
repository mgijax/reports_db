#!/usr/local/bin/python

'''
#
# MRK_Dump2.py 11/16/98
#
# Report:
#       Tab-delimited file of MGI Mouse Markers
#       excluding Withdrawns Symbols
#
# Usage:
#       MRK_Dump2.py
#
# Used by:
#
# Notes:
#
# History:
#
# 05/30/2002	lec
#	- TR 3736; add Marker Type
#
# 01/13/98	lec
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
cmds.append('select m._Marker_key, m.symbol, m.name, m.chromosome, markerType = t.term, offset_str = str(o.offset,10,2) ' + \
	'into #markers ' + \
	'from MRK_Marker m, VOC_Term t, MRK_Offset o ' + \
	'where m._Organism_key = 1 ' + \
	'and m._Marker_Status_key in (1,3) ' + \
	'and m._Marker_Type_key = t._Term_key ' + \
	'and m._Marker_key = o._Marker_key ' + \
	'and o.source = 0')
cmds.append('create index idx1 on #markers(_Marker_key)')
cmds.append('create index idx2 on #markers(symbol)')
db.sql(cmds, None)

results = db.sql('select m.symbol, m.symbol, m.name, m.chromosome, m.markerType, m.offset_str, a.accID ' + \
	'from #markers m, ACC_Accession a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1 ' + \
	'order by m.symbol', 'auto')

for r in results:
	fp.write(r['accID'] + reportlib.TAB + \
	         r['symbol'] + reportlib.TAB + \
		 r['name'] + reportlib.TAB + \
		 r['offset_str'] + reportlib.TAB + \
		 r['chromosome'] + reportlib.TAB + \
		 r['markerType'] + reportlib.CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)

