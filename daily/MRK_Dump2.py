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

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

command = 'select symbol, name, offset_str, chromosome, mgiID, markerType ' + \
	  'from MRK_Mouse_View ' + \
	  'where _Marker_Status_key in (1,3) ' + \
	  'order by symbol'
results = db.sql(command, 'auto')

for r in results:
	fp.write(r['mgiID'] + reportlib.TAB + \
	         r['symbol'] + reportlib.TAB + \
		 r['name'] + reportlib.TAB + \
		 r['offset_str'] + reportlib.TAB + \
		 r['chromosome'] + reportlib.TAB + \
		 r['markerType'] + reportlib.CRT)

reportlib.finish_nonps(fp)

