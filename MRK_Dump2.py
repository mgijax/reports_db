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
# Generated from:
#       Editing Interface Nightly Reports
#
# Used by:
# 	Stan Attenberger at sea@ornl.gov; 
# 	using Java applet to display deletions 
#	next to genetic map of mouse.
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

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

command = 'select symbol, name, offset_str, chromosome, mgiID ' + \
	  'from MRK_Mouse_View ' + \
	  'where _Marker_Status_key = 1 ' + \
	  'order by symbol'
results = db.sql(command, 'auto')

for r in results:
	fp.write(r['mgiID'] + reportlib.TAB + \
	         r['symbol'] + reportlib.TAB + \
		 r['name'] + reportlib.TAB + \
		 r['offset_str'] + reportlib.TAB + \
		 r['chromosome'] + reportlib.CRT)

reportlib.finish_nonps(fp)

