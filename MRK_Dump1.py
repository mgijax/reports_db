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
# Generated from:
#       Editing Interface Nightly Reports
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

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

command = 'select symbol, mgiID ' + \
	  'from MRK_Mouse_View ' + \
	  'where _Marker_Status_key = 1 ' + \
	  'order by symbol'
results = db.sql(command, 'auto')

for r in results:
	fp.write(r['mgiID'] + reportlib.TAB + \
	         r['symbol'] + reportlib.CRT)

reportlib.finish_nonps(fp)

