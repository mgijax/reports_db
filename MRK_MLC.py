#!/usr/local/bin/python

'''
#
# MRK_MLC.py 11/16/98
#
# Report:
#	Tab-delimited file
#       Mouse Markers and their MLC Classes
#       Excluding withdrawns and reserved
#
# Usage:
#       MRK_MLC.py
#
# Generated from:
#       Editing Interface Nightly Reports scripts (nightly_reports)
#
# Used by:
# 	Stan Attenberger at sea@ornl.gov; 
# 	using Java applet to display deletions next to genetic map of mouse.
#
# Notes:
#
# History:
#
# lec	01/13/98
#	- added comments section
#
'''

import sys
import os
import db
import reportlib

CRT = reportlib.CRT

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

command = 'select m.mgiID, m.symbol, l.name ' + \
	  'from MRK_Mouse_View m, MRK_Classes c, MRK_Class l ' + \
	  'where m._Marker_Status_key = 1 ' + \
	  'and m._Marker_key = c._Marker_key ' + \
	  'and c._Class_key = l._Class_key ' + \
	  'order by m.symbol'
results = db.sql(command, 'auto')

for r in results:
	fp.write(r['symbol'] + reportlib.TAB + \
		 r['name'] + reportlib.TAB + \
		 r['mgiID'] + CRT)

reportlib.finish_nonps(fp)

