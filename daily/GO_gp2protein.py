#!/usr/local/bin/python

'''
#
# GO_gp2protein.py (TR 4877, originally TR 3659)
#
# Report:
#       Tab-delimited file
#
# Usage:
#       GO_gp2protein.py
#
# Used by:
#
# Notes:
#
# History:
#	Created for TR3659.  
#	Report the MGI ID, SwissProt ID(s), and GO ID(s) for all mouse
#	markers that have associated SwissProt sequences.

'''
 
import sys
import string
import os
import db
import reportlib

#
# Main
#

TAB = reportlib.TAB
CRT = reportlib.CRT

fp = reportlib.init('gp2protein', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

# Retrieve Markers with GO Annotations that have SP IDs

cmd = 'select distinct m._Marker_key, a.accID, spID = a2.accID ' + \
	'from MRK_Marker m, MRK_Acc_View a, MRK_Acc_View a2 ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1 ' + \
	'and m._Marker_key = a2._Object_key ' + \
	'and a2._LogicalDB_key = 13 ' + \
	'and exists (select 1 from VOC_Annot v ' + \
	'where v._AnnotType_key = 1000 ' + \
	'and m._Marker_key = v._Object_key) ' + \
	'order by a.accID'

results = db.sql(cmd, 'auto')

for r in results:

    fp.write(r['accID'] + reportlib.TAB + 'SP:' + r['spID'] + reportlib.CRT)

reportlib.finish_nonps(fp)

