#!/usr/local/bin/python

'''
#
# Report:
#       Tab-delimited file
#       Mouse Markers and their Interpro Acc Id's
#
# Usage:
#       MRK_InterPro.py
#
#
# Notes:
#
# History:
#
# lec	04/22/2003
#	- TR 3702; InterPro annotations are now made via the VOC Annotation tables
#
# lec	12/31/2002
#	- made a regular nightly report
#
# pm	11/21/2001
#
'''
 
import sys
import os
import db
import mgi_utils
import reportlib

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

# Retrieve MGI Accession number, InterPro accession Id's.

cmd = 'select distinct mgiID = a1.accID, m.symbol, a2.accID ' + \
      'from MRK_Marker m, ACC_Accession a1, VOC_Annot a, ACC_Accession a2 ' + \
      'where m._Organism_key = 1 ' + \
      'and m._Marker_key = a1._Object_key ' + \
      'and a1._MGIType_key = 2 ' + \
      'and a1._LogicalDB_key = 1 ' + \
      'and a1.prefixPart = "MGI:" ' + \
      'and a1.preferred = 1 ' + \
      'and m._Marker_key = a._Object_key ' + \
      'and a._AnnotType_key = 1003 ' + \
      'and a._Term_key = a2._Object_key ' + \
      'and a2._MGIType_key = 13 ' + \
      'and a2._LogicalDB_key = 28 ' + \
      'order by m.symbol '

results = db.sql(cmd, 'auto')

prevMarker = ''

for r in results:

	if prevMarker == r['mgiID']:
		fp.write(' ' + r['accID'])
                prevMarker = r['mgiID']
        else: 
		if prevMarker != '':
	        	fp.write(reportlib.CRT)
		fp.write(r['mgiID'] + reportlib.TAB)
		fp.write(r['symbol'] + reportlib.TAB)
                fp.write(r['accID'])
                prevMarker = r['mgiID']

fp.write(reportlib.CRT)
reportlib.finish_nonps(fp)

