#!/usr/local/bin/python

'''
#
# Report:
#       Tab-delimited file of MGI Mouse Markers including Withdrawns
#	Prints current MGI accession ID of Marker record
#	and list of non-preferred MGI accession IDs of each Marker
#
# Usage:
#       MRK_LocusLink.py
#
# Generated from:
#       Editing Interface Nightly Reports
#
# Used by:
# 	TR 1283 - Donna Maglott at NCBI for LocusLink
#
# Notes:
#
# Splits will retain their MGI Accession IDs.  The 3rd query is to find
# withdrawns which still have a preferred MGI Accession ID.
#
# History:
#
# lec	01/19/2000
#	- created
#
'''
 
import sys
import os
import string
import db
import reportlib

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmd = 'select m._Marker_key, m.symbol, m.name, c._Current_key, offset = str(o.offset,10,2), m.chromosome, a.accID, markerType = t.name, isPrimary = 1 ' + \
  'from MRK_Marker m, MRK_Offset o, MRK_Current c, MRK_Acc_View a, MRK_Types t ' + \
  'where m._Species_key = 1 ' + \
  'and m._Marker_key = o._Marker_key ' + \
  'and o.source = 0 ' + \
  'and m._Marker_key = c._Marker_key ' + \
  'and c._Current_key = a._Object_key ' + \
  'and a.prefixPart = "MGI:" ' + \
  'and a.preferred = 1 ' + \
  'and m._Marker_Type_key = t._Marker_Type_key ' + \
  'union ' + \
  'select m._Marker_key, m.symbol, m.name, c._Current_key, offset = null, m.chromosome, a.accID, markerType = t.name, isPrimary = 0 ' + \
  'from MRK_Marker m, MRK_Current c, MRK_Acc_View a, MRK_Types t ' + \
  'where m._Species_key = 1 ' + \
  'and m._Marker_key = c._Marker_key ' + \
  'and c._Current_key = a._Object_key ' + \
  'and a.prefixPart = "MGI:" ' + \
  'and a.preferred = 0 ' + \
  'and m._Marker_Type_key = t._Marker_Type_key ' + \
  'union '  + \
  'select m._Marker_key, m.symbol, m.name, c._Current_key, offset = null, m.chromosome, a.accID, markerType = t.name, isPrimary = 0 ' + \
  'from MRK_Marker m, MRK_Current c, MRK_Acc_View a, MRK_Types t ' + \
  'where m._Species_key = 1 ' + \
  'and m._Marker_Status_key = 2' + \
  'and m._Marker_key = c._Marker_key ' + \
  'and m._Marker_key = a._Object_key ' + \
  'and a.prefixPart = "MGI:" ' + \
  'and a.preferred = 1 ' + \
  'and m._Marker_Type_key = t._Marker_Type_key ' + \
  'order by _Current_key, isPrimary desc'
results = db.sql(cmd, 'auto')

otherAccIds = []
num = 0

for r in results:
	if r['isPrimary']:

		if len(otherAccIds) > 0 and num > 0:
			fp.write(string.join(otherAccIds, ','))

		if num > 0:
			fp.write(reportlib.CRT)

		fp.write(r['accID'] + reportlib.TAB + \
	         	r['symbol'] + reportlib.TAB + \
		 	r['name'] + reportlib.TAB + \
		 	r['offset'] + reportlib.TAB + \
		 	r['chromosome'] + reportlib.TAB + \
			r['markerType'] + reportlib.TAB)

		otherAccIds = []
		num = num + 1

	else:
		if r['accID'] not in otherAccIds:
			otherAccIds.append(r['accID'])

if len(otherAccIds) > 0:
	fp.write(string.join(otherAccIds, ',') + reportlib.CRT)

reportlib.finish_nonps(fp)

