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
import mgi_utils

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []

cmds.append('select m._Marker_key, m.symbol, m.name, c._Current_key, offset = str(o.offset,10,2), m.chromosome, a.accID, markerType = t.name, isPrimary = 1 ' + \
  'into #markers ' + \
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
  'and m._Marker_Type_key = t._Marker_Type_key')

cmds.append('select m.*, locusID = a.accID ' + \
	'into #markers2 ' + \
	'from #markers m, MRK_Acc_View a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._LogicalDB_key = 24 ' + \
        'union ' + \
	'select m.*, locusID = null ' + \
	'from #markers m ' + \
	'where not exists (select 1 from MRK_Acc_View a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._LogicalDB_key = 24)')

cmds.append('select m.*, otherName = o.name ' + \
	'from #markers2 m, MRK_Other o ' + \
	'where m.isPrimary = 1 ' + \
	'and m._Marker_key *= o._Marker_key ' + \
	'union ' +
	'select m.*, otherName = null ' + \
	'from #markers2 m ' + \
	'where m.isPrimary = 0 ' + \
	'order by _Current_key, isPrimary desc')

results = db.sql(cmds, 'auto')

prevMarker = ''
otherAccIds = []
otherNames = []
locusID = ''
num = 0

for r in results[2]:

	if r['isPrimary']:

		if prevMarker != r['_Marker_key']:

			if num > 0:
				fp.write(string.join(otherAccIds, ',') + reportlib.TAB)
				fp.write(mgi_utils.prvalue(locusID) + reportlib.TAB) 
				fp.write(string.join(otherNames, '|') + reportlib.CRT)

			fp.write(r['accID'] + reportlib.TAB + \
	         		r['symbol'] + reportlib.TAB + \
		 		r['name'] + reportlib.TAB + \
		 		r['offset'] + reportlib.TAB + \
		 		r['chromosome'] + reportlib.TAB + \
				r['markerType'] + reportlib.TAB)

			prevMarker = r['_Marker_key']
			otherAccIds = []
			otherNames = []
			locusID = r['locusID']
			num = num + 1

		if r['otherName'] is not None:
			otherNames.append(r['otherName'])
	else:
		if r['accID'] not in otherAccIds:
			otherAccIds.append(r['accID'])

fp.write(string.join(otherAccIds, ',') + reportlib.TAB)
fp.write(mgi_utils.prvalue(locusID) + reportlib.TAB) 
fp.write(string.join(otherNames, ',') + reportlib.CRT)

reportlib.finish_nonps(fp)

