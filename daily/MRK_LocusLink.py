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
# lec	06/18/2002
#	- rewrote to use dictionaries for locusID, other Acc IDs, other names
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

CRT = reportlib.CRT
TAB = reportlib.TAB

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []

#
# 1. all Approved Marker records
# 2. all Withdrawn Marker records
#

cmds.append('select m._Marker_key, m.symbol, m.name, _Current_key = m._Marker_key, offset = str(o.offset,10,2), ' + \
  'm.chromosome, a.accID, markerType = t.name, isPrimary = 1 ' + \
  'into #markers ' + \
  'from MRK_Marker m, MRK_Offset o, MRK_Acc_View a, MRK_Types t ' + \
  'where m._Species_key = 1 ' + \
  'and m._Marker_Status_key = 1 ' + \
  'and m._Marker_key = o._Marker_key ' + \
  'and o.source = 0 ' + \
  'and m._Marker_key = a._Object_key ' + \
  'and a.prefixPart = "MGI:" ' + \
  'and a.preferred = 1 ' + \
  'and m._Marker_Type_key = t._Marker_Type_key ' + \
  'union '  + \
  'select m._Marker_key, m.symbol, m.name, c._Current_key, offset = str(o.offset,10,2), m.chromosome, ' + \
  'a.accID, markerType = t.name, isPrimary = 0 ' + \
  'from MRK_Marker m, MRK_Offset o, MRK_Current c, MRK_Acc_View a, MRK_Types t ' + \
  'where m._Species_key = 1 ' + \
  'and m._Marker_Status_key = 2 ' + \
  'and m._Marker_key = o._Marker_key ' + \
  'and o.source = 0 ' + \
  'and m._Marker_key = c._Marker_key ' + \
  'and c._Current_key = a._Object_key ' + \
  'and a.prefixPart = "MGI:" ' + \
  'and a.preferred = 1 ' + \
  'and m._Marker_Type_key = t._Marker_Type_key')

cmds.append('create nonclustered index idx_key on #markers(_Marker_key)')

# Get Locus ID for Primary Markers
cmds.append('select m._Marker_key, a.accID ' + \
	'from #markers m, MRK_Acc_View a ' + \
	'where m.isPrimary = 1 ' + \
	'and m._Marker_key = a._Object_key ' + \
	'and a._LogicalDB_key = 24 ')

# Get Secondary MGI Ids for Primary Marker
cmds.append('select m._Marker_key, a.accID ' + \
	'from #markers m, MRK_Acc_View a ' + \
	'where m.isPrimary = 1 ' + \
	'and m._Marker_key = a._Object_key ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 0')

# Get Synonyms for Primary Marker
cmds.append('select m._Marker_key, o.name ' + \
	'from #markers m, MRK_Other o ' + \
	'where m.isPrimary = 1 ' + \
	'and m._Marker_key = o._Marker_key ')

cmds.append('select * from #markers order by _Current_key, isPrimary desc')

results = db.sql(cmds, 'auto')

locusID = {}
for r in results[2]:
	if not locusID.has_key(r['_Marker_key']):
		locusID[r['_Marker_key']] = r['accID']

otherAccId = {}
for r in results[3]:
	if not otherAccId.has_key(r['_Marker_key']):
		otherAccId[r['_Marker_key']] = []
	otherAccId[r['_Marker_key']].append(r['accID'])

otherName = {}
for r in results[4]:
	if not otherName.has_key(r['_Marker_key']):
		otherName[r['_Marker_key']] = []
	otherName[r['_Marker_key']].append(r['name'])

for r in results[5]:

	fp.write(r['accID'] + TAB + \
	 	r['symbol'] + TAB + \
	 	r['name'] + TAB + \
	 	r['offset'] + TAB + \
	 	r['chromosome'] + TAB + \
		r['markerType'] + TAB)

	if r['isPrimary']:
		if otherAccId.has_key(r['_Marker_key']):	
			fp.write(string.joinfields(otherAccId[r['_Marker_key']], '|'))
		fp.write(TAB)

		if locusID.has_key(r['_Marker_key']):	
			fp.write(locusID[r['_Marker_key']])
		fp.write(TAB)

		if otherName.has_key(r['_Marker_key']):	
			fp.write(string.joinfields(otherName[r['_Marker_key']], '|'))
		fp.write(CRT)
	else:
		fp.write(TAB)
		fp.write(TAB)
		fp.write(CRT)

reportlib.finish_nonps(fp)

