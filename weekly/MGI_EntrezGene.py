#!/usr/local/bin/python

'''
#
# Report:
#       Tab-delimited file of MGI Mouse Markers including Withdrawns
#	Prints current MGI accession ID of Marker record
#	and list of non-preferred MGI accession IDs of each Marker
#
# Usage:
#       MGI_EntrezGene.py
#
# Used by:
# 	TR 1283 - Donna Maglott at NCBI for EntrezGene
#
# Notes:
#
# Splits will retain their MGI Accession IDs.  The 3rd query is to find
# withdrawns which still have a preferred MGI Accession ID.
#
# History:
#
# lec	06/18/2002
#	- rewrote to use dictionaries for egID, other Acc IDs, other names
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

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# 1. all Approved Marker records
# 2. all Withdrawn Marker records
#

db.sql('select m._Marker_key, m.symbol, m.name, _Current_key = m._Marker_key, ' + \
  'offset = str(o.offset,10,2), m.chromosome, markerType = t.name, isPrimary = 1, ' + \
  'markerStatus = upper(substring(s.status, 1, 1)) ' + \
  'into #markers ' + \
  'from MRK_Marker m, MRK_Offset o, MRK_Types t, MRK_Status s ' + \
  'where m._Organism_key = 1 ' + \
  'and m._Marker_Status_key in (1,3) ' + \
  'and m._Marker_key = o._Marker_key ' + \
  'and o.source = 0 ' + \
  'and m._Marker_Type_key = t._Marker_Type_key ' + \
  'and m._Marker_Status_key = s._Marker_Status_key ' + \
  'union '  + \
  'select m._Marker_key, m.symbol, m.name, c._Current_key, ' + \
  'offset = str(o.offset,10,2), m.chromosome, markerType = t.name, isPrimary = 0, ' + \
  'markerStatus = upper(substring(s.status, 1, 1)) ' + \
  'from MRK_Marker m, MRK_Offset o, MRK_Current c, MRK_Types t, MRK_Status s ' + \
  'where m._Organism_key = 1 ' + \
  'and m._Marker_Status_key = 2 ' + \
  'and m._Marker_key = o._Marker_key ' + \
  'and o.source = 0 ' + \
  'and m._Marker_key = c._Marker_key ' + \
  'and m._Marker_Type_key = t._Marker_Type_key ' + \
  'and m._Marker_Status_key = s._Marker_Status_key ', None)
db.sql('create nonclustered index idx_key on #markers(_Marker_key)', None)

# MGI ids

results = db.sql('select m._Current_key, a.accID ' + \
	'from #markers m, ACC_Accession a ' + \
  	'where m._Current_key = a._Object_key ' + \
  	'and a._MGIType_key = 2 ' + \
  	'and a.prefixPart = "MGI:" ' + \
  	'and a._LogicalDB_key = 1 ' + \
  	'and a.preferred = 1 ', 'auto')
mgiID = {}
for r in results:
    key = r['_Current_key']
    value = r['accID']
    mgiID[key] = value

# Get EntrezGene ID for Primary Markers
results = db.sql('select m._Marker_key, a.accID ' + \
	'from #markers m, ACC_Accession a ' + \
	'where m.isPrimary = 1 ' + \
	'and m._Marker_key = a._Object_key ' + \
        'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 55 ', 'auto')
egID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    egID[key] = value

# Get Secondary MGI Ids for Primary Marker
results = db.sql('select m._Marker_key, a.accID ' + \
	'from #markers m, ACC_Accession a ' + \
	'where m.isPrimary = 1 ' + \
	'and m._Marker_key = a._Object_key ' + \
        'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
        'and a._LogicalDB_key = 1 ' + \
	'and a.preferred = 0', 'auto')
otherAccId = {}
for r in results:
	if not otherAccId.has_key(r['_Marker_key']):
		otherAccId[r['_Marker_key']] = []
	otherAccId[r['_Marker_key']].append(r['accID'])

# Get Synonyms for Primary Marker
results = db.sql('select m._Marker_key, s.synonym ' + \
	'from #markers m, MGI_Synonym s, MGI_SynonymType st ' + \
	'where m.isPrimary = 1 ' + \
	'and m._Marker_key = s._Object_key ' + \
	'and s._MGIType_key = 2 ' + \
	'and s._SynonymType_key = st._SynonymType_key ' + 
	'and st.synonymType = "exact"', 'auto')
synonym = {}
for r in results:
	key = r['_Marker_key']
	value = r['synonym']
	if not synonym.has_key(key):
		synonym[key] = []
	synonym[key].append(value)

results = db.sql('select * from #markers order by _Current_key, isPrimary desc', 'auto')
for r in results:
	fp.write(mgiID[r['_Current_key']] + TAB + \
	 	r['symbol'] + TAB + \
		r['markerStatus'] + TAB + \
	 	r['name'] + TAB + \
	 	r['offset'] + TAB + \
	 	r['chromosome'] + TAB + \
		r['markerType'] + TAB)

	if r['isPrimary']:
		if otherAccId.has_key(r['_Marker_key']):	
			fp.write(string.joinfields(otherAccId[r['_Marker_key']], '|'))
		fp.write(TAB)

		if egID.has_key(r['_Marker_key']):	
			fp.write(egID[r['_Marker_key']])
		fp.write(TAB)

		if synonym.has_key(r['_Marker_key']):	
			fp.write(string.joinfields(synonym[r['_Marker_key']], '|'))
		fp.write(CRT)
	else:
		fp.write(TAB)
		fp.write(TAB)
		fp.write(CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
