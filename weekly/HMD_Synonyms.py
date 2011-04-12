#!/usr/local/bin/python
"""
HMD_Synonyms.py

# Report:
# 
# Produce a "Gene Dictionary" file in tab-delimited format that includes, for 
# markers with type 'Gene'
#
# MGI ID
# Mouse Symbol
# Mouse Entrez Gene ID
# Human Entrez Gene ID
# Mouse Synonyms
# Human Symbol
#
# Usage:
#       HMD_Synonyms.py
#
# History:
#
# 05/20/2010	lec
#	- TR10049/turn TR9283 into a public report
#
# 2008-09-26	Susan McClatchy
#	- created
"""

import sys
import os
import db
import reportlib
import string


CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE


#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# Retrieve all official/interim mouse genes.
#

db.sql('select distinct m._Marker_key, m.mgiID, m.symbol, s.synonym, ' + \
	'markerkey2 = 0, marker2 = "                              " ' + \
	'into #mousehuman ' + \
	'from MRK_Mouse_View m, MGI_Synonym_MusMarker_View s ' + \
    'where m._Marker_Type_key = 1 ' + \
    'and m._Marker_Status_key in (1,3) ' + \
    'and m._Marker_key *= s._Object_key ' + \
    'and s._SynonymType_key = 1004', None)

db.sql('update #mousehuman ' + \
 	'set synonym = "none" ' + \
 	'from #mousehuman ' + \
 	'where synonym is null', None)

db.sql('create index idx1 on #mousehuman(_Marker_key)', None)

# mouse synonyms

results = db.sql('select distinct _Marker_key, synonym ' + \
	'from #mousehuman', 'auto')

syns = {}

for r in results:
	key = r['_Marker_key']
	value = r['synonym']
	if not syns.has_key(key):
		syns[key] = []
	syns[key].append(value)

#
# Add human ortholog marker keys and symbols where applicable.
#

db.sql('update #mousehuman ' + \
	'set m.markerkey2 = h.markerkey2, ' + \
	'm.marker2 = h.marker2 ' + \
	'from #mousehuman m, HMD_Homology_Pairs_View h ' + \
	'where m._Marker_key = h.markerkey1 ' + \
	'and h.organismkey2 = 2', 'auto')

# mouse entrezgene ids

results = db.sql('select h._Marker_key, a.accID ' + \
	'from #mousehuman h, ACC_Accession a ' + \
	'where h._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 55', 'auto')

mouseEG = {}

for r in results:
	key = r['_Marker_key']
	value = r['accID']
	mouseEG[key] = value

# human entrezgene ids

results = db.sql('select h.markerkey2, a.accID ' + \
	'from #mousehuman h, ACC_Accession a ' + \
	'where h.markerkey2 = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 55', 'auto')

humanEG = {}

for r in results:
	key = r['markerkey2']
	value = r['accID']
	humanEG[key] = value


results = db.sql('select * from #mousehuman h', 'auto')

fp.write('Mouse MGI Accession ID' + TAB)
fp.write('Mouse Marker Symbol' + TAB)
fp.write('Mouse Entrez Gene ID' + TAB)
fp.write('Human Entrez Gene ID' + TAB)
fp.write('Mouse Marker Synonym' + TAB)
fp.write('Human Marker Symbol' + CRT)

for r in results:

	fp.write(r['mgiID'] + TAB + r['symbol'] + TAB)

	if mouseEG.has_key(r['_Marker_key']):
		fp.write(mouseEG[r['_Marker_key']])
	fp.write(TAB)

	if humanEG.has_key(r['markerkey2']):
		fp.write(humanEG[r['markerkey2']])
	fp.write(TAB)

	if syns.has_key(r['_Marker_key']):
		fp.write(string.join(syns[r['_Marker_key']], ','))
	fp.write(TAB)

	fp.write(r['marker2'] + CRT)

reportlib.finish_nonps(fp)
