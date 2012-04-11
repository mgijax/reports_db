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
# 04/11/2012	lec
#	- TR11035/postgres options
#
# 12/28/2011	lec
#	- changed non-ansi-standard query to left outer join
#
# 05/20/2010	lec
#	- TR10049/turn TR9283 into a public report
#
# 2008-09-26	Susan McClatchy
#	- created
"""

import sys
import os
import string
import reportlib

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db


CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE


#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('Mouse MGI Accession ID' + TAB)
fp.write('Mouse Marker Symbol' + TAB)
fp.write('Mouse Entrez Gene ID' + TAB)
fp.write('Human Entrez Gene ID' + TAB)
fp.write('Mouse Marker Synonym' + TAB)
fp.write('Human Marker Symbol' + CRT)

#
# Retrieve all official/interim mouse genes.
#

db.sql('''
	select distinct m._Marker_key, m.mgiID, m.symbol
	into #mousehuman 
	from MRK_Mouse_View m 
    	where m._Marker_Type_key = 1 
	and m._Marker_Status_key in (1,3) 
    	''', None)
db.sql('create index idx1 on #mousehuman(_Marker_key)', None)

# mouse synonyms
results = db.sql('''
	select distinct m._Marker_key, s.synonym 
	from #mousehuman m, MGI_Synonym_MusMarker_View s
        where m._Marker_key = s._Object_key
        and s._SynonymType_key = 1004
	''', 'auto')
syns = {}
for r in results:
	key = r['_Marker_key']
	value = r['synonym']
	if not syns.has_key(key):
		syns[key] = []
	syns[key].append(value)

# mouse entrezgene ids
results = db.sql('''
	select h._Marker_key, a.accID 
	from #mousehuman h, ACC_Accession a 
	where h._Marker_key = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 55
	''', 'auto')
mouseEG = {}
for r in results:
	key = r['_Marker_key']
	value = r['accID']
	mouseEG[key] = value

# human symbols and entrezgene ids
results = db.sql('''
	select distinct m._Marker_key, h.marker2, a.accID
	from #mousehuman m, HMD_Homology_Pairs_View h, ACC_Accession a 
	where m._Marker_key = h.markerkey1 
	and h.organismkey2 = 2
	and h.markerkey2 = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 55
	''', 'auto')
human = {}
humanEG = {}
for r in results:
	key = r['_Marker_key']
	value = r['marker2']
	human[key] = value
	value = r['accID']
	humanEG[key] = value

results = db.sql('select * from #mousehuman order by mgiID', 'auto')
for r in results:

	fp.write(r['mgiID'] + TAB + r['symbol'] + TAB)

	if mouseEG.has_key(r['_Marker_key']):
		fp.write(mouseEG[r['_Marker_key']])
	fp.write(TAB)

	if humanEG.has_key(r['_Marker_key']):
		fp.write(humanEG[r['_Marker_key']])
	fp.write(TAB)

	if syns.has_key(r['_Marker_key']):
		fp.write(string.join(syns[r['_Marker_key']], ','))
	fp.write(TAB)

	if human.has_key(r['_Marker_key']):
		fp.write(human[r['_Marker_key']])
	fp.write(CRT)

reportlib.finish_nonps(fp)
