#!/usr/local/bin/python

'''
#
# MRK_Reference.py 01/23/2003
#
# Report:
#       TR 4454
#	Tab-delimited file of:
#		Mouse Marker MGI ID
#		Symbol
#		Name
#		Synonyms
#		PubMed IDs for References 
#	Sorted by Symbol
#
# Usage:
#       MRK_Reference.py
#
# History:
#
# lec	01/23/2003
#	- created
#
'''
 
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


TAB = reportlib.TAB
CRT = reportlib.CRT

#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# mouse markers
#
db.sql('select a.accID, m._Marker_key, m.symbol, m.name ' + \
	'into #markers ' + \
	'from MRK_Marker m, ACC_Accession a ' + \
	'where m._Organism_key = 1 ' + \
	'and m._Marker_Status_key in (1,3) ' + \
	'and m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1', None)
db.sql('create index idx_marker on #markers(_Marker_key)', None)

#
# references
#
db.sql('select distinct m._Marker_key, r._Refs_key ' + \
	'into #references ' + \
	'from #markers m, MRK_Reference r ' + \
	'where m._Marker_key = r._Marker_key', None)
db.sql('create index idx_refs on #references(_Refs_key)', None)

#
# pub med ids
#
results = db.sql('select r._Marker_key, a.accID ' + \
	'from #references r, ACC_Accession a ' + \
	'where a._MGIType_key = 1 ' + \
	'and r._Refs_key = a._Object_key ' + \
	'and a._LogicalDB_key = 29', 'auto')
pubmed = {}
for r in results:
	key = r['_Marker_key']
	if not pubmed.has_key(key):
		pubmed[key] = []
	pubmed[key].append(r['accID'])

#
# synonyms
#
results = db.sql('select m._Marker_key, s.synonym ' + \
	'from #markers m, MGI_Synonym s, MGI_SynonymType st ' + \
	'where m._Marker_key = s._Object_key ' + \
	'and s._MGIType_key = 2 ' + \
	'and s._SynonymType_key = st._SynonymType_key ' + \
	'and st.synonymType = "exact"', 'auto')
syn = {}
for r in results:
	key = r['_Marker_key']
	if not syn.has_key(key):
		syn[key] = []
	syn[key].append(r['synonym'])

#
# final results
#
results = db.sql('select * from #markers order by symbol', 'auto')
for r in results:

	# The list should include only publications with PubMed identifiers. (per TR)

	if pubmed.has_key(r['_Marker_key']):

		fp.write(r['accID'] + TAB + \
	         	r['symbol'] + TAB + \
		 	r['name'] + TAB)

		if syn.has_key(r['_Marker_key']):
			fp.write(string.joinfields(syn[r['_Marker_key']], '|'))

		fp.write(TAB)
		fp.write(string.joinfields(pubmed[r['_Marker_key']], '|') + CRT)

reportlib.finish_nonps(fp)	# non-postscript file
db.useOneConnection(0)

