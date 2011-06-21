#!/usr/local/bin/python

'''
#
# TR 3655
#
# Report:
#       Tab-delimited file
#
#	1. mouse gene symbol
#	2. mgi acc id
#	3. entrezgene acc id for mouse gene
#	4. gene symbol for human ortholog
#	5. entrezgene acc id for human gene
#	6. nucleotide refseq acc id for mouse gene
#	7. nucleotide refseq acc id for human gene
#	8. protein refseq acc id for mouse gene
#	9. protein refseq acc id for human gene
#	10. sp acc id for mouse gene
#	11. sp acc id for human gene
#	12. if no ref seq or sp id, then GB acc ids for mouse gene
#	13. evidence used to support mouse0human homology
#	14. J-Number(s) for references supporting the orthology
#	15. PubMed ID(s) for references supporting the orthology
#
# Usage:
#       HMD_HumanSequence.py
#
# History:
#
# dbm	04/26/2007
#	- TR 8277; add J-Numbers and PubMed IDs
#
# lec	10/25/2005
#	- MGI 3.5; now loading Human RefSeqs directly into MGI
#
# lec	01/04/2004
#	- TR 5939; EntrezGene->EntrezGene
#
# lec	07/01/2003
#	- TR 4945; added EntrezGene IDs for Mouse and Human
#
# lec	05/07/2002
#	- created
#
'''
 
import sys 
import os
import string
import db
import reportlib

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])

cmds = []

cmds.append('select distinct mouseKey = h1._Marker_key, mouseSym = m1.symbol, ' +
	'humanKey = h2._Marker_key, humanSym = m2.symbol, a.abbrev ' + \
	'into #homology ' +
        'from MRK_Homology_Cache h1, HMD_Homology_Assay ha, HMD_Assay a, ' + \
        'MRK_Homology_Cache h2, ' + \
        'MRK_Marker m1, MRK_Marker m2 ' + \
        'where h1._Organism_key = 1 ' + \
        'and h1._Class_key = h2._Class_key ' + \
        'and h2._Organism_key = 2 ' + \
        'and h1._Marker_key = m1._Marker_key ' + \
        'and h2._Marker_key = m2._Marker_key ' + \
	'and h1._Homology_key = ha._Homology_key ' + \
	'and ha._Assay_key = a._Assay_key')

cmds.append('create nonclustered index index_mouseKey on #homology(mouseKey)')
cmds.append('create nonclustered index index_humanKey on #homology(humanKey)')
cmds.append('create nonclustered index index_humanSym on #homology(humanSym)')

cmds.append('select distinct h.mouseKey, mr.jnumID, mr.pubMedID ' + \
	'into #homologyRef ' + \
	'from #homology h, MRK_Homology_Cache hc, MRK_Reference mr ' + \
	'where h.mouseKey = hc._Marker_key and ' + \
      	'hc._Organism_key = 1 and ' + \
      	'hc._Marker_key = mr._Marker_key and ' + \
      	'hc._Refs_key = mr._Refs_key')

cmds.append('create nonclustered index index_mouseKey on #homologyRef(mouseKey)')

db.sql(cmds, None)

##

# MGI for Mouse
results = db.sql('select a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.mouseKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1', 'auto')
mgiID = {}
for r in results:
	mgiID[r['_Object_key']] = r['accID']

# EntrezGene for Mouse
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.mouseKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 55 ', 'auto')
megID = {}
for r in results:
	megID[r['_Object_key']] = r['accID']


# EntrezGene for Human
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.humanKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 55 ', 'auto')
hegID = {}
for r in results:
	hegID[r['_Object_key']] = r['accID']

# nucleotide RefSeqs for Mouse
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.mouseKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 27 ' + \
	'and a.prefixPart in ("NM_", "XM_")', 'auto')
mNrefseqID = {}
for r in results:
	if not mNrefseqID.has_key(r['_Object_key']):
		mNrefseqID[r['_Object_key']] = []
	mNrefseqID[r['_Object_key']].append(r['accID'])

# nucleotide RefSeqs for Human
results = db.sql('select distinct _Object_key = h.humanKey, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.humanKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 27 ' + \
	'and a.prefixPart in ("NM_", "XM_")', 'auto')
hNrefseqID = {}
for r in results:
	if not hNrefseqID.has_key(r['_Object_key']):
		hNrefseqID[r['_Object_key']] = []
	hNrefseqID[r['_Object_key']].append(r['accID'])

# protein RefSeqs for Mouse
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.mouseKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 27 ' + \
	'and a.prefixPart in ("NP_", "XP_")', 'auto')
mPrefseqID = {}
for r in results:
	if not mPrefseqID.has_key(r['_Object_key']):
		mPrefseqID[r['_Object_key']] = []
	mPrefseqID[r['_Object_key']].append(r['accID'])

# protein RefSeqs for Human
results = db.sql('select distinct _Object_key = h.humanKey, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.humanKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 27 ' + \
	'and a.prefixPart in ("NP_", "XP_")', 'auto')
hPrefseqID = {}
for r in results:
	if not hPrefseqID.has_key(r['_Object_key']):
		hPrefseqID[r['_Object_key']] = []
	hPrefseqID[r['_Object_key']].append(r['accID'])

# SWISSPROT for Mouse
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.mouseKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 13', 'auto')
mspID = {}
for r in results:
	if not mspID.has_key(r['_Object_key']):
		mspID[r['_Object_key']] = []
	mspID[r['_Object_key']].append(r['accID'])

# SWISSPROT for Human
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.humanKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 13', 'auto')
hspID = {}
for r in results:
	if not hspID.has_key(r['_Object_key']):
		hspID[r['_Object_key']] = []
	hspID[r['_Object_key']].append(r['accID'])

# GenBank for Mouse
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.mouseKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 9 ', 'auto')
gbID = {}
for r in results:
	if not gbID.has_key(r['_Object_key']):
		gbID[r['_Object_key']] = []
	gbID[r['_Object_key']].append(r['accID'])

results = db.sql('select distinct mouseKey, abbrev from #homology', 'auto')
habbrev = {}
for r in results:
	if not habbrev.has_key(r['mouseKey']):
		habbrev[r['mouseKey']] = []
	habbrev[r['mouseKey']].append(r['abbrev'])

# J-Numbers
jnumID = {}
results = db.sql('select distinct mouseKey, jnumID from #homologyRef', 'auto')
for r in results:
	if not jnumID.has_key(r['mouseKey']):
		jnumID[r['mouseKey']] = []
	jnumID[r['mouseKey']].append(r['jnumID'])

# Pubmed IDs
pubMedID = {}
results = db.sql('select distinct mouseKey, pubMedID ' + \
	'from #homologyRef ' + \
	'where pubMedID is not null', 'auto')
for r in results:
	if not pubMedID.has_key(r['mouseKey']):
		pubMedID[r['mouseKey']] = []
	pubMedID[r['mouseKey']].append(r['pubMedID'])

results = db.sql('select distinct mouseKey, mouseSym, humanKey, humanSym from #homology', 'auto')

for r in results:

	fp.write(r['mouseSym'] + TAB)
	fp.write(mgiID[r['mouseKey']] + TAB)

	if megID.has_key(r['mouseKey']):
		fp.write(megID[r['mouseKey']])
	fp.write(TAB)

	fp.write(r['humanSym'] + TAB)

	if hegID.has_key(r['humanKey']):
		fp.write(hegID[r['humanKey']])
	fp.write(TAB)

	if mNrefseqID.has_key(r['mouseKey']):
		fp.write(string.join(mNrefseqID[r['mouseKey']], ','))
	fp.write(TAB)

	if hNrefseqID.has_key(r['humanKey']):
		fp.write(string.join(hNrefseqID[r['humanKey']], ','))
	fp.write(TAB)

	if mPrefseqID.has_key(r['mouseKey']):
		fp.write(string.join(mPrefseqID[r['mouseKey']], ','))
	fp.write(TAB)

	if hPrefseqID.has_key(r['humanKey']):
		fp.write(string.join(hPrefseqID[r['humanKey']], ','))
	fp.write(TAB)

	if mspID.has_key(r['mouseKey']):
		fp.write(string.join(mspID[r['mouseKey']], ','))
	fp.write(TAB)

	if hspID.has_key(r['humanKey']):
		fp.write(string.join(hspID[r['humanKey']], ','))
	fp.write(TAB)

	if not mNrefseqID.has_key(r['mouseKey']) \
	   and not mPrefseqID.has_key(r['mouseKey']) \
	   and not mspID.has_key(r['mouseKey']) \
	   and gbID.has_key(r['mouseKey']):
		fp.write(string.join(gbID[r['mouseKey']], ','))
	fp.write(TAB)

	if habbrev.has_key(r['mouseKey']):
		fp.write(string.join(habbrev[r['mouseKey']], ','))
	fp.write(TAB)

	if jnumID.has_key(r['mouseKey']):
		fp.write(string.join(jnumID[r['mouseKey']], ','))
	fp.write(TAB)

	if pubMedID.has_key(r['mouseKey']):
		fp.write(string.join(pubMedID[r['mouseKey']], ','))
	fp.write(CRT)

reportlib.finish_nonps(fp)	# non-postscript file
db.useOneConnection(0)