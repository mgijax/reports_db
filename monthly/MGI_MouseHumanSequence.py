#!/usr/local/bin/python

'''
#
# TR 6973
#
# Report:
#
#	All Mouse Official/Unofficial Markers
#	Fields 1-4 are required
#	All other fields may be blank
#
#       Tab-delimited file
#
#	1. Mouse MGI Accession ID
#	2. Mouse Symbol
#	3. Mouse Name
#	4. Mouse cM
#	5. Mouse EG ID
#	6. NCBI Chr
#	7. NCBI Start Coord
#	8. NCBI End Coord
#	9. NCBI Strand
#	10. Mouse Ensembl ID
#	11. Ensembl Chr
#	12. Ensembl Start Coord
#	13. Ensembl End Coord
#	14. Ensembl Strand
#	15. Mouse GenBank Ids
#	16. Mouse UniGene Ids
#	17. Mouse RefSeq Ids
#	18. Mouse SwissProt Ids
#	19. Mouse InterPro Ids
#	20. Mouse Synonyms
#	21. Human EG ID
#	22. Human Symbol
#	23. Human Name
#	24. Human Chr
#	25. Human RefSeq Ids
#	26. Human Synonyms
#
# Usage:
#       MGI_MouseHumanSequence.py
#
# Notes:
#	- all reports use mgireport directory for output file
#	- all reports use db default of public login
#	- all reports use server/database default of environment
#	- use lowercase for all SQL commands (i.e. select not SELECT)
#	- all public SQL reports require the header and footer
#	- all private SQL reports require the header
#
# History:
#
# lec	07/14/2005
#	- TR 6973
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

coordDisplay = '(%s:%s-%s (%s))'
noneDisplay = TAB
repGenomicKey = 615419
sequenceType = 19
ncbi = 59
ensembl = 60
valueDelimiter = ','

#
# get coordinates
#

def getCoords(logicalDBkey):

    tempCoords = {}

    results = db.sql('select m._Marker_key, mc._Qualifier_key, mc.accID, ' + \
	    'c.chromosome, c.strand, ' + \
	    'startC = convert(int, c.startCoordinate), ' + \
	    'endC = convert(int, c.endCoordinate) ' + \
	        'from #repmarkers m, SEQ_Marker_Cache mc, SEQ_Coord_Cache c ' + \
	        'where m._Marker_key = mc._Marker_key ' + \
	        'and mc._Sequence_key = c._Sequence_key ' + \
	        'and mc._LogicalDB_key = %d ' % (logicalDBkey) , 'auto')

    for r in results:
        key = r['_Marker_key']
        value = r
        qualifier = r['_Qualifier_key']
        tempCoords[key] = value
 
    return tempCoords
    
#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], printHeading = 0, outputdir = os.environ['REPORTOUTPUTDIR'])

#
# select all mouse markers
#

db.sql('select m._Marker_key, m.symbol, m.name ' + \
	'into #markers ' + \
	'from MRK_Marker m ' + \
	'where m._Organism_key = 1 and m._Marker_Status_key in (1,3)', None)

db.sql('create index idx1 on #markers(_Marker_key)', None)

#
# select markers that have genomic coordinates
#

db.sql('select m._Marker_key, m._Sequence_key, c.version into #repmarkers ' + \
	'from SEQ_Marker_Cache m, SEQ_Coord_Cache c ' + \
	'where m._Qualifier_key = %d ' % (repGenomicKey) + \
	'and m._Sequence_key = c._Sequence_key', None)

db.sql('create index idx1 on #repmarkers(_Marker_key)', None)
db.sql('create index idx2 on #repmarkers(_Sequence_key)', None)

# NCBI and Ensembl coordinates

ncbiCoords = getCoords(ncbi)
ensemblCoords = getCoords(ensembl)

#
# select markers that have human orthologs
#

db.sql('select distinct mouseKey = m._Marker_key, ' +
	'humanKey = h2._Marker_key, humanSym = m2.symbol, humanName = m2.name, ' + \
	'humanChr = m2.chromosome + m2.cytogeneticOffset ' + \
	'into #homology ' +
        'from #markers m, HMD_Homology r1, HMD_Homology_Marker h1, ' + \
        'HMD_Homology r2, HMD_Homology_Marker h2, ' + \
        'MRK_Marker m2 ' + \
        'where m._Marker_key = h1._Marker_key ' + \
        'and h1._Homology_key = r1._Homology_key ' + \
        'and r1._Class_key = r2._Class_key ' + \
        'and r2._Homology_key = h2._Homology_key ' + \
        'and h2._Marker_key = m2._Marker_key ' + \
        'and m2._Organism_key = 2 ', None)

db.sql('create nonclustered index index_mouseKey on #homology(mouseKey)', None)
db.sql('create nonclustered index index_humanKey on #homology(humanKey)', None)
db.sql('create nonclustered index index_humanSym on #homology(humanSym)', None)

# cm for Mouse
results = db.sql('select o._Marker_key, o.offset ' + \
	'from #markers m, MRK_Offset o ' + \
	'where m._Marker_key = o._Marker_key ' + \
	'and o.source = 0 ', 'auto') 
mCM = {}
for r in results:
	mCM[r['_Marker_key']] = r['offset']

# MGI for Mouse
results = db.sql('select a._Object_key, a.accID ' + \
	'from #markers m, ACC_Accession a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1', 'auto')
mgiID = {}
for r in results:
	mgiID[r['_Object_key']] = r['accID']

# EntrezGene for Mouse is included in NCBI Coordinates

# EntrezGene for Human
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.humanKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 55 ', 'auto')
hegID = {}
for r in results:
	hegID[r['_Object_key']] = r['accID']

# RefSeq for Mouse
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #markers m, ACC_Accession a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 27 ', 'auto')
mrefseqID = {}
for r in results:
	if not mrefseqID.has_key(r['_Object_key']):
		mrefseqID[r['_Object_key']] = []
	mrefseqID[r['_Object_key']].append(r['accID'])

# RefSeq for Human
results = db.sql('select distinct _Object_key = h.humanKey, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.humanKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 27 ', 'auto')
hrefseqID = {}
for r in results:
	if not hrefseqID.has_key(r['_Object_key']):
		hrefseqID[r['_Object_key']] = []
	hrefseqID[r['_Object_key']].append(r['accID'])

# SWISSPROT for Mouse
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #markers m, ACC_Accession a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key in (13, 41) ', 'auto')
mspID = {}
for r in results:
	if not mspID.has_key(r['_Object_key']):
		mspID[r['_Object_key']] = []
	mspID[r['_Object_key']].append(r['accID'])

# GenBank for Mouse
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #markers m, ACC_Accession a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 9 ', 'auto')
gbID = {}
for r in results:
	if not gbID.has_key(r['_Object_key']):
		gbID[r['_Object_key']] = []
	gbID[r['_Object_key']].append(r['accID'])

# UniGene for Mouse
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #markers m, ACC_Accession a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 23 ', 'auto')
ugID = {}
for r in results:
	if not ugID.has_key(r['_Object_key']):
		ugID[r['_Object_key']] = []
	ugID[r['_Object_key']].append(r['accID'])

# InterPro for Mouse
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #markers m, ACC_Accession a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 28 ', 'auto')
ipID = {}
for r in results:
	if not ipID.has_key(r['_Object_key']):
		ipID[r['_Object_key']] = []
	ipID[r['_Object_key']].append(r['accID'])

# synonyms for Mouse

results = db.sql('select distinct m._Marker_key, s.synonym ' + \
	'from #markers m, MGI_Synonym s  ' + \
	'where m._Marker_key = s._Object_key ' + \
	'and s._MGIType_key = 2 ', 'auto')
mSyn = {}
for r in results:
	if not mSyn.has_key(r['_Marker_key']):
		mSyn[r['_Marker_key']] = []
	mSyn[r['_Marker_key']].append(r['synonym'])

# synonyms for Human

results = db.sql('select distinct h.mouseKey, s.synonym ' + \
	'from #homology h, MGI_Synonym s ' + \
	'where h.humanKey = s._Object_key ' + \
	'and s._MGIType_key = 2 ', 'auto')
hSyn = {}
for r in results:
	if not hSyn.has_key(r['mouseKey']):
		hSyn[r['mouseKey']] = []
	hSyn[r['mouseKey']].append(r['synonym'])

#
# write results
#

results = db.sql('select m._Marker_key, m.symbol, m.name, h.humanKey, h.humanSym, h.humanName, h.humanChr ' + \
	'from #markers m, #homology h ' + \
	'where m._Marker_key = h.mouseKey ' + \
	'union ' + \
	'select m._Marker_key, m.symbol, m.name, null, null, null, null ' + \
	'from #markers m ' + \
	'where not exists (select 1 from #homology h ' + \
	'where m._Marker_key = h.mouseKey) ' + \
	'order by m.symbol', 'auto')

for r in results:

        key = r['_Marker_key']

#	1. Mouse MGI Accession ID
#	2. Mouse Symbol
#	3. Mouse Name
#	4. Mouse cM

	fp.write(mgiID[key] + TAB)
	fp.write(r['symbol'] + TAB)
	fp.write(r['name'] + TAB)
	fp.write(str(mCM[key]) + TAB)

#	5. NCBI EG ID
#	6. NCBI Chr
#	7. NCBI Start Coord
#	8. NCBI End Coord
#	9. NCBI Strand

        if ncbiCoords.has_key(key):
	    c = ncbiCoords[key]
	    fp.write(c['accID'] + TAB)
            fp.write(c['chromosome'] + TAB)
            fp.write(str(c['startC']) + TAB)
            fp.write(str(c['endC']) + TAB)
            fp.write(c['strand'] + TAB)
        else:
	    fp.write(5*noneDisplay)

#	10. Ensembl ID
#	11. Ensembl Chr
#	12. Ensembl Start Coord
#	13. Ensembl End Coord
#	14. Ensembl Strand

        if ensemblCoords.has_key(key):
	    c = ensemblCoords[key]
	    fp.write(c['accID'] + TAB)
            fp.write(c['chromosome'] + TAB)
            fp.write(str(c['startC']) + TAB)
            fp.write(str(c['endC']) + TAB)
            fp.write(c['strand'] + TAB)
        else:
	    fp.write(5*noneDisplay)

#	15. Mouse GenBank Ids

	if gbID.has_key(key):
		fp.write(string.join(gbID[key], valueDelimiter))
	fp.write(TAB)

#	16. Mouse UniGene Ids

	if ugID.has_key(key):
		fp.write(string.join(ugID[key], valueDelimiter))
	fp.write(TAB)

#	17. Mouse RefSeq Ids

	if mrefseqID.has_key(key):
		fp.write(string.join(mrefseqID[key], valueDelimiter))
	fp.write(TAB)

#	18. Mouse SwissProt Ids

	if mspID.has_key(key):
		fp.write(string.join(mspID[key], valueDelimiter))
	fp.write(TAB)

#	19. Mouse InterProt Ids

	if ipID.has_key(key):
		fp.write(string.join(ipID[key], valueDelimiter))
	fp.write(TAB)

#	20. Mouse Synonyms

	if mSyn.has_key(key):
		fp.write(string.join(mSyn[key], valueDelimiter))
	fp.write(TAB)

#	21. Human EG ID

	if r['humanKey'] != None:
		if hegID.has_key(r['humanKey']):
			fp.write(hegID[r['humanKey']])
		fp.write(TAB)

#		22. Human Symbol
#		23. Human Name
#		24. Human Chr

		fp.write(r['humanSym'] + TAB)
		fp.write(r['humanName'] + TAB)
		fp.write(r['humanChr'] + TAB)

#		25. Human RefSeq Ids

		if hrefseqID.has_key(r['humanKey']):
			fp.write(string.join(hrefseqID[r['humanKey']], valueDelimiter))
		fp.write(TAB)

#		26. Human Synonyms

		if hSyn.has_key(key):
			fp.write(string.join(hSyn[key], valueDelimiter))

	else:
		fp.write(5*TAB)

	fp.write(CRT)

reportlib.finish_nonps(fp)	# non-postscript file
db.useOneConnection(0)
