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
#	15. VEGA ID
#	16. VEGA Chr
#	17. VEGA Start Coord
#	18. VEGA End Coord
#	19. VEGA Strand
#	20. Mouse GenBank Ids
#	21. Mouse UniGene Ids
#	22. Mouse RefSeq Ids
#	23. Ensembl Transcript Ids
#	24. Ensembl Protein Ids
#	25. VEGA Transcript Ids
#	26. VEGA Protein Ids
#	27. Mouse SwissProt Ids
#	28. Mouse InterPro Ids
#	29. Mouse Synonyms
#	30. Human EG ID
#	31. Human Symbol
#	32. Human Name
#	33. Human Chr
#	34. Human RefSeq Ids
#	35. Human Synonyms
#
# Usage:
#       MGI_MouseHumanSequence.py
#
# History:
#
# sc    03/12/2010
#       - TR9774 Add Ensembl and VEGA transcripts
#
# lec	07/14/2005
#	- TR 6973
#
'''
 
import sys 
import os
import string
import mgi_utils
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

coordDisplay = '(%s:%s-%s (%s))'
noneDisplay = TAB
repGenomicKey = 615419
sequenceType = 19
ncbi = 59
ensembl = 60
vega = 85
valueDelimiter = ','

# lookup associated transcripts by Ensembl or VEGA genomic ID
genomicToTranscript = {}

# lookup associated proteins by Ensembl or VEGA genomic ID
genomicToProtein = {}

#
# Load Lookups
#

def init():
    global genomicToTranscript, genomicToProtein

    # Transcript ID lookup by Genomic ID
    # transcribed from qualifier
    db.sql('''
	select sa._Sequence_key_1 as transcriptKey, 
	sa._Sequence_key_2 as genomicKey 
	into #transGen 
	from SEQ_Sequence_Assoc sa 
	where sa._Qualifier_key = 5445464
	''', None)  
    db.sql('create index transGene_idx1 on transGen(transcriptKey)', None)
    db.sql('create index transGene_idx2 on #transGen(genomicKey)', None)

    results = db.sql('''
	    select a1.accID as genomicID, a2.accID as transcriptID 
	    from  #transGen t, ACC_Accession a1, ACC_Accession a2 
	    where t.genomicKey = a1._Object_key 
	    and a1._MGIType_key = 19 
	    and a1.preferred = 1 
	    and t.transcriptKey = a2._Object_key 
	    and a2._MGIType_key = 19 
	    and a2.preferred = 1 
	    order by a1.accID
	    ''', 'auto')
    genomicToTranscript = {}
    for r in results:
	#print 'r:%s' % r
	key = r['genomicID']
	value = r['transcriptID']
	if not genomicToTranscript.has_key(key):
	    genomicToTranscript[key] = []
	genomicToTranscript[key].append(value)

    # Protein ID lookup by Genomic ID
    db.sql('''
	select tg.genomicKey, 
	sa._Sequence_key_1 as proteinKey 
	into #protGen 
	from #transGen tg, SEQ_Sequence_Assoc sa 
	where sa._Qualifier_key = 5445465 
	and tg.transcriptKey =  sa._Sequence_key_2
	''', None)

    results = db.sql('''
	    select a1.accID as genomicID, a2.accID as proteinID 
	    from  #protGen t, ACC_Accession a1, ACC_Accession a2 
	    where t.genomicKey = a1._Object_key 
	    and a1._MGIType_key = 19 
	    and a1.preferred = 1 
	    and t.proteinKey = a2._Object_key 
	    and a2._MGIType_key = 19 
	    and a2.preferred = 1 
	    order by a1.accID
	    ''', 'auto')
    genomicToProtein = {}
    for r in results:
	key = r['genomicID']
	value = r['proteinID']
	if not genomicToProtein.has_key(key):
	    genomicToProtein[key] = []
	genomicToProtein[key].append(value)

#
# get coordinates
#

def getCoords(logicalDBkey):

    tempCoords = {}

    results = db.sql('''
	    select m._Marker_key, mc._Qualifier_key, mc.accID, 
	    c.chromosome, c.strand, 
	    convert(int, c.startCoordinate) as startC, 
	    convert(int, c.endCoordinate) as endC 
	        from #repmarkers m, SEQ_Marker_Cache mc, SEQ_Coord_Cache c 
	        where m._Marker_key = mc._Marker_key 
	        and mc._Sequence_key = c._Sequence_key 
	        and mc._LogicalDB_key = %d
		''' % (logicalDBkey) , 'auto')
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
fp = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])
init()

#
# select all mouse markers
#

db.sql('''
	select m._Marker_key, m.symbol, m.name 
	into #markers 
	from MRK_Marker m 
	where m._Organism_key = 1 and m._Marker_Status_key in (1,3)
	''', None)
db.sql('create index markers_idx1 on #markers(_Marker_key)', None)

#
# select markers that have genomic coordinates
#

db.sql('''
	select m._Marker_key, m._Sequence_key, c.version into #repmarkers 
	from SEQ_Marker_Cache m, SEQ_Coord_Cache c 
	where m._Qualifier_key = %d 
	and m._Sequence_key = c._Sequence_key
	''' % (repGenomicKey), None)
db.sql('create index repmarkers_idx1 on #repmarkers(_Marker_key)', None)
db.sql('create index repmarkers_idx2 on #repmarkers(_Sequence_key)', None)

# NCBI, Ensembl and VEGA coordinates

ncbiCoords = getCoords(ncbi)
ensemblCoords = getCoords(ensembl)
vegaCoords = getCoords(vega)

#
# select markers that have human orthologs
#

db.sql('''
	select distinct m._Marker_key as mouseKey, 
	h2._Marker_key as humanKey, 
	m2.symbol as humanSym, 
	m2.name as humanName, 
	m2.chromosome || m2.cytogeneticOffset as humanChr
	into #homology 
        from #markers m, MRK_Homology_Cache h1, MRK_Homology_Cache h2, MRK_Marker m2 
        where m._Marker_key = h1._Marker_key 
        and h1._Class_key = h2._Class_key 
        and h2._Organism_key = 2 
        and h2._Marker_key = m2._Marker_key 
	''', None)
db.sql('create index index_mouseKey on #homology(mouseKey)', None)
db.sql('create index index_humanKey on #homology(humanKey)', None)
db.sql('create index index_humanSym on #homology(humanSym)', None)

# cm for Mouse
results = db.sql('''
	select o._Marker_key, o.offset 
	from #markers m, MRK_Offset o 
	where m._Marker_key = o._Marker_key 
	and o.source = 0 
	''', 'auto') 
mCM = {}
for r in results:
	mCM[r['_Marker_key']] = r['offset']

# MGI for Mouse
results = db.sql('''
	select a._Object_key, a.accID 
	from #markers m, ACC_Accession a 
	where m._Marker_key = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 1 
	and a.prefixPart = "MGI:" 
	and a.preferred = 1
	''', 'auto')
mgiID = {}
for r in results:
	mgiID[r['_Object_key']] = r['accID']

# EntrezGene for Mouse is included in NCBI Coordinates

# EntrezGene for Human
results = db.sql('''
	select distinct a._Object_key, a.accID 
	from #homology h, ACC_Accession a 
	where h.humanKey = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 55 
	''', 'auto')
hegID = {}
for r in results:
	hegID[r['_Object_key']] = r['accID']

# RefSeq for Mouse
results = db.sql('''
	select distinct a._Object_key, a.accID 
	from #markers m, ACC_Accession a 
	where m._Marker_key = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 27 
	''', 'auto')
mrefseqID = {}
for r in results:
	if not mrefseqID.has_key(r['_Object_key']):
		mrefseqID[r['_Object_key']] = []
	mrefseqID[r['_Object_key']].append(r['accID'])

# RefSeq for Human
results = db.sql('''
	select distinct h.humanKey as _Object_key, a.accID 
	from #homology h, ACC_Accession a 
	where h.humanKey = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 27 
	''', 'auto')
hrefseqID = {}
for r in results:
	if not hrefseqID.has_key(r['_Object_key']):
		hrefseqID[r['_Object_key']] = []
	hrefseqID[r['_Object_key']].append(r['accID'])

# SWISSPROT for Mouse
results = db.sql('''
	select distinct a._Object_key, a.accID 
	from #markers m, ACC_Accession a 
	where m._Marker_key = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key in (13, 41) 
	''', 'auto')
mspID = {}
for r in results:
	if not mspID.has_key(r['_Object_key']):
		mspID[r['_Object_key']] = []
	mspID[r['_Object_key']].append(r['accID'])

# GenBank for Mouse
results = db.sql('''
	select distinct a._Object_key, a.accID 
	from #markers m, ACC_Accession a 
	where m._Marker_key = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 9 
	''', 'auto')
gbID = {}
for r in results:
	if not gbID.has_key(r['_Object_key']):
		gbID[r['_Object_key']] = []
	gbID[r['_Object_key']].append(r['accID'])

# UniGene for Mouse
results = db.sql('''
	select distinct a._Object_key, a.accID 
	from #markers m, ACC_Accession a 
	where m._Marker_key = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 23 
	''', 'auto')
ugID = {}
for r in results:
	if not ugID.has_key(r['_Object_key']):
		ugID[r['_Object_key']] = []
	ugID[r['_Object_key']].append(r['accID'])

# InterPro for Mouse
results = db.sql('''
	select distinct a._Object_key, a.accID 
	from #markers m, ACC_Accession a 
	where m._Marker_key = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 28 
	''', 'auto')
ipID = {}
for r in results:
	if not ipID.has_key(r['_Object_key']):
		ipID[r['_Object_key']] = []
	ipID[r['_Object_key']].append(r['accID'])

# synonyms for Mouse

results = db.sql('''
	select distinct m._Marker_key, s.synonym 
	from #markers m, MGI_Synonym s  
	where m._Marker_key = s._Object_key 
	and s._MGIType_key = 2 
	''', 'auto')
mSyn = {}
for r in results:
	if not mSyn.has_key(r['_Marker_key']):
		mSyn[r['_Marker_key']] = []
	mSyn[r['_Marker_key']].append(r['synonym'])

# synonyms for Human

results = db.sql('''
	select distinct h.mouseKey, s.synonym 
	from #homology h, MGI_Synonym s 
	where h.humanKey = s._Object_key 
	and s._MGIType_key = 2 
	''', 'auto')
hSyn = {}
for r in results:
	if not hSyn.has_key(r['mouseKey']):
		hSyn[r['mouseKey']] = []
	hSyn[r['mouseKey']].append(r['synonym'])

#
# write results
#

results = db.sql('''
	(
	select m._Marker_key, m.symbol, m.name, h.humanKey, h.humanSym, h.humanName, h.humanChr 
	from #markers m, #homology h 
	where m._Marker_key = h.mouseKey 
	union 
	select m._Marker_key, m.symbol, m.name, null, null, null, null 
	from #markers m 
	where not exists (select 1 from #homology h 
	where m._Marker_key = h.mouseKey) 
	)
	order by symbol
	''', 'auto')

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
	ensemblID = '' # we need this later
        if ensemblCoords.has_key(key):
	    c = ensemblCoords[key]
	    ensemblID = c['accID'] 
	    fp.write(c['accID'] + TAB)
            fp.write(c['chromosome'] + TAB)
            fp.write(str(c['startC']) + TAB)
            fp.write(str(c['endC']) + TAB)
            fp.write(c['strand'] + TAB)
        else:
	    fp.write(5*noneDisplay)

#	15. VEGA ID
#	16. VEGA Chr
#	17. VEGA Start Coord
#	18. VEGA End Coord
#	19. VEGA Strand
	vegaID = '' # we need this later
        if vegaCoords.has_key(key):
	    c = vegaCoords[key]
	    vegaID = c['accID']
	    fp.write(vegaID + TAB)
            fp.write(c['chromosome'] + TAB)
            fp.write(str(c['startC']) + TAB)
            fp.write(str(c['endC']) + TAB)
            fp.write(c['strand'] + TAB)
        else:
	    fp.write(5*noneDisplay)

#	20. Mouse GenBank Ids

	if gbID.has_key(key):
		fp.write(string.join(gbID[key], valueDelimiter))
	fp.write(TAB)

#	21. Mouse UniGene Ids

	if ugID.has_key(key):
		fp.write(string.join(ugID[key], valueDelimiter))
	fp.write(TAB)

#	22. Mouse RefSeq Ids

	if mrefseqID.has_key(key):
		fp.write(string.join(mrefseqID[key], valueDelimiter))
	fp.write(TAB)
#	23. Ensembl Transcript IDs
#	24. Ensembl Protein IDs
	if ensemblID != '':
	    if genomicToTranscript.has_key(ensemblID):
		fp.write(string.join(genomicToTranscript[ensemblID], valueDelimiter))
	    fp.write(TAB)
	    if genomicToProtein.has_key(ensemblID):
		fp.write(string.join(genomicToProtein[ensemblID], valueDelimiter))
	    fp.write(TAB)
	else:
	     fp.write(2*noneDisplay)

#       25. VEGA Transcript IDs
#       26. VEGA Protein IDs
        if vegaID != '':
            if genomicToTranscript.has_key(vegaID):
                fp.write(string.join(genomicToTranscript[vegaID], valueDelimiter))
            fp.write(TAB)
            if genomicToProtein.has_key(vegaID):
                fp.write(string.join(genomicToProtein[vegaID], valueDelimiter))
            fp.write(TAB)
        else:
             fp.write(2*noneDisplay)

#	27. Mouse SwissProt Ids

	if mspID.has_key(key):
		fp.write(string.join(mspID[key], valueDelimiter))
	fp.write(TAB)

#	28. Mouse InterProt Ids

	if ipID.has_key(key):
		fp.write(string.join(ipID[key], valueDelimiter))
	fp.write(TAB)

#	29. Mouse Synonyms

	if mSyn.has_key(key):
		fp.write(string.join(mSyn[key], valueDelimiter))
	fp.write(TAB)

#	30. Human EG ID

	if r['humanKey'] != None:
		if hegID.has_key(r['humanKey']):
			fp.write(hegID[r['humanKey']])
		fp.write(TAB)

#		31. Human Symbol
#		32. Human Name
#		33. Human Chr

		fp.write(r['humanSym'] + TAB)
		fp.write(r['humanName'] + TAB)
		fp.write(mgi_utils.prvalue(r['humanChr']) + TAB)

#		34. Human RefSeq Ids

		if hrefseqID.has_key(r['humanKey']):
			fp.write(string.join(hrefseqID[r['humanKey']], valueDelimiter))
		fp.write(TAB)

#		35. Human Synonyms

		if hSyn.has_key(key):
			fp.write(string.join(hSyn[key], valueDelimiter))

	else:
		fp.write(5*TAB)

	fp.write(CRT)

reportlib.finish_nonps(fp)	# non-postscript file
db.useOneConnection(0)
