#!/usr/local/bin/python

'''
# Name: MGI_OMIM.py
#
# (formerly:  MGI_PhenoOMIM.py  "TR9308: OMIM report for Jeremy Miller")
#
# Lists mouse markers and human markers with disease annotations, plus the
# homology concordance between them.  The latest modifications stem from the
# scrum-bob N-to-M orthology (n2mO) project, User Story 107, Sprint 7.
# 
# Report columns:
#
# 1) OMIM Disease ID
# 2) OMIM Disease Name
# 3) HomoloGene ID (for the homology class of the marker, if one)
# 4) Common Organism Name
# 5) NCBI Taxon ID (for the organism)
# 6) Symbol (for the marker)
# 7) EntrezGene ID
# 8) Mouse MGI ID (mouse markers only)
#
# Sort by: 2, 3, 5 (above)
#
# History:
#
# 04/17/2013 - jsb - rewrote for n-to-m homologies
#
# 07/05/2011	lec
#	- TR10770/fixed bugs
#
# 05/12/2011	lec
#	- TR9308/add MGI ids/add to public reports
#
# 2008-11-06	Susan McClatchy
#	- created
'''

import sys
import os
import string
import reportlib
import symbolsort

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
	db.setAutoTranslate(False)
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db


CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

OMIM_GENOTYPE = 1005		# from VOC_AnnotType
NOT_QUALIFIER = 1614157		# from VOC_Term
TERM_MGITYPE = 13		# from ACC_MGIType
DRIVER_NOTE = 1034		# from MGI_NoteType
GT_ROSA = 37270			# MRK_Marker for 'Gt(ROSA)26Sor'
MOUSE = 1			# MGI_Organism

###--- Functions ---###

def resolve (key, keys):
	# Postgres can return fieldnames in a different case, so we need this
	# method to map from a case sensitive key name to a non-case sensitive
	# key name

	if key in keys:			# as-is?
		return key

	if key.lower() in keys:		# lowercase?
		return key.lower()

	keymap = {}			# different mixed case?
	for k in keys:
		keymap[k.lower()] = k

	if keymap.has_key(key.lower()):
		return keymap[key.lower()]

	raise 'error', 'Cannot find key (%s) in keys %s' % (key, str(keys))

def buildCache (results, keyName, valueName):
	cache = {}

	if results:
		keys = results[0].keys()

		keyName = resolve(keyName, keys)
		valueName = resolve(valueName, keys)

		for row in results:
			cache[row[keyName]] = row[valueName]
	return cache

def lookup (cache, key):
	if not cache.has_key(key):
		return None
	return cache[key]

def customSort (a, b):
	# custom sort function for global 'results' list of tuples

	# smart alpha sort on disease term
	termCmp = symbolsort.nomenCompare(a[0], b[0])
	if termCmp != 0:
		return termCmp

	# smart alpha sort on HomoloGene ID
	if a[1] == None:
		if b[1] != None:
			return 1
	elif b[1] == None:
		return -1
	else:
		hgCmp = symbolsort.nomenCompare(a[1], b[1])
		if hgCmp != 0:
			return hgCmp

	# smart alpha sort on NCBI Taxon ID
	taxCmp = symbolsort.nomenCompare(a[2], b[2])
	if taxCmp != 0:
		return taxCmp

	# shouldn't get to this point, but if we do, just sort by the tuple
	# as a whole to have a deterministic fourth criteria

	return cmp(a, b)

def maskNulls (s):
	if s == None:
		return ''
	return str(s)

###--- Main Program ---###

db.useOneConnection(1)

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('OMIM Disease ID' + TAB)
fp.write('OMIM Disease Name' + TAB)
fp.write('HomoloGene ID' + TAB)
fp.write('Common Organism Name' + TAB)
fp.write('NCBI Taxon ID' + TAB)
fp.write('Symbol' + TAB)
fp.write('EntrezGene ID' + TAB)
fp.write('Mouse MGI ID' + CRT)

#
# cache organism names and taxon IDs
#
organismResults = db.sql('''select o._Organism_key, o.commonName, a.accID
	from MGI_Organism o
	left outer join ACC_Accession a on (
		o._Organism_key = a._Object_key
		and a._MGIType_key = 20
		and a._LogicalDB_key = 32)''', 'auto')

commonNames = buildCache(organismResults, '_Organism_key', 'commonName')
taxonIDs = buildCache(organismResults, '_Organism_key', 'accID') 

#
# cache OMIM disease terms and IDs
#
termResults = db.sql('''select t._Term_key, t.term, a.accID
	from VOC_Term t, ACC_Accession a
	where t._Vocab_key = 44
		and t._Term_key = a._Object_key
		and a._MGIType_key = 13
		and a._LogicalDB_key = 15''', 'auto')

diseaseTerms = buildCache(termResults, '_Term_key', 'term')
diseaseIDs = buildCache(termResults, '_Term_key', 'accID')

#
# cache the HomoloGene ID for each marker, where available
#
homoloGeneResults = db.sql('''select mc.clusterID, mcm._Marker_key
	from MRK_Cluster mc, MRK_ClusterMember mcm, VOC_Term vt
	where mc._Cluster_key = mcm._Cluster_key
		and mc._ClusterSource_key = vt._Term_key
		and vt.term = 'HomoloGene' ''', 'auto')

homoloGeneIDs = buildCache(homoloGeneResults, '_Marker_key', 'clusterID')

#
# cache the EntrezGene ID for each marker, where available
#
entrezGeneResults = db.sql('''select m._Marker_key, a.accID
	from MRK_Marker m, ACC_Accession a
	where m._Marker_key = a._Object_key
		and a._MGIType_key = 2
		and a.preferred = 1
		and a._LogicalDB_key = 55
		and a.private = 0''', 'auto')

entrezGeneIDs = buildCache(entrezGeneResults, '_Marker_key', 'accID')

#
# cache the MGI ID for each mouse marker, where available
#
mgiResults = db.sql('''select m._Marker_key, a.accID
	from MRK_Marker m, ACC_Accession a
	where m._Marker_key = a._Object_key
		and a._MGIType_key = 2
		and a.preferred = 1
		and a._LogicalDB_key = 1
		and a.private = 0''', 'auto')

mgiIDs = buildCache(mgiResults, '_Marker_key', 'accID')

#
# cache the symbol for each marker, where available
#
symbolResults = db.sql('''select m._Marker_key, m.symbol
	from MRK_Marker m
	where m._Marker_Status_key in (1,3)''', 'auto')

symbols = buildCache(symbolResults, '_Marker_key', 'symbol')

# 
# get the disease associations for human markers:
#	1. these have no NOT annotations
#	2. these are annotated directly to markers
#
humanResults = db.sql('''select distinct va._Object_key as _Marker_key,
		va._Term_key,
		mm._Organism_key
	from VOC_Annot va, MRK_Marker mm
	where va._AnnotType_key = 1006
		and mm._Marker_key = va._Object_key''', 'auto') 

# each item in 'results' consists of:
# (disease name, HomoloGene ID, Taxon ID, term key, marker key, organism key)

results = []

for row in humanResults:
	termKey = row['_Term_key']
	markerKey = row['_Marker_key']
	organismKey = row['_Organism_key']

	results.append ( (lookup(diseaseTerms, termKey),
		lookup(homoloGeneIDs, markerKey),
		lookup(taxonIDs, organismKey),
		termKey, 
		markerKey,
		organismKey) )

#
# get the disease associations for mouse markers.  These annotations are
# pulled up from genotype to allele to marker, wherever there is a valid path
# from the marker to the genotype.  We exclude paths involving:
#	1. recombinase alleles (ones with driver notes)
#	2. wild-type alleles
#	3. complex, not conditional genotypes
#	4. complex, not conditional genotypes with transgenes
#	5. marker Gt(ROSA)
#	6. annotations with NOT qualifiers
#

# First we need to identify the complex-not-conditional genotypes.  To do this
# we work through a series of queries...

allQuery = '''select _Genotype_key from gxd_genotype'''

homozygoteQuery = '''select _Genotype_key
	from gxd_allelepair
	where _Compound_key = 847167
		and _PairState_key = 847138'''
		
heterozygoteQuery = '''select _Genotype_key
	from gxd_allelepair
	where _Compound_key = 847167
		and _PairState_key = 847137'''

transgeneQuery = '''select distinct g._Genotype_key
	from gxd_allelegenotype g, all_allele a
	where g._Allele_key = a._Allele_key
		and a._Allele_Type_key in (847127, 847128, 847129, 2327160)'''

complexQuery = '''select _Genotype_key, count(1)
	from gxd_allelepair
	group by _Genotype_key
	having count(1) > 1'''

conditionalQuery = '''select _Genotype_key
	from gxd_genotype
	where isConditional = 1'''

# queries and their respective abbreviations; later queries take
# precedence when setting the type for a genotype
queries = [ (allQuery, 'ot'),
		(homozygoteQuery, 'hm'),
		(heterozygoteQuery, 'ht'),
		(transgeneQuery, 'tg'),
		(complexQuery, 'cx'),
		(conditionalQuery, 'cn') ]

genotypeTypes = {}		# genotype key -> two letter type abbreviation

for (query, abbrev) in queries:
	rows = db.sql(query, 'auto')

	for row in rows:
		genotypeTypes[row['_Genotype_key']] = abbrev

# Now having identified the type of each genotype, we can proceed with pulling
# OMIM annotations up to markers, remembering the six rules above.

cmd = '''select distinct gag._Marker_key,
		gag._Genotype_key,
		va._Term_key
	from gxd_genotype gg,
		gxd_allelegenotype gag,
		voc_annot va,
		all_allele a
	where gg._Genotype_key = gag._Genotype_key
		and gg._Genotype_key = va._Object_key
		and va._AnnotType_key = %d
		and va._Qualifier_key != %d
		and gag._Allele_key = a._Allele_key
		and a.isWildType = 0
		and not exists (select 1 from MGI_Note mn
			where mn._NoteType_key = %d
			and mn._Object_key = gag._Allele_key)
		and gag._Marker_key != %d''' % (
	OMIM_GENOTYPE, NOT_QUALIFIER, DRIVER_NOTE, GT_ROSA)

rows = db.sql(cmd, 'auto')


# remember -- each item in 'results' consists of:
# (disease name, HomoloGene ID, Taxon ID, term key, marker key, organism key)

alreadySeen = {}		# (marker key, term key) = 1

for row in rows:
	termKey = row['_Term_key']
	markerKey = row['_Marker_key']
	genotypeKey = row['_Genotype_key']

	# skip complex-not-conditional genotypes
	if genotypeTypes[genotypeKey] == 'cx':
		continue

	pair = (markerKey, termKey)
	if alreadySeen.has_key(pair):
		continue
	alreadySeen[pair] = 1

	results.append ( (lookup(diseaseTerms, termKey),
		lookup(homoloGeneIDs, markerKey),
		lookup(taxonIDs, MOUSE),
		termKey, 
		markerKey,
		MOUSE) )

results.sort(customSort)

# Output fields:
# 1) OMIM Disease ID
# 2) OMIM Disease Name
# 3) HomoloGene ID (for the homology class of the marker, if one)
# 4) Common Organism Name
# 5) NCBI Taxon ID (for the organism)
# 6) Symbol (for the marker)
# 7) EntrezGene ID
# 8) Mouse MGI ID (mouse markers only)
#
# Sorted by: 2, 3, 5 (above)

for (disease, hgID, taxonID, termKey, markerKey, organismKey) in results:
	diseaseID = lookup(diseaseIDs, termKey)
	organism = lookup(commonNames, organismKey)
	symbol = lookup(symbols, markerKey)
	egID = lookup(entrezGeneIDs, markerKey)
	mgiID = lookup(mgiIDs, markerKey)

	for item in [ diseaseID, disease, hgID, organism, taxonID, symbol,
		egID, mgiID ]:
			fp.write (maskNulls(item) + TAB)
	fp.write(CRT)

db.useOneConnection(0)
reportlib.finish_nonps(fp)

