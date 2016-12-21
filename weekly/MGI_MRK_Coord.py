#!/usr/local/bin/python

'''
#
# Report:
#       Tab-delimited file of MGI Genome Coordinates
#	(TR 6573)
#
#   1. MGI accession id
#   2. marker type
#   3. feature type
#   4. symbol
#   5. name
#   6. chromosome
#   7. gene start
#   8. gene end
#   9. gene strand
#  10. provider
#  11. display name for provider
#
# Usage:
#       MGI_MRK_Coord.py
#
# History:
#
# 9/30/2015	lec
#	- removed un-used variables
#
# 10/09/2013 - lec - TR11502/add _Organims_key to version query
#
# 07/26/2012 - jsb - initial addition
#
'''
 
import sys
import os
import mgi_utils
import reportlib
import db

db.setTrace()

CRT = reportlib.CRT
TAB = reportlib.TAB

noneDisplay = 'null'

def getFeatureTypeCache():
    # get the feature type for each marker

    results = db.sql ('''select distinct _Marker_key, directTerms
    	from MRK_MCV_Cache''', 'auto')

    features = {}		# features[key] = [ type 1, ..., type n ]

    for row in results:
	    key = row['_Marker_key']
	    if not features.has_key(key):
		    features[key] = [ row['directTerms'] ]
	    else:
		    features[key].append (row['directTerms'])
    return features

def getProviderCache():
    # get the provider of the coordinates for each marker (not always the same
    # as the display name)

    # sequence-based coordinates first
    # and mlc.strand = scc.strand

    results = db.sql ('''select mlc._Marker_key, mcc.name
	from MRK_Location_Cache mlc,
		SEQ_Marker_Cache smc,
		SEQ_Coord_Cache scc,
		MAP_Coordinate mc,
		MAP_Coord_Collection mcc
	where mlc._Organism_key = 1
		and mlc._Marker_key = smc._Marker_key
		and smc._Qualifier_key = 615419
		and smc._Sequence_key = scc._Sequence_key
		and mlc.startCoordinate = scc.startCoordinate
		and mlc.endCoordinate = scc.endCoordinate
		and mlc.provider = mcc.abbreviation
		and scc._Map_key = mc._Map_key
		and mc._Collection_key = mcc._Collection_key''', 'auto')

    providers = {}	# providers[marker key] = provider name

    for row in results:
	    providers[row['_Marker_key']] = row['name']

    # then coordinates directly attached to markers
    # and mlc.strand = mcf.strand

    results = db.sql ('''select mlc._Marker_key, mcc.name
    	from MRK_Location_Cache mlc,
		MAP_Coord_Feature mcf,
		MAP_Coordinate mc,
		MAP_Coord_Collection mcc
	where mlc._Organism_key = 1
		and mlc.startCoordinate = mcf.startCoordinate
		and mlc.endCoordinate = mcf.endCoordinate
		and mlc._Marker_key = mcf._Object_key
		and mcf._MGIType_key = 2
		and mlc.provider = mcc.abbreviation
		and mcf._Map_key = mc._Map_key
		and mc._Collection_key = mcc._Collection_key''', 'auto')

    for row in results:
	    providers[row['_Marker_key']] = row['name']
    return providers

def getCoordinates():
    # get the coordinates for each marker

    results = db.sql ('''select mlc._Marker_key,
    		mt.name as markerType,
		mm.symbol,
		mm.name,
		mlc.chromosome,
		mlc.genomicChromosome,
		mlc.startCoordinate::int as startCoordinate,
		mlc.endCoordinate::int as endCoordinate,
		mlc.strand,
		mlc.provider as displayName,
		aa.accID
	from MRK_Location_Cache mlc,
		MRK_Marker mm,
		MRK_Types mt,
		ACC_Accession aa
	where mlc._Marker_key = mm._Marker_key
		and mm._Marker_Type_key = mt._Marker_Type_key
		and mlc._Marker_key = aa._Object_key
		and aa._MGIType_key = 2
		and aa.private = 0
		and aa.preferred = 1
		and aa._LogicalDB_key = 1
		and aa.prefixPart = 'MGI:'
		and mlc._Organism_key = 1
		and mlc.startCoordinate is not null''', 'auto')

    return results

#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('1. MGI Marker Accession ID' + TAB)
fp.write('2. Marker Type' + TAB)
fp.write('3. Feature Type' + TAB)
fp.write('4. Marker Symbol' + TAB)
fp.write('5. Marker Name' + TAB)
fp.write('6. Chromosome' + TAB)
fp.write('7. Start Coordinate' + TAB)
fp.write('8. End Coordinate' + TAB)
fp.write('9. Strand' + TAB)
fp.write('10. Genome Build' + TAB)
fp.write('11. Provider Collection' + TAB)
fp.write('12. Provider Display' + CRT)

# load data into memory

coords = getCoordinates()
providers = getProviderCache()
featureTypes = getFeatureTypeCache()

results = db.sql ('''select distinct version
	from MRK_Location_Cache
	where _Organism_key = 1 and version is not null''', 'auto')

if results:
	genomeBuild = results[0]['version']
else:
	genomeBuild = 'Unknown'

# walk through the markers to build our output file

for r in coords:
    key = r['_Marker_key']

    if featureTypes.has_key(key):
	    featureType = ';'.join(featureTypes[key])
    else:
	    featureType = noneDisplay

    if providers.has_key(key):
	    provider = providers[key]
    else:
	    provider = noneDisplay

    fp.write(r['accID'] + TAB)
    fp.write(r['markerType'] + TAB)
    fp.write(featureType + TAB)
    fp.write(r['symbol'] + TAB)
    fp.write(r['name'] + TAB)

    # prefer to display genomic chromosome (associated with coordinates) 
    # rather than genetic chromosome (associated with cM / cytoband)
    if r['genomicChromosome']:
	fp.write(r['genomicChromosome'] + TAB)
    else:
	fp.write(r['chromosome'] + TAB)

    fp.write(str(r['startCoordinate']) + TAB)
    fp.write(str(r['endCoordinate']) + TAB)
    fp.write(mgi_utils.prvalue(r['strand']) + TAB)
    fp.write(genomeBuild + TAB)
    fp.write(provider + TAB)
    fp.write(r['displayName'] + TAB)

    fp.write(CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
