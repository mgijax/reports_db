#!/usr/local/bin/python

'''
#
# Report:
#       Tab-delimited file of MGI Genome Coordinates
#	(TR 6573)
#
#   1. MGI accession id
#   2. marker type
#   3. symbol
#   4. name
#   5. representative genome chromosome
#   6. representative genome start
#   7. representative genome end
#   8. representative genome strand
#   9. represenative genome build
#  10. entrez gene id
#  11. NCBI gene chromosome
#  12. NCBI gene start
#  13. NCBI gene end
#  14. NCBI gene strand
#  15. ensembl gene id
#  16. ensembl gene chromosome
#  17. ensembl gene start
#  18. ensembl gene end
#  19. ensembl gene strand
#
# Usage:
#       MGI_Coordinate.py
#
# Used by:
# 	
#
# Notes:
#
# History:
#
# lec	02/17/2005
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

coordDisplay = '(%s:%s-%s (%s))'
noneDisplay = 'null' + TAB
repGenomicKey = 615419
sequenceType = 19
ncbi = 59
ensembl = 60

def getCoords(logicalDBkey):
    global repCoords

    tempCoords = {}

    results = db.sql('select m._Marker_key, mc._Qualifier_key, a.accID, ' + \
	    'c.chromosome, c.strand, ' + \
	    'startC = convert(int, c.startCoordinate), ' + \
	    'endC = convert(int, c.endCoordinate) ' + \
	        'from #markers m, SEQ_Marker_Cache mc, SEQ_Coord_Cache c, ACC_Accession a ' + \
	        'where m._Marker_key = mc._Marker_key ' + \
	        'and mc._Sequence_key = c._Sequence_key ' + \
	        'and mc._Sequence_key = a._Object_key ' + \
	        'and a._MGIType_key = %d ' % (sequenceType) + \
	        'and a._LogicalDB_key = %d ' % (logicalDBkey) , 'auto')

    for r in results:
        key = r['_Marker_key']
        value = r
        qualifier = r['_Qualifier_key']

	if qualifier == repGenomicKey:
	    repCoords[key] = value

        tempCoords[key] = value
 
    return tempCoords
    
#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

# markers that have coordinates
# (representative genomic)

db.sql('select m._Marker_key, m._Sequence_key, c.version into #repmarkers ' + \
	'from SEQ_Marker_Cache m, SEQ_Coord_Cache c ' + \
	'where m._Qualifier_key = %d ' % (repGenomicKey) + \
	'and m._Sequence_key = c._Sequence_key', None)

db.sql('create index idx1 on #repmarkers(_Marker_key)', None)
db.sql('create index idx2 on #repmarkers(_Sequence_key)', None)

# rest of marker information

db.sql('select m._Marker_key, a.accID, a.numericPart, mm.symbol, mm.name, markerType = t.name ' + 
	'into #markers ' + \
	'from #repmarkers m, ACC_Accession a, MRK_Marker mm, MRK_Types t ' + \
	'where m._Marker_key = a._Object_key ' + \
  	'and a._MGIType_key = 2 ' + \
  	'and a.prefixPart = "MGI:" ' + \
  	'and a._LogicalDB_key = 1 ' + \
  	'and a.preferred = 1 ' + \
	'and m._Marker_key = mm._Marker_key ' + \
	'and mm._Marker_Type_key = t._Marker_Type_key', None)

db.sql('create index idx1 on #markers(_Marker_key)', None)

# get genome build version

results = db.sql('select distinct version from #repmarkers', 'auto')
genomeBuild = results[0]['version']

# get coordinates

ncbiCoords = {}
ensemblCoords = {}
repCoords = {}
ncbiCoords = getCoords(ncbi)
ensemblCoords = getCoords(ensembl)

# process results

results = db.sql('select * from #markers order by numericPart', 'auto')

for r in results:
    key = r['_Marker_key']

    fp.write(r['accID'] + TAB)
    fp.write(r['markerType'] + TAB)
    fp.write(r['symbol'] + TAB)
    fp.write(r['name'] + TAB)

    rep = repCoords[key]
    fp.write(rep['chromosome'] + TAB)
    fp.write(str(rep['startC']) + TAB)
    fp.write(str(rep['endC']) + TAB)
    fp.write(rep['strand'] + TAB)
    fp.write(genomeBuild + TAB)

#    fp.write(coordDisplay % (rep['chromosome'], rep['startC'], rep['endC'], rep['strand']) + TAB)

    # NCBI coordinates

    if ncbiCoords.has_key(key):
	c = ncbiCoords[key]
	fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
#        fp.write(coordDisplay % (c['chromosome'], c['startC'], c['endC'], c['strand']) + TAB)
    else:
	fp.write(5*noneDisplay)

    # Ensembl coordinates

    if ensemblCoords.has_key(key):
	c = ensemblCoords[key]
	fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
#        fp.write(coordDisplay % (c['chromosome'], c['startC'], c['endC'], c['strand']) + TAB)
    else:
	fp.write(5*noneDisplay)

    fp.write(CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)