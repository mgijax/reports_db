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
#   5. representative genome id
#   6. representative genome chromosome
#   7. representative genome start
#   8. representative genome end
#   9. representative genome strand
#  10. represenative genome build
#  11. entrez gene id
#  12. NCBI gene chromosome
#  13. NCBI gene start
#  14. NCBI gene end
#  15. NCBI gene strand
#  16. Ensembl gene id
#  17. Ensembl gene chromosome
#  18. Ensembl gene start
#  19. Ensembl gene end
#  20. Ensembl gene strand
#  21. UniSTS gene start
#  22. UniSTS gene end
#
# Usage:
#       MGI_Coordinate.py
#
# History:
#
# 04/21/2006 lec
#	- added UniSTS coordinates
#
# 02/17/2005 lec
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
vega = 85
vegaprovider = 'VEGA Gene Model'
ncbi = 59
ncbiprovider = 'NCBI Gene Model'
ensembl = 60
ensemblprovider = 'Ensembl Gene Model'
unists = 80
unistsprovider = 'NCBI UniSTS'
qtl = 1
qtlprovider = 'MGI QTL'
mirbase = 83
mirbaseprovider = 'miRBase'

def getCoords(logicalDBkey, provider):

    global repCoords

    tempCoords = {}

    # we're assuming that markers may have a VEGA, NCBI and/or Ensembl coordinate
    # OR a UniSTS coordinate but not both

    if logicalDBkey in [vega, ncbi, ensembl]:
        results = db.sql('select m._Marker_key, a.accID, ' + \
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
            tempCoords[key] = value
 
    # UniSTS, QTL, mirBASE

    else:
        results = db.sql('select m._Marker_key, ' + \
	        'c.chromosome, c.strand, ' + \
	        'startC = convert(int, c.startCoordinate), ' + \
	        'endC = convert(int, c.endCoordinate) ' + \
	            'from #markers m, MRK_Location_Cache c ' + \
	            'where m._Marker_key = c._Marker_key ' + \
	            'and c.provider = "%s" ' % (provider), 'auto')

        for r in results:
            key = r['_Marker_key']
            value = r
            tempCoords[key] = value

    # get the representative coordinates

    results = db.sql('select m._Marker_key, a.accID, ' + \
	        'c.chromosome, c.strand, ' + \
	        'startC = convert(int, c.startCoordinate), ' + \
	        'endC = convert(int, c.endCoordinate), genomeBuild = c.version ' + \
	            'from #markers m, SEQ_Marker_Cache mc, SEQ_Coord_Cache c, ACC_Accession a ' + \
	            'where m._Marker_key = mc._Marker_key ' + \
		    'and mc._Qualifier_key = 615419 ' + \
		    'and mc._Sequence_key = c._Sequence_key ' + \
	            'and mc._Sequence_key = a._Object_key ' + \
	            'and a._MGIType_key = %d ' % (sequenceType) + \
	            'and a._LogicalDB_key = %d ' % (logicalDBkey) , 'auto')

    for r in results:
        key = r['_Marker_key']
        value = r
        repCoords[key] = value

    return tempCoords
    
#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('MGI accession id' + TAB)
fp.write('marker type' + TAB)
fp.write('marker symbol' + TAB)
fp.write('marker name' + TAB)
fp.write('representative genome id' + TAB)
fp.write('representative genome chromosome' + TAB)
fp.write('representative genome start' + TAB)
fp.write('representative genome end' + TAB)
fp.write('representative genome strand' + TAB)
fp.write('representative genome build' + TAB)
fp.write('Entrez gene id' + TAB)
fp.write('NCBI gene chromosome' + TAB)
fp.write('NCBI gene start' + TAB)
fp.write('NCBI gene end' + TAB)
fp.write('NCBI gene strand' + TAB)
fp.write('Ensembl gene id' + TAB)
fp.write('Ensembl gene chromosome' + TAB)
fp.write('Ensembl gene start' + TAB)
fp.write('Ensembl gene end' + TAB)
fp.write('Ensembl gene strand' + TAB)
fp.write('VEGA gene id' + TAB)
fp.write('VEGA gene chromosome' + TAB)
fp.write('VEGA gene start' + TAB)
fp.write('VEGA gene end' + TAB)
fp.write('VEGA gene strand' + TAB)
fp.write('UniSTS gene start' + TAB)
fp.write('UniSTS gene end' + CRT)

# all active markers

db.sql('select m._Marker_key, a.accID, a.numericPart, m.symbol, m.name, markerType = t.name ' + 
	'into #markers ' + \
	'from MRK_Marker m, MRK_Types t, ACC_Accession a ' + \
	'where m._Organism_key = 1 ' + \
	'and m._Marker_Status_key in (1,3) ' + \
	'and m._Marker_key = a._Object_key ' + \
  	'and a._MGIType_key = 2 ' + \
  	'and a.prefixPart = "MGI:" ' + \
  	'and a._LogicalDB_key = 1 ' + \
  	'and a.preferred = 1 ' + \
	'and m._Marker_Type_key = t._Marker_Type_key', None)

db.sql('create index idx1 on #markers(_Marker_key)', None)

# get coordinates

vegaCoords = {}
ncbiCoords = {}
ensemblCoords = {}
unistsCoords =  {}
qtlCoords = {}
mirbaseCoords = {}
repCoords = {}

vegaCoords = getCoords(vega, vegaprovider)
ncbiCoords = getCoords(ncbi, ncbiprovider)
ensemblCoords = getCoords(ensembl, ensemblprovider)
unistsCoords = getCoords(unists, unistsprovider)
qtlCoords = getCoords(qtl, qtlprovider)
mirbaseCoords = getCoords(mirbase, mirbaseprovider)

# process results

results = db.sql('select * from #markers order by numericPart', 'auto')

for r in results:
    key = r['_Marker_key']

    fp.write(r['accID'] + TAB)
    fp.write(r['markerType'] + TAB)
    fp.write(r['symbol'] + TAB)
    fp.write(r['name'] + TAB)

    # representative coordinate

    if repCoords.has_key(key):
        c = repCoords[key]
        fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(mgi_utils.prvalue(c['strand']) + TAB)
        fp.write(c['genomeBuild'] + TAB)
    else:
	fp.write(5*noneDisplay)

    # NCBI coordinate

    if ncbiCoords.has_key(key):
	c = ncbiCoords[key]
	fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
    else:
	fp.write(5*noneDisplay)

    # Ensembl coordinate

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

    # VEGA coordinate

    if vegaCoords.has_key(key):
	c = vegaCoords[key]
	fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
    else:
	fp.write(5*noneDisplay)

    # UniSTS coordinate

    if unistsCoords.has_key(key):
	c = unistsCoords[key]
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
    else:
	fp.write(2*noneDisplay)

    fp.write(CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
