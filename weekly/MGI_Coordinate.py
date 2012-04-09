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
#
#   5. representative genome id
#   6. representative genome chromosome
#   7. representative genome start
#   8. representative genome end
#   9. representative genome strand
#  10. represenative genome build
#
#  11. entrez gene id
#  12. NCBI gene chromosome
#  13. NCBI gene start
#  14. NCBI gene end
#  15. NCBI gene strand
#
#  16. Ensembl gene id
#  17. Ensembl gene chromosome
#  18. Ensembl gene start
#  19. Ensembl gene end
#  20. Ensembl gene strand
#
#  21. VEGA gene id
#  22. VEGA gene chromosome
#  23. VEGA gene start
#  24. VEGA gene end
#  25. VEGA gene strand
#
#  26. UniSTS gene chromosome
#  27. UniSTS gene start
#  28. UniSTS gene end
#
#  29. MGI QTL gene chromosome
#  30. MGI QTL gene start
#  31. MGI QTL gene end
#
#  32. miRBase gene id
#  33. miRBase gene chromosome
#  34. miRBase gene start
#  35. miRBase gene end
#  36. miRBase gene strand
#
#  37. Roopenian gene start
#  38. Roopenian gene end
#
# Usage:
#       MGI_Coordinate.py
#
# History:
#
# 02/15/2011	lec
#	- TR 10597/QTL/miRBase null counts are wrong
# 05/20/2010
#	- TR 10178/add miRBase gene id/chromosome/gene strand
#
# 07/30/2008 lec
#	- TR9152/Steve Grubb; add UniSTS chromosome
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
import db
import reportlib
import mgi_utils

CRT = reportlib.CRT
TAB = reportlib.TAB

coordDisplay = '(%s:%s-%s (%s))'
noneDisplay = 'null' + TAB
repGenomicKey = 615419
sequenceType = 19
mgiMarkerType = 2
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
roopenian = 1
roopenianprovider = 'Roopenian STS'

def getCoords(logicalDBkey, provider):

    global repCoords

    tempCoords = {}

    # we're assuming that markers may have a VEGA, NCBI and/or Ensembl coordinate
    # OR a UniSTS coordinate but not both
    #
    # VEGA, NCBI, Ensembl
    # contains accID, chromosome, strand, startC, endC
    #

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
 
    # miRBase
    # contains accid, chromosome, strand, startC, endC

    elif logicalDBkey in [mirbase]:

        results = db.sql('select m._Marker_key, a.accID, ' + \
                'c.chromosome, c.strand, ' + \
                'startC = convert(int, c.startCoordinate), ' + \
                'endC = convert(int, c.endCoordinate) ' + \
                    'from #markers m, MRK_Location_Cache c, ACC_Accession a ' + \
                    'where m._Marker_key = c._Marker_key ' + \
                    'and c.provider = "%s" ' % (provider) + \
		    'and m._Marker_key = a._Object_key ' + \
		    'and a._MGIType_key = %d ' % (mgiMarkerType) + \
		    'and a._LogicalDB_Key = %d ' % (logicalDBkey) , 'auto')

        for r in results:
            key = r['_Marker_key']
            value = r
            tempCoords[key] = value

    # UniSTS, QTL, Roopenian
    # contains chromosome, strand, startC, endC

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

fp.write('1. MGI accession id' + TAB)
fp.write('2. marker type' + TAB)
fp.write('3. marker symbol' + TAB)
fp.write('4. marker name' + TAB)

fp.write('5. representative genome id' + TAB)
fp.write('6. representative genome chromosome' + TAB)
fp.write('7. representative genome start' + TAB)
fp.write('8. representative genome end' + TAB)
fp.write('9. representative genome strand' + TAB)
fp.write('10. representative genome build' + TAB)

fp.write('11. Entrez gene id' + TAB)
fp.write('12. NCBI gene chromosome' + TAB)
fp.write('13. NCBI gene start' + TAB)
fp.write('14. NCBI gene end' + TAB)
fp.write('15. NCBI gene strand' + TAB)

fp.write('16. Ensembl gene id' + TAB)
fp.write('17. Ensembl gene chromosome' + TAB)
fp.write('18. Ensembl gene start' + TAB)
fp.write('19. Ensembl gene end' + TAB)
fp.write('20. Ensembl gene strand' + TAB)

fp.write('21. VEGA gene id' + TAB)
fp.write('22. VEGA gene chromosome' + TAB)
fp.write('23. VEGA gene start' + TAB)
fp.write('24. VEGA gene end' + TAB)
fp.write('25. VEGA gene strand' + TAB)

fp.write('26. UniSTS gene chromosome' + TAB)
fp.write('27. UniSTS gene start' + TAB)
fp.write('28. UniSTS gene end' + TAB)

fp.write('29. MGI QTL gene chromosome' + TAB)
fp.write('30. MGI QTL gene start' + TAB)
fp.write('31. MGI QTL gene end' + TAB)

fp.write('32. miRBase gene id' + TAB)
fp.write('33. miRBase gene chromosome' + TAB)
fp.write('34. miRBase gene start' + TAB)
fp.write('35. miRBase gene end' + TAB)
fp.write('36. miRBase gene strand' + TAB)

fp.write('37. Roopenian STS gene start' + TAB)
fp.write('38. Roopenian STS gene end' + CRT)

# all active markers

db.sql('''select m._Marker_key, a.accID, a.numericPart, m.symbol, m.name, markerType = t.name 
        into #markers 
        from MRK_Marker m, MRK_Types t, ACC_Accession a 
        where m._Organism_key = 1 
        and m._Marker_Status_key in (1,3) 
        and m._Marker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a.prefixPart = "MGI:" 
        and a._LogicalDB_key = 1 
        and a.preferred = 1 
        and m._Marker_Type_key = t._Marker_Type_key
	''', None)

db.sql('create index idx1 on #markers(_Marker_key)', None)

# get coordinates

vegaCoords = {}
ncbiCoords = {}
ensemblCoords = {}
unistsCoords =  {}
qtlCoords = {}
mirbaseCoords = {}
roopenianCoords = {}
repCoords = {}

vegaCoords = getCoords(vega, vegaprovider)
ncbiCoords = getCoords(ncbi, ncbiprovider)
ensemblCoords = getCoords(ensembl, ensemblprovider)
unistsCoords = getCoords(unists, unistsprovider)
qtlCoords = getCoords(qtl, qtlprovider)
mirbaseCoords = getCoords(mirbase, mirbaseprovider)
roopenianCoords = getCoords(roopenian, roopenianprovider)

# process results

results = db.sql('select * from #markers order by numericPart', 'auto')

for r in results:
    key = r['_Marker_key']

    # 1-4
    fp.write(r['accID'] + TAB)
    fp.write(r['markerType'] + TAB)
    fp.write(r['symbol'] + TAB)
    fp.write(r['name'] + TAB)

    # representative coordinate (5-10)

    if repCoords.has_key(key):
        c = repCoords[key]
        fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(mgi_utils.prvalue(c['strand']) + TAB)
        fp.write(c['genomeBuild'] + TAB)
    else:
        fp.write(6*noneDisplay)

    # NCBI coordinate (11-15)

    if ncbiCoords.has_key(key):
        c = ncbiCoords[key]
        fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
    else:
        fp.write(5*noneDisplay)

    # Ensembl coordinate (16-20)

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

    # VEGA coordinate (21-25)

    if vegaCoords.has_key(key):
        c = vegaCoords[key]
        fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
    else:
        fp.write(5*noneDisplay)

    # UniSTS coordinate (26-28)

    if unistsCoords.has_key(key):
        c = unistsCoords[key]
        fp.write(str(c['chromosome']) + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
    else:
        fp.write(3*noneDisplay)

    # MGI QTL coordinate (29-31)

    if qtlCoords.has_key(key):
        c = qtlCoords[key]
        fp.write(str(c['chromosome']) + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
    else:
        fp.write(3*noneDisplay)

    # miRBase coordinate (32-36)

    if mirbaseCoords.has_key(key):
        c = mirbaseCoords[key]
        fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
    else:
        fp.write(5*noneDisplay)

    # Roopenian coordinate (37-38)

    if roopenianCoords.has_key(key):
        c = roopenianCoords[key]
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
    else:
        fp.write(2*noneDisplay)

    fp.write(CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
