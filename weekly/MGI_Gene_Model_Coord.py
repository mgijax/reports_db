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
#   5. genome build
#
#   6. entrez gene id
#   7. NCBI gene chromosome
#   8. NCBI gene start
#   9. NCBI gene end
#  10. NCBI gene strand
#
#  11. Ensembl gene id
#  12. Ensembl gene chromosome
#  13. Ensembl gene start
#  14. Ensembl gene end
#  15. Ensembl gene strand
#
#  16. VEGA gene id
#  17. VEGA gene chromosome
#  18. VEGA gene start
#  19. VEGA gene end
#  20. VEGA gene strand
#
# Usage:
#       MGI_gene_model_coord.py
#
# History:
#
# 07/26/2012 - jsb - built from now-defunct MGI_Coordinate.py by omitting
#	columns to the right of the VEGA strand (part of the Coordinates for
#	any Marker - C4AM - project)
#
'''
 
import sys
import os
import mgi_utils
import reportlib
import db

db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE()

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

genomeBuild = None

def getCoords(logicalDBkey, provider):

    global genomeBuild

    tempCoords = {}

    # we're assuming that markers may have a VEGA, NCBI and/or Ensembl coordinate
    # OR a UniSTS coordinate but not both
    #
    # VEGA, NCBI, Ensembl
    # contains accID, chromosome, strand, startC, endC
    #

    if logicalDBkey in [vega, ncbi, ensembl]:
        results = db.sql('''
		select m._Marker_key, a.accID, 
                       c.chromosome, c.strand, 
                       convert(int, c.startCoordinate) as startC, 
                       convert(int, c.endCoordinate) as endC 
                    from #markers m, SEQ_Marker_Cache mc, SEQ_Coord_Cache c, ACC_Accession a 
                    where m._Marker_key = mc._Marker_key 
                    and mc._Sequence_key = c._Sequence_key 
                    and mc._Sequence_key = a._Object_key 
                    and a._MGIType_key = %d
                    and a._LogicalDB_key = %d 
		 ''' % (sequenceType, logicalDBkey) , 'auto')

        for r in results:
            key = r['_Marker_key']
            value = r
            tempCoords[key] = value
    else:
	raise 'BadValueError', 'Unknown logicalDBkey in getCoords: %s' % \
		logicalDBkey
 
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

fp.write('5. genome build' + TAB)

fp.write('6. Entrez gene id' + TAB)
fp.write('7. NCBI gene chromosome' + TAB)
fp.write('8. NCBI gene start' + TAB)
fp.write('9. NCBI gene end' + TAB)
fp.write('10. NCBI gene strand' + TAB)

fp.write('11. Ensembl gene id' + TAB)
fp.write('12. Ensembl gene chromosome' + TAB)
fp.write('13. Ensembl gene start' + TAB)
fp.write('14. Ensembl gene end' + TAB)
fp.write('15. Ensembl gene strand' + TAB)

fp.write('16. VEGA gene id' + TAB)
fp.write('17. VEGA gene chromosome' + TAB)
fp.write('18. VEGA gene start' + TAB)
fp.write('19. VEGA gene end' + TAB)
fp.write('20. VEGA gene strand' + CRT)

# all active markers

db.sql('''select m._Marker_key, a.accID, a.numericPart, m.symbol, m.name, t.name as markerType
        into #markers 
        from MRK_Marker m, MRK_Types t, ACC_Accession a 
        where m._Organism_key = 1 
        and m._Marker_Status_key in (1,3) 
        and m._Marker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a.prefixPart = 'MGI:' 
        and a._LogicalDB_key = 1 
        and a.preferred = 1 
        and m._Marker_Type_key = t._Marker_Type_key
	''', None)

db.sql('create index idx1 on #markers(_Marker_key)', None)

# get coordinates

vegaCoords = {}
ncbiCoords = {}
ensemblCoords = {}

vegaCoords = getCoords(vega, vegaprovider)
ncbiCoords = getCoords(ncbi, ncbiprovider)
ensemblCoords = getCoords(ensembl, ensemblprovider)

# get genome build for mouse markers (avoid human markers, which would have
# a different build number)

results = db.sql ('''select distinct version
	from MRK_Location_Cache
	where version is not null
	and _Organism_key = 1''', 'auto')
if results:
	genomeBuild = results[0]['version']
else:
	genomeBuild = 'Unknown'

# process results

results = db.sql('select * from #markers order by numericPart', 'auto')

for r in results:
    key = r['_Marker_key']

    # if no gene models, skip this marker
    if not (ncbiCoords.has_key(key) or ensemblCoords.has_key(key) or \
	vegaCoords.has_key(key)):
	    continue
    
    # 1-4
    fp.write(r['accID'] + TAB)
    fp.write(r['markerType'] + TAB)
    fp.write(r['symbol'] + TAB)
    fp.write(r['name'] + TAB)

    # genome build (5)
    fp.write(genomeBuild + TAB)

    # NCBI coordinate (6-10)

    if ncbiCoords.has_key(key):
        c = ncbiCoords[key]
        fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
    else:
        fp.write(5*noneDisplay)

    # Ensembl coordinate (11-15)

    if ensemblCoords.has_key(key):
        c = ensemblCoords[key]
        fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
    else:
        fp.write(5*noneDisplay)

    # VEGA coordinate (16-20)

    if vegaCoords.has_key(key):
        c = vegaCoords[key]
        fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
    else:
        fp.write(5*noneDisplay)

    fp.write(CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
