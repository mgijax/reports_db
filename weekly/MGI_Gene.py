#!/usr/local/bin/python

'''
#
# Report:
#       Tab-delimited file of MGI Genes
#	(TR 11968)
#
# Usage:
#       MGI_gene_model_coord.py
#
# History:
#
# 04/16/2015 - sc - created
#
'''
 
import sys
import os
import mgi_utils
import reportlib
import string
import pg_db
db = pg_db
db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE()

CRT = reportlib.CRT
TAB = reportlib.TAB

noneDisplay = 'null' + TAB
sequenceType = 19
mgiMarkerType = 2
vega = 85
vegaprovider = 'VEGA Gene Model'
ncbi = 59
ncbiprovider = 'NCBI Gene Model'
ensembl = 60
ensemblprovider = 'Ensembl Gene Model'

ccdsDict = {}
ccds = 136

# {clusterKey: [ [], [] ], ...}
hgncClusterDict = {}

# {humanMarkerKey:hgncID, ...}
hgncDict = {}

# {mouseMarkerKey:[human HGNC IDs], ...}
mouseToHGNCDict = {}

hgnc = 64

homologeneDict = {}
homologene = 81

def getCoords(logicalDBkey, provider):


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

def loadLookups():
    global hgncDict, mouseToHGNCDict, hgncClusterDict, ccdsDict, homologeneDict

    # load lookup mapping human marker key to its HGNC ID    
    results = db.sql('''select accID as hgncID, _Object_key as hMarkerKey
	from ACC_Accession
	where _LogicalDB_key = 64
	and _MGIType_key = 2''', 'auto')

    for r in results:
	hgncDict[r['hMarkerKey']] = r['hgncID']

    # load lookup mapping hgnc cluster IDs to their members (mouse and human)
    # type: homology
    # source: HGNC
    results = db.sql('''select cm._Cluster_key, cm._Marker_key, m._Organism_key
	from MRK_Cluster c, MRK_ClusterMember cm, MRK_Marker m
	where c._ClusterType_key = 9272150
	and c._ClusterSource_key = 13437099
	and c._Cluster_key = cm._Cluster_key
	and cm._Marker_key = m._Marker_key''', 'auto')
    for r in results:
	clusterKey = r['_Cluster_key']
	markerKey = r['_Marker_key']
	orgKey = r['_Organism_key']
	if not hgncClusterDict.has_key(clusterKey):
	   hgncClusterDict[clusterKey] = [[],[]]
	if orgKey == 1:
	    hgncClusterDict[clusterKey][0].append(markerKey)
   	else:
	    hgncClusterDict[clusterKey][1].append(markerKey)

    # load lookup mapping mouse marker key to HGNC IDs of its human homologs
    for clusterKey in hgncClusterDict.keys():
	mouseList = hgncClusterDict[clusterKey][0]
	humanList = hgncClusterDict[clusterKey][1]
	hgncList = []
	for hKey in humanList:
	    if hgncDict.has_key(hKey):
		hgncList.append(hgncDict[hKey])
	if hgncList == []:
	    # no HGNC IDS
	    continue
	for mKey in mouseList:
	    mouseToHGNCDict[mKey] = hgncList
	    
    results = db.sql('''select _Object_key as markerKey, accID
        from ACC_Accession
        where _MGIType_key = 2
        and preferred = 1
        and _LogicalDB_key = %s''' % ccds, 'auto')
    for r in results:
	mKey = r['markerKey']
	id =  r['accID']
	if not ccdsDict.has_key(mKey):
	    ccdsDict[mKey] = []
	ccdsDict[mKey].append(id)

    results = db.sql('''select _Object_key as markerKey, accID
        from ACC_Accession
        where _MGIType_key = 2
        and preferred = 1
        and _LogicalDB_key = %s''' % homologene, 'auto')
    for r in results:
        homologeneDict[r['markerKey']] = r['accID']

	
     
#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('MGI Accession ID' + TAB)
fp.write('Marker Symbol' + TAB)
fp.write('Marker Name' + TAB)
fp.write('Feature Type' + TAB)

fp.write('EntrezGene ID' + TAB)
fp.write('NCBI Gene chromosome' + TAB)
fp.write('NCBI Gene start' + TAB)
fp.write('NCBI Gene end' + TAB)
fp.write('NCBI Gene strand' + TAB)

fp.write('Ensembl Gene ID' + TAB)
fp.write('Ensembl Gene chromosome' + TAB)
fp.write('Ensembl Gene start' + TAB)
fp.write('Ensembl Gene end' + TAB)
fp.write('Ensembl Gene strand' + TAB)

fp.write('VEGA Gene ID' + TAB)
fp.write('VEGA Gene chromosome' + TAB)
fp.write('VEGA Gene start' + TAB)
fp.write('VEGA Gene end' + TAB)
fp.write('VEGA Gene strand' + TAB)
fp.write('CCDS IDs' + TAB)
fp.write('HGNC ID' + TAB)
fp.write('HomoloGene ID' + CRT)
# all active markers

db.sql('''select m._Marker_key, a.accID, a.numericPart, m.symbol, m.name, t.term as featureType
        into #markers
        from MRK_Marker m, ACC_Accession a, VOC_Annot v, VOC_Term t
        where m._Organism_key = 1
        and m._Marker_Status_key in (1,3)
        and m._Marker_Type_key = 1
        and m._Marker_key = a._Object_key
        and a._MGIType_key = 2
        and a.prefixPart = 'MGI:'
        and a._LogicalDB_key = 1
        and a.preferred = 1
        and m._Marker_key = v._Object_key
        and v._AnnotType_key = 1011
	and v._Term_key = t._Term_key
	''', None)

db.sql('create index idx1 on #markers(_Marker_key)', None)

# get coordinates

vegaCoords = {}
ncbiCoords = {}
ensemblCoords = {}

vegaCoords = getCoords(vega, vegaprovider)
ncbiCoords = getCoords(ncbi, ncbiprovider)
ensemblCoords = getCoords(ensembl, ensemblprovider)

loadLookups()

# process results

results = db.sql('select * from #markers order by numericPart', 'auto')

for r in results:
    key = r['_Marker_key']

    # 1-3
    fp.write(r['accID'] + TAB)
    fp.write(r['symbol'] + TAB)
    fp.write(r['name'] + TAB)
    fp.write(r['featureType'] + TAB)
    # NCBI coordinate (4-8)

    if ncbiCoords.has_key(key):
        c = ncbiCoords[key]
        fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
    else:
        fp.write(5*noneDisplay)

    # Ensembl coordinate (9-13)

    if ensemblCoords.has_key(key):
        c = ensemblCoords[key]
        fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
    else:
        fp.write(5*noneDisplay)

    # VEGA coordinate (14-18)

    if vegaCoords.has_key(key):
        c = vegaCoords[key]
        fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
    else:
        fp.write(5*noneDisplay)

    # CCDS (19)
    if ccdsDict.has_key(key):
	ids = string.join(ccdsDict[key], ',')
	fp.write(ids + TAB)
    else:
	fp.write(noneDisplay)

    # HGNC (20)
    if mouseToHGNCDict.has_key(key):
        ids = string.join(mouseToHGNCDict[key], ',')
        fp.write(ids + TAB)
    else:
        fp.write(noneDisplay)

    # HomoloGene (21)
    if homologeneDict.has_key(key):
        id = homologeneDict[key]
        fp.write(id + TAB)
    else:
        fp.write(noneDisplay)

    fp.write(CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
