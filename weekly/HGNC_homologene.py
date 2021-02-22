
'''
#
# Report:
#       Tab-delimited file of MGI Genes
#	(TR 11968)
#
# Usage:
#       HGNC_homologene.py
#
# History:
#
# sc    02/22/2021
#       TR13349 - B39 project. Update to use alliance direct homology
#               (was using Homologene and HGNC)
#
# 11/24/2015 - TR11968
#	updated the homology
#
# 04/16/2015 - sc - created
#
'''
 
import sys
import os
import mgi_utils
import reportlib
import string
import db
import Set

db.setTrace()

CRT = reportlib.CRT
TAB = reportlib.TAB

noneDisplay = 'null' + TAB
sequenceType = 19
mgiMarkerType = 2
ncbi = 59
ncbiprovider = 'NCBI Gene Model'
ensembl = 60
ensemblprovider = 'Ensembl Gene Model'

ccdsDict = {}
ccds = 136

# mouse marker feature types
# {markerKey: featureTypeKey, ...}
featureTypeDict = {}

# {mouseKey: [humanKeys], ...} for from MRK_Cluster*
homologyDict = {}

# HGNC {humanKey: hgncId, ...}
hgncIdDict = {}
# HG {humanKey: [hgId1, ...], ...} 
hgeneIdDict = {}

hgnc = 64
homologene = 81

def getCoords(logicalDBkey, provider):


    tempCoords = {}

    # we're assuming that markers may have a NCBI and/or Ensembl coordinate
    # OR a UniSTS coordinate but not both
    #
    # NCBI, Ensembl
    # contains accID, chromosome, strand, startC, endC
    #

    if logicalDBkey in [ncbi, ensembl]:
        results = db.sql('''
                select m._Marker_key, a.accID, 
                       c.chromosome, c.strand, 
                       c.startCoordinate::int as startC, 
                       c.endCoordinate::int as endC 
                    from markers m, SEQ_Marker_Cache mc, SEQ_Coord_Cache c, ACC_Accession a 
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
        raise ValueError('Unknown logicalDBkey in getCoords: %s' % \
                logicalDBkey)
 
    return tempCoords

def loadLookups():
    global homologyDict, hgncIdDict, hgeneIdDict, ccdsDict, featureTypeDict

    # load lookup mapping mouse marker key to human marker keys from both HGNC and HG clusters
    # while we are at it, load mapping of human markers to their HGene ids
    results = db.sql('''select cm._Cluster_key, cm._Marker_key, m._Organism_key, c.clusterId
        from MRK_Cluster c, MRK_ClusterMember cm, MRK_Marker m
        where c._ClusterType_key = 9272150
        and c._ClusterSource_key = 75885739
        and c._Cluster_key = cm._Cluster_key
        and cm._Marker_key = m._Marker_key''', 'auto')

    tempClusterDict = {}
    for r in results:
        clusterKey = r['_Cluster_key']
        markerKey = r['_Marker_key']
        orgKey = r['_Organism_key']
        hgeneId = r['clusterId']
        if clusterKey not in tempClusterDict:
           tempClusterDict[clusterKey] = [[],[]]
        if orgKey == 1:
            tempClusterDict[clusterKey][0].append(markerKey)
             # if this is a homologene cluster (aka has a clusterId (HGNC is null))
            if hgeneId != None:
                if markerKey not in hgeneIdDict:
                    hgeneIdDict[markerKey]= []
                hgeneIdDict[markerKey].append(hgeneId)
        else:
            tempClusterDict[clusterKey][1].append(markerKey)

    # load lookup mapping mouse marker key to the set of all human markers it has homology with
    for clusterKey in list(tempClusterDict.keys()):
        mouseList = tempClusterDict[clusterKey][0]
        humanList = tempClusterDict[clusterKey][1]
        hgncList = []
        for mKey in mouseList:
            if mKey not in homologyDict:
                homologyDict[mKey] = []
            homologyDict[mKey] =  homologyDict[mKey] + humanList

    # load lookup mapping human markers to their HGNC IDs
    results = db.sql('''select accid, _Object_key
        from ACC_Accession
        where _MGIType_key = 2 
        and _LogicalDB_key = 64''', 'auto')
    for r in results:
        hgncIdDict[r['_Object_key']] = r['accid']

    # load ccds lookup
    results = db.sql('''select _Object_key as markerKey, accID
        from ACC_Accession
        where _MGIType_key = 2
        and preferred = 1
        and _LogicalDB_key = %s''' % ccds, 'auto')
    for r in results:
        mKey = r['markerKey']
        id =  r['accID']
        if mKey not in ccdsDict:
            ccdsDict[mKey] = []
        ccdsDict[mKey].append(id)

    # load feature type lookup
    results = db.sql('''select m._Marker_key, t.term 
        from MRK_Marker m, VOC_Annot v, VOC_Term t
        where m._Organism_key = 1
        and m._Marker_Status_key = 1
        and m._Marker_Type_key in (1,7)
        and m._Marker_Type_key = 1
        and m._Marker_key = v._Object_key
        and v._AnnotType_key = 1011
        and v._Term_key = t._Term_key''', 'auto')

    for r in results:
        markerKey = r['_Marker_key']
        featureType = r['term']
        featureTypeDict[markerKey] = featureType
        
     
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

fp.write('CCDS IDs' + TAB)
fp.write('HGNC ID' + TAB)
fp.write('HomoloGene ID' + CRT)

# all active genes/pseudogenes
db.sql('''select m._Marker_key, a.accID, a.numericPart, m.symbol, m.name
        into temporary table markers
        from MRK_Marker m, ACC_Accession a
        where m._Organism_key = 1
        and m._Marker_Status_key = 1
        and m._Marker_Type_key in (1,7)
        and m._Marker_key = a._Object_key
        and a._MGIType_key = 2
        and a.prefixPart = 'MGI:'
        and a._LogicalDB_key = 1
        and a.preferred = 1
        ''', None)

db.sql('create index idx1 on markers(_Marker_key)', None)

# get coordinates

ncbiCoords = {}
ensemblCoords = {}

ncbiCoords = getCoords(ncbi, ncbiprovider)
ensemblCoords = getCoords(ensembl, ensemblprovider)

loadLookups()

# process results

results = db.sql('select * from markers order by numericPart', 'auto')

for r in results:
    key = r['_Marker_key']

    # 1-4
    fp.write(r['accID'] + TAB)
    fp.write(r['symbol'] + TAB)
    fp.write(r['name'] + TAB)
    if key in featureTypeDict:
        fp.write(featureTypeDict[key] + TAB)
    else:
        fp.write(noneDisplay)

    # NCBI coordinate (5-9)

    if key in ncbiCoords:
        c = ncbiCoords[key]
        fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
    else:
        fp.write(5*noneDisplay)

    # Ensembl coordinate (10-14)

    if key in ensemblCoords:
        c = ensemblCoords[key]
        fp.write(c['accID'] + TAB)
        fp.write(c['chromosome'] + TAB)
        fp.write(str(c['startC']) + TAB)
        fp.write(str(c['endC']) + TAB)
        fp.write(c['strand'] + TAB)
    else:
        fp.write(5*noneDisplay)

    # CCDS (20)
    if key in ccdsDict:
        ids = '|'.join(ccdsDict[key])
        fp.write(ids + TAB)
    else:
        fp.write(noneDisplay)

    # HGNC and HomoloGene(21-22)
    # If the mouse key has no HG or HGNC homology leave blank
    if key not in homologyDict:
        fp.write(noneDisplay)
        fp.write(noneDisplay)
    else:
        # list of human markers this mouse has homology with
        humanList = homologyDict[key]
        hgncList = []
        hgeneList = []

        # get the HGene ID associated with  the human homologs for this mouse marker
        if key in hgeneIdDict:
            hgeneList = hgeneIdDict[key]

        # get the HGNC IDs associated with the human homologs for this mouse marker
        for hKey in humanList:
            if hKey in hgncIdDict:
                hgncList.append(hgncIdDict[hKey])

        # create sets, there could be dups
        hgncSet = set(hgncList)
        hgeneSet = set(hgeneList)
        # write out each set pipe delimited
        if hgncSet:
            fp.write('|'.join(hgncSet) + TAB)
        else:
            fp.write(noneDisplay)
        if hgeneSet: 
            fp.write('|'.join(hgeneSet) + TAB)
        else:
            fp.write(noneDisplay)
            
    fp.write(CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
