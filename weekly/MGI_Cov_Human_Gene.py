"""

# Report: MGI_Cov_Human_Gene.py
# 
# http://wts.informatics.jax.org/searches/tr.detail.cgi?TR_Nr=TR13459
# 
# Filter the Alliance human disease file DISEASE-ALLIANCE_HUMAN.tsv  DISEASE-ALLIANCE_HUMAN_37.tsv
# and pull in gene symbol and MGI ID from MGD
#
# Output format: (1-7 from file, 8, 9 from database)
# 1. DBObjectID (HGNC ID, input file col 4)
# 2. DBObjectSymbol (Human Gene Symbol, input file col 5))
# 3. AssociationType (input file col 6)
# 4. DOID (input file col 7)
# 5. DOtermName  (input file col 8)
# 6. Reference  (input file col 14)
# 7. Source  (input file col 16
# 8. MGI gene ID
# 9. MGI gene symbol (ortholog of human gene in 1)
#
# skip comment lines(#) and column header
# skip if HGNC ID not in MGD
# skip if association type not in ['is_implicated_in', 'is_not_implicated_in']
# pull in alliance clustered mouse orthology if it exists
#
# Usage:
#       MGI_Cov_Human_Gene.py
#
# History:
#
#  sc  02/02/2021
#       - created
#
"""

import sys
import os
import string
import reportlib
import db

db.setTrace()

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Lookups
#
DOIDIncludeList = ['DOID:0080599','DOID:0080600','DOID:0080848', 'DOID:0080642', 'DOID:0080711', 'DOID:2945']
assocTypeIncludeList = ['is_implicated_in', 'is_not_implicated_in']

# {HGNC ID: ['mgiID|symbol, ...], ...}
hgncIdToMouseHomDict = {}

# input file
inFile = os.getenv('ALLIANCE_HUMAN_FILE')

fpIn = None
fpOut = None

def init():
    global hgncIdToMouseHomDict, fpIn, fpOut 

    fpOut = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

    # create Alliance Human Disease file descriptor
    fpIn = open(inFile, 'r') 

    # create HGNC ID to Mouse Homology Lookup
    db.sql('''-- get hgnc id and alliance cluster to which they belong
        select a.accid as hgncID, cm._cluster_key
        into temporary table clusterKeys
        from acc_accession a, mrk_cluster c, mrk_clustermember cm
        where c._clustersource_key = 75885740 -- alliance clustered
        and c._cluster_key = cm._cluster_key
        and cm._marker_key = a._object_key
        and a._mgitype_key = 2
        and a._LogicalDB_key = 64 -- hgnc''', None)

    db.sql('''create index idx1 on clusterKeys(_Cluster_key)''', None)

    #results = db.sql('''select a.accid, ck.hgncID, m.symbol, m._marker_key, m._organism_key, cm.*
    results = db.sql('''select a.accid, ck.hgncID, m.symbol
        from clusterKeys ck, MRK_ClusterMember cm, mrk_marker m, acc_accession a
        where ck._Cluster_key = cm._Cluster_key
        and cm._marker_key = m._marker_key
        and m._organism_key = 1
        and m._marker_key = a._object_key
        and a._mgitype_key = 2
        and a._logicalDB_key = 1
        and a.preferred = 1
        and a.prefixPart = 'MGI:'
        order by ck.hgncID, cm._cluster_key, cm.sequenceNum ''', 'auto')

    for r in results:
        hgncID = r['hgncID']
        mouseHomology = '%s|%s' % (r['accid'], r['symbol'])
        if hgncID not in hgncIdToMouseHomDict:
            hgncIdToMouseHomDict[hgncID] = []
        hgncIdToMouseHomDict[hgncID].append(mouseHomology)
#
# Main
#

init()
    
# Print the header record for the report.
#
fpOut.write('#\n')
fpOut.write('# This report contains human genes associated with coronavirus research references from the Alliance of Genome Resources.\n')
fpOut.write('# Each line contains an individual human gene, mouse ortholog and reference combination.\n')
fpOut.write('#\n')
fpOut.write('#1. DBObjectID (HGNC ID)\n')
fpOut.write('#2. DBObjectSymbol (Human Gene Symbol)\n')
fpOut.write('#3. AssociationType\n')
fpOut.write('#4. DOID\n')
fpOut.write('#5. DOtermName\n')
fpOut.write('#6. Reference\n')
fpOut.write('#7. Source\n')
fpOut.write('#8. MGI gene ID\n')
fpOut.write('#9. MGI gene symbol (ortholog of human gene in 1)\n')
fpOut.write('#\n')

for line in fpIn.readlines():
    if str.find(line, '#') == 0 or str.find(line, 'Taxon') == 0: # ignore comments, and header
        continue
    
    line = str.strip(line)
    tokens = str.split(line, TAB)
    DBObjectID = tokens[3]
    DBObjectSymbol = tokens[4]
    AssociationType = tokens[5]
    DOID = tokens[6]
    DOtermName = tokens[7]
    Reference = tokens[13]
    Source = tokens[15]
    symbol = '' # default if HGNC ID not in database
    mgiID = ''  # as we want to include these in the report
    if DOID not in DOIDIncludeList:
        continue
    if AssociationType not in assocTypeIncludeList:
        continue
    if DBObjectID not in hgncIdToMouseHomDict:
        fpOut.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (DBObjectID, TAB, DBObjectSymbol, TAB, AssociationType, TAB, DOID, TAB, DOtermName, TAB, Reference, TAB, Source, TAB, symbol, TAB, mgiID, CRT))
        continue
    mouseHomologyList = hgncIdToMouseHomDict[DBObjectID]
    for homology in mouseHomologyList:
         symbol, ID = str.split(homology, '|')
         fpOut.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (DBObjectID, TAB, DBObjectSymbol, TAB, AssociationType, TAB, DOID, TAB, DOtermName, TAB, Reference, TAB, Source, TAB, symbol, TAB, ID, CRT))

fpIn.close()
fpOut.close()

reportlib.finish_nonps(fpOut)
