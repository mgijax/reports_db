#!/usr/local/bin/python

'''
#
# HMD_HumanPhenotype.py 
#
# Report:
#     A tab-delimited report of all mouse-human orthologous gene pairs, 
#     with all MP IDs associated with any allele of the mouse gene.
#
#     Output fields: 
#     
#     Human Gene Symbol 
#     Human EntrezGene ID 
#     Mouse Gene Symbol 
#     MGI Accession ID 
#     MP IDs associated with any allele of the mouse gene (space- or comma-separated) 
#     Homologene ID
#     HGNC yes/no
#
#     Sorted by:
#     Human Gene Symbol. 
#    
# Usage:
#       HMD_HumanPhenotype.py
#
# History:
#
# lec	03/26/2019
#	- TR13037/use rollup rules (_annottype_key = 1015) instead of (_annottype_key = 1002)
#
# lec	07/06/2017
#	- TR12615/exclude obsolete MP terms
#
# sc	11/14/2016
#	- TR12404 Use Hybrid Homology, report  Homologene cluster ID and HGNC y/n
#
# sc	04/19/2013
#	- N2MO; update to use MRK_Cluster* tables and deal w/
#	  N-M human/mouse homologies; added homologeneID to report
#
# lec	04/23/2008
#	- TR 8963; only select primary MP Header Ids
#
# lec	01/04/2004
#	- TR 5939; EntrezGene->EntrezGene
#
# lec	08/27/2003
#	- TR 5082; make a public report
#
# dm    07/14/03
#       -Inital Creation
# 
'''
 
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

def createDict(results, keyField, valueField):

        d = {}
        for r in results:
                key = r[keyField]
                value = r[valueField]
                if not d.has_key(key):
                    d[key] = ' '
		if d[key].find(value) == -1:
		    d[key] = string.join([d[key],value,])
        return d


#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# create list of HGNC mouse keys
hgncMouseList = []
results = db.sql('''select mcm._Marker_Key
    from MRK_Cluster mc, MRK_ClusterMember mcm, MRK_Marker m
    where mc._ClusterSource_key = 13437099
    and mc._Cluster_key = mcm._Cluster_key
    and mcm._Marker_key = m._Marker_key
    and m._Organism_key = 1''', 'auto')
for r in results:
    if r['_Marker_key'] in hgncMouseList:
	print 'mouse in multi HGNC clusters: %s' % r['_Marker_key']
    hgncMouseList.append(r['_Marker_key'])

# create dictionary of HG mouse keys -> HG ID
hgMouseDict = {}
results = db.sql('''select mcm._Marker_Key, mc.clusterID
    from MRK_Cluster mc, MRK_ClusterMember mcm, MRK_Marker m
    where mc._ClusterSource_key = 9272151
    and mc._Cluster_key = mcm._Cluster_key
    and mcm._Marker_key = m._Marker_key
    and m._Organism_key = 1''', 'auto')
for r in results:
    key = r['_Marker_key']
    id = r['clusterId']
    if key in hgMouseDict:
	print 'mouse in multi HG clusters: %s %s %s' % (key, id, hgMouseDict[key])
    hgMouseDict[key] = id

#
# Get mouse to human orthologous marker pair's symbols and keys
#
# Hybrid = 13764519 HomoloGene = 9272151
db.sql('''select c.clusterID, cm.*
    into temporary table mouse
    from MRK_Cluster c, MRK_ClusterMember cm, MRK_Marker m
    where c._ClusterType_key = 9272150
    and c._ClusterSource_key = 13764519   
    and c._Cluster_key = cm._Cluster_key
    and cm._Marker_key = m._Marker_key
    and m._Organism_key = 1''', None)

db.sql('create index mouse_idx1 on mouse(_Cluster_key)', None)

db.sql('''select cm.*
    into temporary table human
    from MRK_Cluster c, MRK_ClusterMember cm, MRK_Marker m
    where c._ClusterType_key = 9272150
    and c._ClusterSource_key = 13764519   
    and c._Cluster_key = cm._Cluster_key
    and cm._Marker_key = m._Marker_key
    and m._Organism_key = 2''', None)

db.sql('create index human_idx1 on human(_Cluster_key)', None)

db.sql('''
	select distinct hm._Marker_key as mouseKey, 
			m1.symbol as mouseSym, 
	                hh._Marker_key as humanKey, 
			m2.symbol as humanSym 
	into temporary table homology 
	from mouse hm, human hh, 
	MRK_Marker m1, MRK_Marker m2 
	where hm._Cluster_key = hh._Cluster_key
	and hm._Marker_key = m1._Marker_key 
	and hh._Marker_key = m2._Marker_key 
	''', None)

db.sql('create index index_mouseKey on homology(mouseKey)', None)
db.sql('create index index_humanKey on homology(humanKey)', None)
db.sql('create index index_humanSym on homology(humanSym)', None)

#
# Get the MGI IDs for the mouse markers
#

results = db.sql('''
	select a._Object_key, a.accID 
	from homology h, ACC_Accession a 
	where a._Object_key = h.mouseKey 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 1 
	and a.prefixPart = 'MGI:' 
	and a.preferred = 1 
	''', 'auto')
mmgi = createDict(results, '_Object_key', 'accID')

#
# Get the EntrezGene for the human markers
#

results = db.sql('''
	select a._Object_key, a.accID 
        from homology h, ACC_Accession a 
        where a._Object_key = h.humanKey 
        and a._MGIType_key = 2 
        and a._LogicalDB_key = 55
	''', 'auto')
hlocus = createDict(results, '_Object_key', 'accID')

#
# Get the MP Header Terms for the mouse markers using 1015/MP->Marker (derivied) annotations
# find the SetMember/ancestor of MP/descendent
#
# only select primary ids and non-obsolete ids
#
results = db.sql('''
	select distinct h.mouseKey, a.accID 
        from homology h, VOC_Annot v, MGI_SetMember sm, DAG_Closure dc, ACC_Accession a
	where h.mouseKey = v._Object_key
        and v._AnnotType_key = 1015
        and v._Term_key = dc._DescendentObject_key
        and dc._AncestorObject_key = sm._Object_key
        and sm._Set_key = 1051
        and sm._Object_key = a._Object_key
        and a._MGIType_key = 13
        and a.preferred = 1
        and a._LogicalDB_key = 34
	''', 'auto')
mpheno = createDict(results, 'mouseKey', 'accID')

results = db.sql('select * from homology order by humanSym, mouseSym', 'auto')

for r in results:

    hasHGNC = 'no'
    clusterID = None

    mouseKey = r['mouseKey']

    if mouseKey in  hgncMouseList:
        hasHGNC = 'yes'

    if mouseKey in hgMouseDict:
        clusterID = hgMouseDict[mouseKey]

    fp.write(r['humanSym'] + TAB)

    if hlocus.has_key(r['humanKey']):
        fp.write(hlocus[r['humanKey']] + TAB)
    else:
        fp.write(TAB)

    if clusterID == None:
	fp.write(TAB)
    else:
	fp.write(clusterID + TAB)

    fp.write(hasHGNC +  TAB)
    fp.write(r['mouseSym'] + TAB)

    if mmgi.has_key(mouseKey):
        fp.write(mmgi[mouseKey] + TAB)

    if mpheno.has_key(mouseKey):
        fp.write(mpheno[mouseKey] + TAB)
    else:
        fp.write(TAB)

    fp.write(CRT)

reportlib.finish_nonps(fp)
