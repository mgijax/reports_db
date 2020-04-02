
'''
#
# MGI_MarkerNames.py
#
# Report:
#       Tab-delimited file
#	TR 10591
# Usage:
#       MGI_MarkerNames
#
#	output:
#	
#	field 1: MGI ID (primary)
#	field 2: Symbol
#	field 3: Marker Type
#	field 4: Feature Type (|-delimited): blank if no feature
#	field 5-?: Synonyms, Accession id
#
# Used by:
#       Matt Hibbs (Matt.Hibbs@jax.org)
#	Al Simons (Al.Simons@jax.org)
#
# History:
#
# lec	05/13/2015
#	- TR11652/convert to MRK_Cluster
#
# lec	01/19/2012
#	- include MRK_History + MGI_Synonym
#
# lec	01/18/2012
#	- use sequence cache to select both primary and secondary accession ids
#	- add non-mouse genes
#
# lec	01/11/2012
#	- select all synonyms
#
# lec	02/22/2011
#	- TR10591/new
#
'''

import sys
import os
import string
import mgi_utils
import reportlib
import db

db.setTrace()

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# all official mouse markers

db.sql('''select m._Marker_key, m.symbol, m.name, a.accID, a.numericPart, t.name as markerType
        into temporary table markers 
        from MRK_Marker m, ACC_Accession a, MRK_Types t
        where m._Organism_key = 1
        and m._Marker_Status_key = 1
        and m._Marker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a._LogicalDB_key = 1 
        and a.prefixPart = 'MGI:' 
        and a.preferred = 1
        and m._Marker_Type_key = t._Marker_Type_key
        ''', None)

db.sql('create index idx1 on markers(_Marker_key)', None)

#
# marker feature type
#
mgiFeature = {}
results = db.sql('''
           select m._Marker_key, t.term
           from markers m, VOC_Annot a, VOC_Term t
           where m._Marker_key = a._Object_key
           and a._AnnotType_key = 1011
           and a._Term_key = t._Term_key
        ''', 'auto')
for r in results:
    key = r['_Marker_key']
    value = r['term']
    if key not in mgiFeature:
        mgiFeature[key] = []
    mgiFeature[key].append(value)

#
# synonyms
#
results = db.sql('''select m._Marker_key, s.synonym 
        from markers m, MGI_Synonym s, MGI_SynonymType st 
        where m._Marker_key = s._Object_key 
        and s._MGIType_key = 2 
        and s._SynonymType_key = st._SynonymType_key 
        union
        select distinct m._Marker_key, h.history
        from markers m, MRK_History_View h
        where m._Marker_key = h._Marker_key
        ''', 'auto')
mgiSyn = {}
for r in results:
    key = r['_Marker_key']
    value = r['synonym']
    if key not in mgiSyn:
        mgiSyn[key] = []
    mgiSyn[key].append(value)

#
# non-mouse symbols
#
results = db.sql('''
        select distinct m._Marker_key, mm.symbol || '|' || o.commonName as synonym
        from markers m,
                MRK_Cluster mc, MRK_ClusterMember mcm, MRK_ClusterMember mcm2, 
                MRK_Marker mm, MGI_Organism o
        where mc._ClusterSource_key = 9272151
        and mc._Cluster_key = mcm._Cluster_key
        and mcm._Marker_key = m._Marker_key
        and mcm._Cluster_key = mcm2._Cluster_key
        and mcm2._Marker_key = mm._Marker_key
        and mm._Organism_key in (2,40)
        and mm._Organism_key = o._Organism_key
        ''', 'auto')
mgiNonMouse = {}
for r in results:
    key = r['_Marker_key']
    value = r['synonym']
    if key not in mgiNonMouse:
        mgiNonMouse[key] = []
    mgiNonMouse[key].append(value)

results = db.sql('''
        select distinct m._Marker_key, s.synonym || '|' || o.commonName as synonym
        from markers m, MRK_Cluster mc, MRK_ClusterMember mcm, MRK_ClusterMember mcm2, 
                MRK_Marker mm, MGI_Organism o, MGI_Synonym s
        where mc._ClusterSource_key = 9272151
        and mc._Cluster_key = mcm._Cluster_key
        and mcm._Marker_key = m._Marker_key
        and mcm._Cluster_key = mcm2._Cluster_key
        and mcm2._Marker_key = mm._Marker_key
        and mm._Organism_key in (2,40)
        and mm._Organism_key = o._Organism_key
        and mm._Marker_key = s._Object_key
        and s._MGIType_key = 2
        ''', 'auto')
for r in results:
    key = r['_Marker_key']
    value = r['synonym']
    if key not in mgiNonMouse:
        mgiNonMouse[key] = []
    mgiNonMouse[key].append(value)

# all data from ACC_Accession

results = db.sql('''select distinct m._Marker_key, a.accID
      from markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      union
      select distinct s._Marker_key, s.accID
      from markers m, SEQ_Marker_Cache s
      where m._Marker_key = s._Marker_key 
      ''', 'auto')
mgiID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if key not in mgiID:
        mgiID[key] = []
    mgiID[key].append(value)

# process

results = db.sql('select * from markers order by numericPart', 'auto')

for r in results:
        key = r['_Marker_key']

        fp.write(r['accID'] + reportlib.TAB)
        fp.write(r['symbol'] + reportlib.TAB)
        fp.write(r['markerType'] + reportlib.TAB)

        if key in mgiFeature:
            fp.write('|'.join(mgiFeature[key]))
        fp.write(reportlib.TAB)

        fp.write(r['name'] + reportlib.TAB)

        if key in mgiSyn:
            fp.write(reportlib.TAB.join(mgiSyn[key]) + reportlib.TAB)

        if key in mgiNonMouse:
            fp.write(reportlib.TAB.join(mgiNonMouse[key]) + reportlib.TAB)

        if key in mgiID:
            fp.write(reportlib.TAB.join(mgiID[key]))
        fp.write(reportlib.CRT)

reportlib.finish_nonps(fp)
