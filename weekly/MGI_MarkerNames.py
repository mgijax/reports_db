#!/usr/local/bin/python

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
# lec	02/22/2011
#	- TR10591/new
#
'''

import sys
import os
import string
import db
import mgi_utils
import reportlib

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# all official/interim mouse markers

db.sql('''select m._Marker_key, m.symbol, m.name, a.accID, a.numericPart, markerType = t.name
        into #markers 
        from MRK_Marker m, ACC_Accession a, MRK_Types t
        where m._Organism_key = 1
        and m._Marker_Status_key in (1,3) 
        and m._Marker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a._LogicalDB_key = 1 
        and a.prefixPart = "MGI:" 
	and a.preferred = 1
	and m._Marker_Type_key = t._Marker_Type_key
	''', None)

db.sql('create index idx1 on #markers(_Marker_key)', None)

#
# marker feature type
#
mgiFeature = {}
results = db.sql('''
           select m._Marker_key, t.term
           from #markers m, VOC_Annot a, VOC_Term t
           where m._Marker_key = a._Object_key
           and a._AnnotType_key = 1011
           and a._Term_key = t._Term_key
        ''', 'auto')
for r in results:
    key = r['_Marker_key']
    value = r['term']
    if not mgiFeature.has_key(key):
        mgiFeature[key] = []
    mgiFeature[key].append(value)

#
# synonyms
#
results = db.sql('''select m._Marker_key, s.synonym 
	from #markers m, MGI_Synonym s, MGI_SynonymType st 
	where m._Marker_key = s._Object_key 
	and s._MGIType_key = 2 
	and s._SynonymType_key = st._SynonymType_key 
	and st.synonymType = "exact"
	''', 'auto')
mgiSyn = {}
for r in results:
    key = r['_Marker_key']
    value = r['synonym']
    if not mgiSyn.has_key(key):
        mgiSyn[key] = []
    mgiSyn[key].append(value)

# all data from ACC_Accession

results = db.sql('''select distinct m._Marker_key, a.accID
      from #markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      ''', 'auto')
mgiID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if not mgiID.has_key(key):
        mgiID[key] = []
    mgiID[key].append(value)

# process

results = db.sql('select * from #markers order by numericPart', 'auto')

for r in results:
	key = r['_Marker_key']

	fp.write(r['accID'] + reportlib.TAB)
        fp.write(r['symbol'] + reportlib.TAB)
        fp.write(r['markerType'] + reportlib.TAB)

	if mgiFeature.has_key(key):
            fp.write(string.join(mgiFeature[key], '|'))
        fp.write(reportlib.TAB)

	fp.write(r['name'] + reportlib.TAB)

	if mgiSyn.has_key(key):
	    fp.write(string.join(mgiSyn[key], reportlib.TAB) + reportlib.TAB)

	if mgiID.has_key(key):
	    fp.write(string.join(mgiID[key], reportlib.TAB))
	fp.write(reportlib.CRT)

reportlib.finish_nonps(fp)
