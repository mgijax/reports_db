#!/usr/local/bin/python

'''
#
# MRK_List.py
#
# Report:
#       TR 9120 (original)
#
#	Tab-delimited version of existing report, MRK_List.sql
#	"Genetic Marker List (sorted alphabetically/includes withdrawns)"
#
# Usage:
#       MRK_List.py
#
# History:
#
#
# sc    06/27/2014
#       - TR11560 Feature Relationships project
#       - exclude feature type 'mutation defined region'
#
# lec	04/23/2012
#	- TR11035/Postgres cleanup; merge 
#		MRK_List1, MRK_List2, MRK_List3, MRK_List4,
#		MRK_Synonym, MRK_Dump1, MRK_Dump2 into one report
#
# jer	07/01/2008
#	- created
#
'''
 
import sys 
import os
import string
import mgi_utils
import reportlib
import db

db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE(False)

TAB = reportlib.TAB
CRT = reportlib.CRT

db.useOneConnection(1)
fp1 = reportlib.init('MRK_List1', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp2 = reportlib.init('MRK_List2', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

headers = [
    "MGI Accession ID",
    "Chr", 
    "cM Position", 
    "genome coordinate start",
    "genome coordinate end",
    "strand",
    "Marker Symbol", 
    "Status", 
    "Marker Name", 
    "Marker Type",
    "Feature Type",
    "Marker Synonyms (pipe-separated)"]

fp1.write(TAB.join(headers) + CRT)
fp2.write(TAB.join(headers) + CRT)

#
# select all mouse markers
#
#    and symbol = 'Kit'
db.sql('''
    select m.*, 
    upper(substring(s.status, 1, 1)) as markerstatus,
    mlc.genomicChromosome,
    t.name as markertype,
        case
        when o.cmoffset >= 0 then to_char(o.cmOffset, \'999.99\')
        when o.cmoffset = -999.0 then '       N/A'
        when o.cmoffset = -1.0 then '  syntenic'
        end as cmposition
    into temporary table markers
    from MRK_Marker m, MRK_Status s, MRK_Types t, MRK_Offset o,
    	MRK_Location_Cache mlc
    where m._organism_key = 1
    and m._marker_key = o._marker_key
    and m._marker_key = mlc._marker_key
    and o.source = 0
    and m._marker_type_key = t._marker_type_key
    and m._marker_status_key = s._marker_status_key
    ''', None)
db.sql('create index markers_idx1 on markers(_Marker_key)', None)
db.sql('create index markers_idx2 on markers(symbol)', None)

#
# coordinates
#
results = db.sql('''	
    select m._marker_key,
           c.strand, 
	   c.startCoordinate::int as startC,
	   c.endCoordinate::int as endC
    from markers m, MRK_Location_Cache c
    where m._marker_key = c._marker_key
	''', 'auto')
coords = {}
for r in results:
    key = r['_marker_key']
    value = r
    if not coords.has_key(key):
	coords[key] = []
    coords[key].append(value)

#
# feature types
#
results = db.sql('''
	select m._marker_key, s.term, s._MCVTerm_key
        from markers m, MRK_MCV_Cache s 
        where m._marker_key = s._marker_key 
        and s.qualifier = 'D'
	''', 'auto')
# {markerKey:[list of feature types]
featureTypes = {}
# {markerKey:[list of feature type keys]
featureTypeByKey = {}
for r in results:
        markerKey = r['_marker_key']
        value = r['term']
	termKey = r['_MCVTerm_key']
        if not featureTypes.has_key(markerKey):
                featureTypes[markerKey] = []
        featureTypes[markerKey].append(value)
	if not featureTypeByKey.has_key(markerKey):
                featureTypeByKey[markerKey] = []
        featureTypeByKey[markerKey].append(termKey)
#
# synonyms
#
results = db.sql('''	
    select m._marker_key, s.synonym
    from markers m, MGI_Synonym s, MGI_SynonymType st
    where m._marker_key = s._object_key
    and s._mgitype_key = 2
    and s._synonymtype_key = st._synonymtype_key
    and st.synonymtype = 'exact'
	''', 'auto')
synonyms = {}
for r in results:
    key = r['_marker_key']
    value = r['synonym']
    if not synonyms.has_key(key):
	synonyms[key] = []
    synonyms[key].append(value)

#
# main report
#

query1 = '''
    select a.accid, 
	   m._marker_key, 
	   m.chromosome, 
	   m.genomicChromosome,
	   m.cmposition, 
	   m.symbol, 
	   m.markerstatus, 
	   m.name, 
	   m.markertype
    from markers m, ACC_Accession a
    where m._marker_status_key in (1, 3)
    and m._marker_key = a._object_key
    and a._mgitype_key = 2
    and a.prefixpart = 'MGI:'
    and a._logicaldb_key = 1
    and a.preferred = 1
    '''

query2 = '''
    select 'NULL', 
           m._marker_key, 
           m.chromosome, 
	   m.genomicChromosome,
           m.cmposition, 
           m.symbol, 
           m.markerstatus, 
           m.name, 
           m.markertype
    from markers m
    where m._marker_status_key = 2
    '''

#
# include withdrawns
#
results = db.sql('(%s union %s) order by symbol' % (query1, query2), 'auto')
for r in results:

    key = r['_marker_key']
    # if the marker's feature type is not 
    # 'mutation defined region', key=11928467 write out to the report

    # default feature type
    fTypes = ''
    if featureTypes.has_key(key):
	mcvKeyList = featureTypeByKey[key]
	if 11928467 in mcvKeyList:
	    continue
	else:
	    fTypes = (string.join(featureTypes[key],'|'))

    fp1.write(mgi_utils.prvalue(r['accid']) + TAB)

    if r['genomicChromosome']:
    	fp1.write(r['genomicChromosome'] + TAB)
    else:
    	fp1.write(r['chromosome'] + TAB)

    fp1.write(r['cmposition'] + TAB)

    if coords.has_key(key):
	fp1.write(mgi_utils.prvalue(coords[key][0]['startC']) + TAB)
	fp1.write(mgi_utils.prvalue(coords[key][0]['endC']) + TAB)
	fp1.write(mgi_utils.prvalue(coords[key][0]['strand']) + TAB)
    else:
	fp1.write(TAB + TAB + TAB)

    fp1.write(r['symbol'] + TAB)
    fp1.write(r['markerstatus'] + TAB)
    fp1.write(r['name'] + TAB)
    fp1.write(r['markertype'] + TAB)

    fp1.write(fTypes + TAB)

    if synonyms.has_key(key):
	fp1.write(string.join(synonyms[key],'|'))
    fp1.write(CRT)

#
# do not include withdrawns
#
results = db.sql('%s order by symbol' % (query1), 'auto')
for r in results:

    key = r['_marker_key']

    # if the marker's feature type is not
    # 'mutation defined region', key=11928467 write out to the report

    # default ft for withdrawn markers
    fTypes = ''
    if featureTypes.has_key(key):
        mcvKeyList = featureTypeByKey[key]
        if 11928467 in mcvKeyList:
            continue
        else:
            fTypes = (string.join(featureTypes[key],'|'))

    fp2.write(mgi_utils.prvalue(r['accid']) + TAB)

    if r['genomicChromosome']:
    	fp2.write(r['genomicChromosome'] + TAB)
    else:
    	fp2.write(r['chromosome'] + TAB)

    fp2.write(r['cmposition'] + TAB)

    if coords.has_key(key):
	fp2.write(mgi_utils.prvalue(coords[key][0]['startC']) + TAB)
	fp2.write(mgi_utils.prvalue(coords[key][0]['endC']) + TAB)
	fp2.write(mgi_utils.prvalue(coords[key][0]['strand']) + TAB)
    else:
	fp2.write(TAB + TAB + TAB)

    fp2.write(r['symbol'] + TAB)
    fp2.write(r['markerstatus'] + TAB)
    fp2.write(r['name'] + TAB)
    fp2.write(r['markertype'] + TAB)

    if featureTypes.has_key(key):
	fp2.write(fTypes)
    fp2.write(TAB)

    if synonyms.has_key(key):
	fp2.write(string.join(synonyms[key],'|'))
    fp2.write(CRT)

reportlib.finish_nonps(fp1)	# non-postscript file
reportlib.finish_nonps(fp2)	# non-postscript file
db.useOneConnection(0)

