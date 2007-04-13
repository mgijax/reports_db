#!/usr/local/bin/python

'''
#
# MRK_GeneTrap.py 
#
# Report:
#       Tab-delimited file
#       Mouse Genes and their Gene Trap associations
#       TR8112
# Usage:
#       MRK_GeneTrap.py
#
# Used by:
#       Those establishing relationships between MGI Markers
#       and Gene Trap cell line ids
#
# Notes:
#
# History:
#
# sc    04/12/2007 - created
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

# all official/interim mouse markers that have at least one Gene Trap

db.sql('select m._Marker_key, m.symbol, m.name, m.chromosome, ' + \
        'o.offset, markerStatus = upper(substring(s.status, 1, 1)), markerType = t.name ' + \
        'into #markers ' + \
        'from MRK_Marker m, MRK_Offset o, MRK_Status s, MRK_Types t ' + \
        'where m._Organism_key = 1 ' + \
        'and m._Marker_Status_key in (1,3) ' + \
        'and m._Marker_key = o._Marker_key ' + \
        'and o.source = 0 ' + \
        'and m._Marker_Status_key = s._Marker_Status_key ' + \
        'and m._Marker_Type_key = t._Marker_Type_key ' + \
        'and exists (select 1 from ACC_Accession a where m._Marker_key = a._Object_key ' + \
        'and a._MGIType_key = 2 and a._LogicalDB_key in ' + \
        '(select _Object_key ' +
        'from MGI_SetMember ' + \
        'where _Set_key = 1023) )', None)
db.sql('create index idx1 on #markers(_Marker_key)', None)
db.sql('create index idx2 on #markers(symbol)', None)

# MGI ids

results = db.sql('select distinct m._Marker_key, a.accID ' + \
      'from #markers m, ACC_Accession a ' + \
      'where m._Marker_key = a._Object_key ' + \
      'and a._MGIType_key = 2 ' + \
      'and a._LogicalDB_key = 1 ' + \
      'and a.prefixPart = "MGI:" ' + \
      'and a.preferred = 1', 'auto')
mgiID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    mgiID[key] = value

# Gene Trap ids

results = db.sql('select distinct m._Marker_key, a.accID ' + \
      'from #markers m, ACC_Accession a ' + \
      'where m._Marker_key = a._Object_key ' + \
      'and a._MGIType_key = 2 ' + \
      'and a._LogicalDB_key in ' + \
      '(select _Object_key ' +        
        'from MGI_SetMember ' + \
        'where _Set_key = 1023 )', 'auto')
gtID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if not gtID.has_key(key):
        gtID[key] = []
    gtID[key].append(value)

# process

results = db.sql('select * from #markers order by symbol', 'auto')

for r in results:
	key = r['_Marker_key']

	if r['offset'] == -1.0:
		offset = 'syntenic'
	elif r['offset'] == -999.0:
		offset = 'N/A'
	else:
		offset = str(r['offset'])

	fp.write(mgiID[key] + reportlib.TAB + \
	       	 r['symbol'] + reportlib.TAB + \
	       	 r['markerStatus'] + reportlib.TAB + \
	         r['markerType'] + reportlib.TAB + \
	         r['name'] + reportlib.TAB + \
	         offset + reportlib.TAB + \
	         r['chromosome'] + reportlib.TAB)

	if gtID.has_key(key):
		fp.write(string.join(gtID[key], ' '))
	fp.write(reportlib.TAB)

	fp.write(reportlib.CRT)

reportlib.finish_nonps(fp)
