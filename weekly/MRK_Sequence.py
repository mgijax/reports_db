#!/usr/local/bin/python

'''
#
# MRK_Sequence.py 03/02/99
#
# Report:
#       Tab-delimited file
#       Mouse Markers and their Nucleotide Sequence Accession numbers.
#	and their UniGene Accession numbers (TR 1631).
#
# Usage:
#       MRK_Sequence.py
#
# Used by:
#       Those establishing relationships between MGI Markers
#	and Nucleotide Sequences.
#
# Notes:
#
# History:
#
# lec	01/27/2005
#	- TR 6529
#
# lec	12/18/2003
#	- TR 5440; added Marker Type
#
# lec	05/30/1999
#	- TR 1631; add UniGene Accession numbers
#
# lec	03/02/1999
#	- TR 130; now use direct Marker-Seq ID relationships
#
# lec	01/18/1999
#	- missing Segment records where relationship is null
#
# lec	01/13/98
#	- added comments section
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

# all official/interim mouse markers that have at least one Sequence ID

cmds = []
cmds.append('select m._Marker_key, m.symbol, m.name, m.chromosome, ' + \
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
	'and a._MGIType_key = 2 and a._LogicalDB_key in (9, 27) and a.prefixPart not in ("XP_", "NP_"))')
cmds.append('create index idx1 on #markers(_Marker_key)')
cmds.append('create index idx2 on #markers(symbol)')
db.sql(cmds, None)

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

# GenBank ids

results = db.sql('select distinct m._Marker_key, a.accID ' + \
      'from #markers m, ACC_Accession a ' + \
      'where m._Marker_key = a._Object_key ' + \
      'and a._MGIType_key = 2 ' + \
      'and a._LogicalDB_key = 9', 'auto')
gbID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if not gbID.has_key(key):
	gbID[key] = []
    gbID[key].append(value)

# UniGene ids

results = db.sql('select distinct m._Marker_key, a.accID ' + \
      'from #markers m, ACC_Accession a ' + \
      'where m._Marker_key = a._Object_key ' + \
      'and a._MGIType_key = 2 ' + \
      'and a._LogicalDB_key = 23', 'auto')
ugID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if not ugID.has_key(key):
	ugID[key] = []
    ugID[key].append(value)

# RefSeq ids

results = db.sql('select distinct m._Marker_key, a.accID ' + \
      'from #markers m, ACC_Accession a ' + \
      'where m._Marker_key = a._Object_key ' + \
      'and a._MGIType_key = 2 ' + \
      'and a._LogicalDB_key = 27 ' + \
      'and a.prefixPart not in ("XP_", "NP_")', 'auto')
rsID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if not rsID.has_key(key):
	rsID[key] = []
    rsID[key].append(value)

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

	if gbID.has_key(key):
		fp.write(string.join(gbID[key], ' '))
	fp.write(reportlib.TAB)

	if ugID.has_key(key):
		fp.write(string.join(ugID[key], ' '))
	fp.write(reportlib.TAB)

	if rsID.has_key(key):
		fp.write(string.join(rsID[key], ' '))
	fp.write(reportlib.CRT)

reportlib.finish_nonps(fp)
