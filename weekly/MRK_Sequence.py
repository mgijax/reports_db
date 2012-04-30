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
# sc	03/12/2010
#	- TR9774 Add Ensembl and VEGA transcripts
#
# lec	10/10/2006
#	- only include Markers that have at least one sequence.
#
# lec	09/19/2006
#	- make sure there are no deleted sequences in the report
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
import mgi_utils
import reportlib

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db


#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# deleted sequences

db.sql('select s._Sequence_key into #deleted from SEQ_Sequence s where s._SequenceStatus_key = 316343', None)
db.sql('create index deleted_idx1 on #deleted(_Sequence_key)', None)

db.sql('''
    select a.accID, a._LogicalDB_key 
    into #deletedIDs from #deleted d, ACC_Accession a 
    where d._Sequence_key = a._Object_key 
    and a._MGIType_key = 19
    ''', None)
db.sql('create index deletedIDS_idx1 on #deletedIDs(accID)', None)
db.sql('create index deletedIDS_idx2 on #deletedIDs(_LogicalDB_key)', None)

# all official/interim mouse markers that have at least one Sequence ID

db.sql('''
	select m._Marker_key, m.symbol, m.name, m.chromosome, 
	o.offset, upper(substring(s.status, 1, 1)) as markerStatus, t.name as markerType
	into #markers 
	from MRK_Marker m, MRK_Offset o, MRK_Status s, MRK_Types t 
	where m._Organism_key = 1 
	and m._Marker_Status_key in (1,3) 
	and m._Marker_key = o._Marker_key 
	and o.source = 0 
	and m._Marker_Status_key = s._Marker_Status_key 
	and m._Marker_Type_key = t._Marker_Type_key 
	and exists (select 1 from ACC_Accession a where m._Marker_key = a._Object_key 
	and a._MGIType_key = 2 and a._LogicalDB_key in (9, 27, 131, 133) and a.prefixPart not in ("XP_", "NP_"))
	''', None)
db.sql('create index markers_idx1 on #markers(_Marker_key)', None)
db.sql('create index markers_idx2 on #markers(symbol)', None)

# MGI ids

results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from #markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 1 
      and a.prefixPart = 'MGI:' 
      and a.preferred = 1
      ''', 'auto')
mgiID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    mgiID[key] = value

# GenBank ids

results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from #markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 9 
      and not exists (select 1 from #deletedIDs d where a.accID = d.accID and a._LogicalDB_key = d._LogicalDB_key)
      ''', 'auto')
gbID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if not gbID.has_key(key):
	gbID[key] = []
    gbID[key].append(value)

# UniGene ids

results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from #markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 23 
      and not exists (select 1 from #deletedIDs d where a.accID = d.accID and a._LogicalDB_key = d._LogicalDB_key)
      ''', 'auto')
ugID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if not ugID.has_key(key):
	ugID[key] = []
    ugID[key].append(value)

# RefSeq ids

results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from #markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 27 
      and a.prefixPart not in ("XP_", "NP_") 
      and not exists (select 1 from #deletedIDs d where a.accID = d.accID and a._LogicalDB_key = d._LogicalDB_key)
      ''', 'auto')
rsID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if not rsID.has_key(key):
	rsID[key] = []
    rsID[key].append(value)

# Ensembl transript IDs
results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from #markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 133 
      and not exists (select 1 from #deletedIDs d 
      where a.accID = d.accID 
      and a._LogicalDB_key = d._LogicalDB_key)
      ''', 'auto')
ensID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if not ensID.has_key(key):
        ensID[key] = []
    ensID[key].append(value)

# VEGA transcript IDs
results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from #markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 131 
      and not exists (select 1 from #deletedIDs d 
      where a.accID = d.accID 
      and a._LogicalDB_key = d._LogicalDB_key)
      ''', 'auto')
vegaID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if not vegaID.has_key(key):
        vegaID[key] = []
    vegaID[key].append(value)

# process

results = db.sql('select * from #markers order by symbol', 'auto')

for r in results:
	key = r['_Marker_key']
	symbol = r['symbol']

	if not gbID.has_key(key) and not ugID.has_key(key) and not rsID.has_key(key):
	    print 'not gb', symbol
	    continue

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
	fp.write(reportlib.TAB)

	if vegaID.has_key(key):
		fp.write(string.join(vegaID[key], ' '))
	fp.write(reportlib.TAB)

        if ensID.has_key(key):
                fp.write(string.join(ensID[key], ' '))
	fp.write(reportlib.CRT)

reportlib.finish_nonps(fp)

