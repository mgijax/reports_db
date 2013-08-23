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
# lec	05/14/2009
#	- TR9405/gene trap less filling
#
# sc    04/12/2007 - created
#
'''

import sys
import os
import string
import mgi_utils
import reportlib
import pg_db
db = pg_db
db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE()

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# all official/interim mouse markers that have at least one Gene Trap

db.sql('''
	select m._Marker_key, m.symbol, m.name, m.chromosome, 
        o.offset, 
	upper(substring(s.status, 1, 1)) as markerStatus, 
	t.name as markerType 
        into #markers 
        from MRK_Marker m, MRK_Offset o, MRK_Status s, MRK_Types t 
        where m._Organism_key = 1 
        and m._Marker_Status_key in (1,3) 
        and m._Marker_key = o._Marker_key 
        and o.source = 0 
        and m._Marker_Status_key = s._Marker_Status_key 
        and m._Marker_Type_key = t._Marker_Type_key 
	and exists (select 1 from ALL_Marker_Assoc am, ALL_Allele a 
	where m._Marker_key = am._Marker_key 
	and am._Allele_key = a._Allele_key 
        and a.isMixed = 0 
	and a._Allele_Type_key = 847121)
	''', None)
db.sql('create index idx1 on #markers(_Marker_key)', None)
db.sql('create index idx2 on #markers(symbol)', None)

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

# Mutant Cell Line for gene traps

results = db.sql('''
	select distinct m._Marker_key, c.cellLine 
        from #markers m, ALL_Marker_Assoc am, ALL_Allele a, ALL_Allele_CellLine ac, ALL_Cellline c 
        where m._Marker_key = am._Marker_key 
        and am._Allele_key  = a._Allele_key 
        and a._Allele_Type_key = 847121 
        and a.isMixed = 0 
        and a._Allele_key = ac._Allele_key 
        and ac._MutantCellLine_key = c._CellLine_key 
	''', 'auto')
#      and c.cellLine not in ("Not Specified", "Not Applicable")', 'auto')
mutantCellLine = {}
for r in results:
    key = r['_Marker_key']
    value = r['cellLine']
    if not mutantCellLine.has_key(key):
        mutantCellLine[key] = []
    mutantCellLine[key].append(value)

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

	if mutantCellLine.has_key(key):
		fp.write(string.join(mutantCellLine[key], ' '))
	fp.write(reportlib.TAB)

	fp.write(reportlib.CRT)

reportlib.finish_nonps(fp)
