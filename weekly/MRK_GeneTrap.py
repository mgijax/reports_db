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
# lec	05/12/2015
#	- TR12020/use ALL_Allele._Marker_key
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
import db

db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE(False)

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# all official/interim mouse markers that have at least one Gene Trap

db.sql('''
	select m._Marker_key, m.symbol, m.name, m.chromosome, 
        o.cmoffset, 
	upper(substring(s.status, 1, 1)) as markerStatus, 
	t.name as markerType 
        into temporary table markers 
        from MRK_Marker m, MRK_Offset o, MRK_Status s, MRK_Types t 
        where m._Organism_key = 1 
        and m._Marker_Status_key = 1 
        and m._Marker_key = o._Marker_key 
        and o.source = 0 
        and m._Marker_Status_key = s._Marker_Status_key 
        and m._Marker_Type_key = t._Marker_Type_key 
	and exists (select 1 from ALL_Allele a 
	where m._Marker_key = a._Marker_key
	and a.isMixed = 0 
	and a._Allele_Type_key = 847121)
	''', None)
db.sql('create index idx1 on markers(_Marker_key)', None)
db.sql('create index idx2 on markers(symbol)', None)

# MGI ids

results = db.sql('''
	select distinct m._Marker_key, a.accID 
        from markers m, ACC_Accession a 
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
        from markers m, ALL_Allele a, ALL_Allele_CellLine ac, ALL_Cellline c 
        where m._Marker_key = a._Marker_key 
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

results = db.sql('select * from markers order by symbol', 'auto')

for r in results:
	key = r['_Marker_key']

	if r['cmoffset'] == -1.0:
		cmoffset = 'syntenic'
	elif r['cmoffset'] == -999.0:
		cmoffset = 'N/A'
	else:
		cmoffset = str(r['cmoffset'])

	fp.write(mgiID[key] + reportlib.TAB + \
	       	 r['symbol'] + reportlib.TAB + \
	       	 r['markerStatus'] + reportlib.TAB + \
	         r['markerType'] + reportlib.TAB + \
	         r['name'] + reportlib.TAB + \
	         cmoffset + reportlib.TAB + \
	         r['chromosome'] + reportlib.TAB)

	if mutantCellLine.has_key(key):
		fp.write(string.join(mutantCellLine[key], ' '))
	fp.write(reportlib.TAB)

	fp.write(reportlib.CRT)

reportlib.finish_nonps(fp)
