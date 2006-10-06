#!/usr/local/bin/python

'''
#
# MRK_GXD.py 06/20/2002
#
# Report:
#       Tab-delimited file of MGI Genes with associated GXD Index records
#
# Usage:
#       MRK_GXD.py
#
# Used by:
#	Deanna Church of NCBI
#
# Notes:
#
# History:
#
# lec	06/20/2002
#	- TR 3798
#
'''
 
import sys
import os
import string
import db
import reportlib

TAB = reportlib.TAB
CRT = reportlib.CRT

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

db.sql('select _Marker_key, _Refs_key into #gxd from GXD_Index', None)
db.sql('create index idx1 on #gxd(_Marker_key)', None)
db.sql('create index idx2 on #gxd(_Refs_key)', None)

results = db.sql('select a.accID, g._Marker_key ' + \
	'from #gxd g, ACC_Accession a ' + \
	'where g._Refs_key = a._Object_key ' + \
	'and a._MGIType_key = 1 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "J:" ' + \
	'and a.preferred = 1', 'auto')

refs = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if not refs.has_key(key):
	refs[key] = []
    refs[key].append(value) 

results = db.sql('select distinct a.accID, g._Marker_key, m.symbol ' + \
	'from #gxd g, MRK_Marker m, ACC_Accession a ' + \
	'where g._Marker_key = m._Marker_key ' + \
	'and m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.preferred = 1', 'auto')

for r in results:
	fp.write(r['accID'] + TAB + r['symbol'] + TAB + string.join(refs[r['_Marker_key']], ',') + CRT)

reportlib.finish_nonps(fp)

