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

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []

cmds.append('select distinct _Marker_key into #gxd from GXD_Index')
cmds.append('create index idx1 on #gxd(_Marker_key)')
cmds.append('select a.accID, m.symbol ' + \
	'from #gxd g, MRK_Marker m, ACC_Accession a ' + \
	'where g._Marker_key = m._Marker_key ' + \
	'and m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.preferred = 1')

results = db.sql(cmds, 'auto')

for r in results[-1]:
	fp.write(r['accID'] + TAB + r['symbol'] + CRT)

reportlib.finish_nonps(fp)

