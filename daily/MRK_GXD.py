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

cmd = 'select distinct a.accID, m.symbol ' + \
	'from MRK_Marker m, GXD_Index g, ACC_Accession a ' + \
	'where m._Marker_key = g._Marker_key ' + \
	'and m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.preferred = 1'

results = db.sql(cmd, 'auto')

for r in results:
	fp.write(r['accID'] + TAB + r['symbol'] + CRT)

reportlib.finish_nonps(fp)

