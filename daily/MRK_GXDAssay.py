#!/usr/local/bin/python

'''
#
# MRK_GXDAssay.py 08/28/2002
#
# Report:
#       Tab-delimited file of MGI Genes with associated GXD Assay records
#
# Usage:
#       MRK_GXDAssay.py
#
# Used by:
#	Deanna Church of NCBI
#
# Notes:
#
# History:
#
# lec	08/28/2002
#	- TR 4027
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
	'from MRK_Marker m, GXD_Assay g, ACC_Accession a ' + \
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

