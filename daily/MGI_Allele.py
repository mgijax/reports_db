#!/usr/local/bin/python

'''
#
# MGI_Allele.py 04/01/2004
#
# Report:
#       Tab-delimited file
#       Allele Nomenclature
#
#	Fields:
#		MGI ID of Allele
#		Allele Symbol
#		Allele Name
#		PubMed ID of Original Reference
#		Gene Symbol
#		RefSeq ID of Gene
#
# Usage:
#       MGI_Allele.py
#
# Used by:
#
# Notes:
#
# History:
#
# lec	04/01/2004
#	- TR 5637
#
'''
 
import sys
import os
import string
import db
import reportlib

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []

# Retrieve all Approved Alleles

cmds.append('select a._Allele_key, a._Marker_key, a.symbol, a.name, marker = m.symbol ' + \
	'into #alleles ' + \
	'from ALL_Allele a, MRK_Marker m ' + \
	'where a._Allele_Status_key = 4 ' + \
	'and a._Marker_key = m._Marker_key')

# Retrieve MGI Accession number

cmds.append('select s._Allele_key, a.accID ' + \
	'from #alleles s, ALL_Acc_View a ' + \
	'where s._Allele_key = a._Object_key ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1')

# Retrieve PubMed IDs for Original Reference

cmds.append('select s._Allele_key, b.accID ' + \
	'from #alleles s, ALL_Reference r, BIB_Acc_View b ' + \
	'where s._Allele_key = r._Allele_key ' + \
	'and r._RefsType_key = 1 ' + \
	'and r._Refs_key = b._Object_key ' + \
	'and b._LogicalDB_key = 29 ')

# Retrieve RefSeq ID for Gene

cmds.append('select s._Marker_key, a.accID ' + \
	'from #alleles s, MRK_Acc_View a ' + \
	'where s._Marker_key = a._Object_key ' + \
	'and a._LogicalDB_key = 27 ')

cmds.append('select * from #alleles order by marker, symbol')

results = db.sql(cmds, 'auto')

mgiIDs = {}
for r in results[-4]:
	mgiIDs[r['_Allele_key']] = r['accID']
	
pubIDs = {}
for r in results[-3]:
	pubIDs[r['_Allele_key']] = r['accID']
	
refIDs = {}
for r in results[-2]:
	refIDs[r['_Marker_key']] = r['accID']
	
for r in results[-1]:
	fp.write(mgiIDs[r['_Allele_key']] + reportlib.TAB + \
		 r['symbol'] + reportlib.TAB + \
		 r['name'] + reportlib.TAB)

	if pubIDs.has_key(r['_Allele_key']):
		fp.write(pubIDs[r['_Allele_key']])
	fp.write(reportlib.TAB)

	fp.write(r['marker'] + reportlib.TAB)

	if refIDs.has_key(r['_Marker_key']):
		fp.write(refIDs[r['_Marker_key']])
	fp.write(reportlib.TAB)

	fp.write(reportlib.CRT)

reportlib.finish_nonps(fp)

