#!/usr/local/bin/python

'''
#
# MGI_PhenotypicAllele.py 05/04/2004
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
#		MGI ID of Gene
#		Gene Symbol
#		RefSeq ID of Gene
#
# Usage:
#       MGI_PhenotypicAllele.py
#
# Used by:
#
# Notes:
#
# History:
#
# lec	05/04/2004
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

fp.write('#\n#Allele symbols are usually of the form xxx<yy>, where the <> enclose the part of the symbol that is superscripted.\n#\n')

cmds = []

# Retrieve all Approved Alleles
# exclude wild types
# exclude QTLs

cmds.append('select a._Allele_key, a._Marker_key, a.symbol, a.name, marker = m.symbol ' + \
	'into #alleles ' + \
	'from ALL_Allele a, MRK_Marker m ' + \
	'where a._Allele_Status_key = 4 ' + \
	'and a.symbol not like "%<+>" ' + \
	'and a._Marker_key = m._Marker_key ' + \
	'and m._Marker_Type_key != 6')

# Retrieve MGI Accession number for Allele

cmds.append('select s._Allele_key, a.accID ' + \
	'from #alleles s, ACC_Accession a ' + \
	'where s._Allele_key = a._Object_key ' + \
	'and a._MGIType_key = 11 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1')

# Retrieve MGI Accession number for Marker

cmds.append('select s._Marker_key, a.accID ' + \
	'from #alleles s, ACC_Accession a ' + \
	'where s._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1')

# Retrieve PubMed IDs for Original Reference

cmds.append('select s._Allele_key, b.accID ' + \
	'from #alleles s, ALL_Reference r, ACC_Accession b ' + \
	'where s._Allele_key = r._Allele_key ' + \
	'and r._RefsType_key = 1 ' + \
	'and r._Refs_key = b._Object_key ' + \
	'and b._MGIType_key = 1 ' + \
	'and b._LogicalDB_key = 29 ')

# Retrieve RefSeq ID for Gene

cmds.append('select s._Marker_key, a.accID ' + \
	'from #alleles s, MRK_Acc_View a ' + \
	'where s._Marker_key = a._Object_key ' + \
	'and a._LogicalDB_key = 27 ')

cmds.append('select * from #alleles order by marker, symbol')

results = db.sql(cmds, 'auto')

amgiIDs = {}
for r in results[-5]:
	amgiIDs[r['_Allele_key']] = r['accID']
	
mmgiIDs = {}
for r in results[-4]:
	mmgiIDs[r['_Marker_key']] = r['accID']

pubIDs = {}
for r in results[-3]:
	pubIDs[r['_Allele_key']] = r['accID']
	
refIDs = {}
for r in results[-2]:
	refIDs[r['_Marker_key']] = r['accID']
	
for r in results[-1]:
	fp.write(amgiIDs[r['_Allele_key']] + reportlib.TAB + \
		 r['symbol'] + reportlib.TAB + \
		 r['name'] + reportlib.TAB)

	if pubIDs.has_key(r['_Allele_key']):
		fp.write(pubIDs[r['_Allele_key']])
	fp.write(reportlib.TAB)

	# if Transgene, do not print gene ID or gene symbol
	if string.find(r['symbol'], 'Tg(') < 0:
		fp.write(mmgiIDs[r['_Marker_key']] + reportlib.TAB + \
			r['marker'] + reportlib.TAB)
	else:
		fp.write(reportlib.TAB + reportlib.TAB)

	if refIDs.has_key(r['_Marker_key']):
		fp.write(refIDs[r['_Marker_key']])
	fp.write(reportlib.TAB)

	fp.write(reportlib.CRT)

reportlib.finish_nonps(fp)

