#!/usr/local/bin/python

'''
#
# HMD_HumanDisease.py
#
# Report:
#
#     This report was originally a custom SQL request (TR 8437) for
#     Stephen Campbell (Stephen.J.Campbell@pfizer.com).
# 
#     Output fields (tab-delimited):
# 
#     1) Allele ID
#     2) Allele Symbol
#     3) Entrez Gene ID
#     4) PubMed ID
#     5) OMIM ID
#     6) OMIM Term
#     7) MP ID
#     8) MP Term
# 
# Usage:
#     HMD_HumanDisease.py
#
# History:
#
# dbm	8/31/2007
#	-  Create from custom SQL request written by smc (TR 8437)
#
'''
 
import os
import sys
import reportlib
import db
import string


CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# select all OMIM annotations
#

db.sql('select distinct a._Object_key, a._Term_key, a.term, a.qualifier, e._Refs_key ' + 
	'into #omim ' + \
	'from VOC_Annot_View a, VOC_Evidence e ' + \
	'where a._AnnotType_key = 1005 ' + \
	'and a._Annot_key = e._Annot_key', None)
db.sql('create index idx1 on #omim(_Object_key)', None)
db.sql('create index idx2 on #omim(_Term_key)', None)
db.sql('create index idx3 on #omim(_Refs_key)', None)

#
# resolve OMIM ids
#
results = db.sql('select distinct m._Term_key, a.accID from #omim m, ACC_Accession a ' + \
	'where m._Term_key = a._Object_key ' + \
	'and a._MGIType_key = 13 ' + \
	'and a.preferred = 1', 'auto')
omimID = {}
for r in results:
	key = r['_Term_key']
	value = r['accID']
	omimID[key] = value

#
# resolve MP ids
#
results = db.sql('select distinct m._Object_key, a.accID ' + \
	'from #omim m, VOC_Annot_View a ' + \
	'where m._Object_key = a._Object_key ' + \
	'and a._AnnotType_key = 1002', 'auto')
MPID = {}
for r in results:
	key = r['_Object_key']
	value = r['accID']
	if not MPID.has_key(key):
		MPID[key] = []
	MPID[key].append(value)

#
# resolve MP terms
#
results = db.sql('select distinct m._Object_key, a.term ' + \
	'from #omim m, VOC_Annot_View a ' + \
	'where m._Object_key = a._Object_key ' + \
	'and a._AnnotType_key = 1002', 'auto')
MPTerm = {}
for r in results:
	key = r['_Object_key']
	value = r['term']
	if not MPTerm.has_key(key):
		MPTerm[key] = []
	MPTerm[key].append(value)

#
# resolve Pubmed References for OMIM annotations
#
results = db.sql('select distinct m._Object_key, m._Term_key, a.accID ' + \
	'from #omim m, ACC_Accession a ' + \
	'where m._Refs_key = a._Object_key ' + \
	'and a._LogicalDB_key = 29 ', 'auto')
	
omimRef = {}
for r in results:
	key = `r['_Object_key']` + `r['_Term_key']`
	value = r['accID']
	if not omimRef.has_key(key):
		omimRef[key] = []
	omimRef[key].append(value)


#
# resolve allele IDs
#
results = db.sql('select distinct m._Object_key, a.accID ' + \
	'from #omim m, GXD_AlleleGenotype g, ACC_Accession a, ALL_Allele l ' + \
	'where m._Object_key = g._Genotype_key ' + \
	'and g._Allele_key = a._Object_key ' + \
	'and a._MGIType_key = 11 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1 ' + \
	'and g._Allele_key = l._Allele_key ' + \
	'and l.isWildType = 0', 'auto')
alleleID = {}
for r in results:
	key = r['_Object_key']
	value = r['accID']
	alleleID[key] = value

#
# resolve allele symbols
#
results = db.sql('select distinct m._Object_key, a.symbol ' + \
	'from #omim m, GXD_AlleleGenotype g, ALL_Allele a ' + \
	'where m._Object_key = g._Genotype_key ' + \
	'and g._Allele_key = a._Allele_key ' + \
	'and a.isWildType = 0', 'auto')
symbol = {}
for r in results:
	key = r['_Object_key']
	value = r['symbol']
	symbol[key] = value

#
# resolve Marker ID for Entrez Gene ID
#
results = db.sql('select distinct m._Object_key, a.accID ' + \
	'from #omim m, GXD_AlleleGenotype g, ACC_Accession a ' + \
	'where m._Object_key = g._Genotype_key ' + \
	'and g._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 55', 'auto')
omimMarker = {}
for r in results:
	key = r['_Object_key']
	value = r['accID']
	if not omimMarker.has_key(key):
		omimMarker[key] = []
	omimMarker[key].append(value)

#
# process results
#
results = db.sql('select * from #omim order by term, _Object_key', 'auto')

fp.write('Allele ID' + TAB + 'Allele Symbol' + TAB + 'Entrez Gene ID' + TAB + 'PubMed ID' + TAB)
fp.write('OMIM ID' + TAB + 'OMIM Term' + TAB + 'MP ID'+ TAB + 'MP Term' + CRT + CRT)

for r in results:

	genotype = r['_Object_key']
	term = r['_Term_key']
	refKey = `genotype` + `term`


	if alleleID.has_key(genotype):
		fp.write(alleleID[genotype] + TAB + \
		symbol[genotype])
		if omimMarker.has_key(genotype):
			fp.write(TAB + string.join(omimMarker[genotype], ',') + TAB)
		else:
			fp.write(TAB + TAB)
		if omimRef.has_key(refKey):
			fp.write(string.join(omimRef[refKey], ','))
		fp.write(TAB + omimID[term] + TAB + r['term'] + TAB)
		if MPID.has_key(genotype):
			fp.write(string.join(MPID[genotype], ','))
		fp.write(TAB)	
		if MPTerm.has_key(genotype):
			fp.write(string.join(MPTerm[genotype], ','))
		fp.write(CRT)

reportlib.finish_nonps(fp)	# non-postscript file

