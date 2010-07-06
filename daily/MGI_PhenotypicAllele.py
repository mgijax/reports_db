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
#		Allele Type
#		PubMed ID of Original Reference
#		MGI ID of Gene
#		Gene Symbol
#		RefSeq ID of Gene
#		Ensembl ID
#		MP IDs of MP annotations
#		Synonyms
#
# Usage:
#       MGI_PhenotypicAllele.py
#
# History:
#
# lec	07/06/2010
#	- TR 10280/add synonyms
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

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('#\n#Allele symbols are usually of the form xxx<yy>, where the <> enclose the part of the symbol that is superscripted.\n#\n')
fp.write('#Transgene insertions, with symbols of the form Tg(aaa)##bbb, are included in this listing, but notably have no corresponding gene marker.\n#\n')
fp.write('#For details of nomenclature rules, see http://www.informatics.jax.org/mgihome/nomen/index.shtml\n#\n')

# Retrieve all Approved Alleles
# exclude wild types
# exclude QTLs

db.sql('select a._Allele_key, a._Marker_key, a.symbol, a.name, alleleType = t2.term, marker = m.symbol ' + \
	'into #alleles ' + \
	'from ALL_Allele a, VOC_Term t1, VOC_Term t2, MRK_Marker m ' + \
	'where a._Allele_Status_key = t1._Term_key ' + \
	'and t1.term = "Approved" ' + \
	'and a.isWildType = 0 ' + \
	'and a._Allele_Type_key = t2._Term_key ' + \
	'and a._Marker_key = m._Marker_key ' + \
	'and m._Marker_Type_key != 6', None)
db.sql('create index idx1 on #alleles(_Allele_key)', None)

# Retrieve MGI Accession number for Allele

results = db.sql('select s._Allele_key, a.accID ' + \
	'from #alleles s, ACC_Accession a ' + \
	'where s._Allele_key = a._Object_key ' + \
	'and a._MGIType_key = 11 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1', 'auto')
amgiIDs = {}
for r in results:
	amgiIDs[r['_Allele_key']] = r['accID']
	

# Retrieve MGI Accession number for Marker

results = db.sql('select s._Marker_key, a.accID ' + \
	'from #alleles s, ACC_Accession a ' + \
	'where s._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1', 'auto')
mmgiIDs = {}
for r in results:
	mmgiIDs[r['_Marker_key']] = r['accID']

# Retrieve PubMed IDs for Original Reference

results = db.sql('select s._Allele_key, b.accID ' + \
	'from #alleles s, MGI_Reference_Assoc r, MGI_RefAssocType rt, ACC_Accession b ' + \
	'where s._Allele_key = r._Object_key ' + \
	'and r._MGIType_key = 11 ' + \
	'and r._RefAssocType_key = rt._RefAssocType_key ' + \
	'and rt.assocType = "Original" ' + \
	'and r._Refs_key = b._Object_key ' + \
	'and b._MGIType_key = 1 ' + \
	'and b._LogicalDB_key = 29 ', 'auto')
pubIDs = {}
for r in results:
	pubIDs[r['_Allele_key']] = r['accID']

# Retrieve RefSeq ID for Gene

results = db.sql('select s._Marker_key, a.accID ' + \
	'from #alleles s, ACC_Accession a ' + \
	'where s._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 27 ' + \
	'and a.prefixPart in ("NM_", "XM_")', 'auto')
refIDs = {}
for r in results:
	refIDs[r['_Marker_key']] = r['accID']
	
# Retrieve Ensembl Gene Model ID for Gene

results = db.sql('select s._Marker_key, a.accID ' + \
	'from #alleles s, ACC_Accession a ' + \
	'where s._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 60 ', 'auto')
ensemblIDs = {}
for r in results:
	ensemblIDs[r['_Marker_key']] = r['accID']
	
# Retrieve MP IDs for MP annotations

results = db.sql('select distinct s._Allele_key, a.accID ' + \
	'from #alleles s, GXD_AlleleGenotype ga, VOC_AnnotHeader na, ACC_Accession a ' + \
	'where s._Allele_key = ga._Allele_key ' + \
	'and ga._Genotype_key = na._Object_key ' + \
	'and na._AnnotType_key = 1002 ' + \
	'and na._Term_key = a._Object_key ' + \
	'and a._MGIType_key = 13 ' + \
	'and a.preferred = 1', 'auto')

mpIDs = {}
for r in results:
	if not mpIDs.has_key(r['_Allele_key']):
		mpIDs[r['_Allele_key']] = []
	mpIDs[r['_Allele_key']].append(r['accID'])
	
# Retrieve Synonyms
results = db.sql('select s._Allele_key, ss.synonym ' + \
        'from #alleles s, MGI_Synonym ss, MGI_SynonymType st ' + \
        'where s._Allele_key = ss._Object_key ' + \
        'and ss._MGIType_key = 11 ' + \
        'and ss._SynonymType_key = st._SynonymType_key ', 'auto')
        #'and st.synonymType = "exact"', 'auto')
synonym = {}
for r in results:
        key = r['_Allele_key']
        value = r['synonym']
        if not synonym.has_key(key):
                synonym[key] = []
        synonym[key].append(value)

#
# Main
#

results = db.sql('select * from #alleles order by marker, symbol', 'auto')

for r in results:

	fp.write(amgiIDs[r['_Allele_key']] + reportlib.TAB + \
		 r['symbol'] + reportlib.TAB + \
		 r['name'] + reportlib.TAB + \
		 r['alleleType'] + reportlib.TAB)

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

	if ensemblIDs.has_key(r['_Marker_key']):
		fp.write(ensemblIDs[r['_Marker_key']])
	fp.write(reportlib.TAB)

	if mpIDs.has_key(r['_Allele_key']):
		fp.write(string.join(mpIDs[r['_Allele_key']], ','))
	fp.write(reportlib.TAB)

	if synonym.has_key(r['_Allele_key']):
		fp.write(string.join(synonym[r['_Allele_key']], '|'))
	fp.write(reportlib.CRT)

reportlib.finish_nonps(fp)

