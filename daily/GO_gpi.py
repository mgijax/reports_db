#!/usr/local/bin/python

'''
#
# GO_gpi.py
#
# Report:
#       see: wiki/mediawiki/index.php/sw:GOload 
#	contains link to GO/GPI format
#
# Output format:
#
#   1. DB_Object_ID
#   2. DB_Object_Symbol                   
#   3. DB_Object_Name          (optional)
#   4. DB_Object_Synonym(s)    (optional)
#   5. DB_ObjecT_Type
#   6. Taxon                              tax:10090
#   7. Parent_Object_ID        (optional) 
#	if DB_Object_ID = Isoform, then MGI:id of the Isoform
#   8. DB_Xref(s)              (optional) 
#	if DB_Object_ID = Isoform, then UniProtKB:id of the Isoform
#   9. Gene_Product_Properties (optional)
#
# History:
#
# 06/27/2016	lec
#	- TR12349/12345/GPAD/GPI
#
'''

import sys
import os
import mgi_utils
import reportlib
import db

db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE(False)

SPECIES = 'taxon:10090'
DBTYPE_MARKER = 'gene'
DBTYPE_ISOFORM = 'protein'
DBTYPE_RNA = 'transcript'

TAB = reportlib.TAB
CRT = reportlib.CRT

db.useOneConnection(1)

fp = reportlib.init('mgi_association', fileExt = '.gpi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('!gpi-version: 1.0\n')
fp.write('!date: %s $\n' % (mgi_utils.date("%m/%d/%Y")))
fp.write('!\n')
fp.write('! from Mouse Genome Database (MGD) & Gene Expression Database (GXD)\n')
fp.write('!\n')
fp.write('! 1. DB_Object_ID\n')
fp.write('!     MGI:MGI:xxxxx (gene), PR:xxxx (protein), RNA (transcript)\n')
fp.write('! 2. DB_Object_Symbol\n')   
fp.write('! 3. DB_Object_Name\n')
fp.write('! 4. DB_Object_Synonym(s)\n')
fp.write('! 5. DB_ObjecT_Type\n')
fp.write('! 6. Taxon (tax:10090)\n')
fp.write('! 7. Parent_Object_ID\n')
fp.write('! 	if DB_Object_ID = Isoform, then MGI:id of the Isoform\n')
fp.write('! 8. DB_Xref(s)\n')
fp.write('!     if DB_Object_ID = Isoform, then UniProtKB:id of the Isoform\n')
fp.write('! 9. Gene_Product_Properties\n')
fp.write('!\n')

#
# markers
#

#
# markers:marker->isoform
#

# isoform == 1 marker
isoformByMarker = {}

# marker >= 1 isoform
markerByIsoform = {}

results = db.sql('''
	select 'MGI:' || a1.accID as markerID, a2.accID as prID
	from VOC_Annot v, ACC_Accession a1, ACC_Accession a2
	where v._AnnotType_key = 1019
	and v._Object_key = a1._Object_key
	and a1._MGIType_key = 2
	and a1._LogicalDB_key = 1
	and a1.prefixPart = 'MGI:'
	and a1.preferred = 1
	and v._Term_key = a2._Object_key
	and a2._MGIType_key = 13
   	''', 'auto')

for r in results:

   	byMarker = r['markerID']
	byIsoform = r['prID']

	if byMarker not in isoformByMarker:
		isoformByMarker[byMarker] = []
	isoformByMarker[byMarker].append(byIsoform)

	if byIsoform not in markerByIsoform:
		markerByIsoform[byIsoform] = []
	markerByIsoform[byIsoform].append(byMarker)

#
# markers:synonyms
#

markerSynonyms = {}

results = db.sql('''
	select 'MGI:' || a.accID as accID, s.synonym
	from ACC_Accession a, MRK_Marker m, MGI_Synonym s
	where a._MGIType_key = 2
	and a._LogicalDB_key = 1
	and a.prefixPart = 'MGI:'
	and a.preferred = 1
	and a._Object_key = m._Marker_key
	and m._Marker_Type_key = 1
	and m._Marker_Status_key in (1,3)
	and m._Organism_key = 1
	and m._Marker_key = s._Object_key
	and s._MGIType_key = 2
	and s._SynonymType_key = 1004
   	''', 'auto')
for r in results:
	key = r['accID']
	value = r['synonym']
	if key not in markerSynonyms:
		markerSynonyms[key] = []
	markerSynonyms[key].append(value)

#
# markers:terms
#

results = db.sql('''
	select 'MGI:' || a.accID as accID, m.symbol, m.name
	from ACC_Accession a, MRK_Marker m
	where a._MGIType_key = 2
	and a._LogicalDB_key = 1
	and a.prefixPart = 'MGI:'
	and a.preferred = 1
	and a._Object_key = m._Marker_key
	and m._Marker_Type_key = 1
	and m._Marker_Status_key in (1,3)
	and m._Organism_key = 1
   	''', 'auto')

for r in results:
	marker = r['accID']

	fp.write(r['accID'] + TAB)
	fp.write(r['symbol'] + TAB)
	fp.write(r['name'] + TAB)

	if marker in markerSynonyms:
		fp.write("|".join(markerSynonyms[marker]))
	fp.write(TAB)

	fp.write(DBTYPE_MARKER + TAB)
	fp.write(SPECIES + TAB)
	fp.write(TAB)
	fp.write(TAB)
	fp.write(CRT)

#
# end markers
#

#
# isoforms
#

#
# isoforms:synonyms
#

isoformSynonyms = {}

results = db.sql('''
	select a.accID, s.synonym
	from ACC_Accession a, VOC_Term t, MGI_Synonym s
	where t._Vocab_key = 112
	and t._term_key = a._Object_key and a._LogicalDB_key = 183
	and a._LogicalDB_key = 183
	and t._Term_key = s._Object_key and s._MGIType_key = 13
   	''', 'auto')
for r in results:
	key = r['accID']
	value = r['synonym']
	if key not in isoformSynonyms:
		isoformSynonyms[key] = []
	isoformSynonyms[key].append(value)

#
# isoforms:uniprot ids
#

isoformUniProt = {}

results = db.sql('''
	select a1.accID as prID, p.value
	from VOC_Annot v, ACC_Accession a1, VOC_Evidence ve, VOC_Evidence_Property p
	where v._AnnotType_key = 1019
	and v._Term_key = a1._Object_key
	and a1._MGIType_key = 13
	and v._Annot_key = ve._Annot_key
	and ve._AnnotEvidence_key = p._AnnotEvidence_key
   	''', 'auto')
for r in results:
	key = r['prID']
	value = r['value']
	if key not in isoformUniProt:
		isoformUniProt[key] = []
	isoformUniProt[key].append(value)

#
# isoforms:terms
#

results = db.sql('''
	select a.accID, t.term as symbol, s.note as name
	from ACC_Accession a, VOC_Term t, VOC_Text s
	where t._Vocab_key = 112
	and t._term_key = a._Object_key 
	and a._LogicalDB_key = 183
	and t._Term_key = s._Term_key
   	''', 'auto')

for r in results:
	isoform = r['accID']

	fp.write(isoform + TAB)
	fp.write(r['symbol'] + TAB)
	fp.write(r['name'] + TAB)

	if isoform in isoformSynonyms:
		fp.write("|".join(isoformSynonyms[isoform]))
	fp.write(TAB)

	fp.write(DBTYPE_ISOFORM + TAB)
	fp.write(SPECIES + TAB)

	# Parent_Object_ID/MGI id
	if isoform in markerByIsoform:
		fp.write(markerByIsoform[isoform][0])
	fp.write(TAB)

	# DB_Xref/UniProtKB id
	if isoform in isoformUniProt:
		fp.write(isoformUniProt[isoform][0])
	fp.write(TAB)

	fp.write(CRT)

#
# end isoforms
#

#
# RNAs that are associated with at most 1 marker
#

results = db.sql('''
	WITH rna AS
	(
	select distinct sm.accID as rnaID, sm._Marker_key, sm._LogicalDB_key
        	from SEQ_Marker_Cache sm
        	where sm._SequenceType_key = 316346
        	and sm._Marker_Type_key = 1 
        	and sm._Organism_key = 1 
	)
	select rna.rnaID, a.accID as markerID, rna._LogicalDB_key
	from rna, ACC_Accession a
	where rnaID in (select rnaID from rna group by rnaID having count(*) = 1)
        and rna._Marker_key = a._Object_key
        and a._MGIType_key = 2
        and a._LogicalDB_key = 1
        and a.preferred = 1
	order by rna._LogicalDB_key
   	''', 'auto')

for r in results:

	ldb = r['_LogicalDB_key']
	if ldb == 9:
		ldbName = 'EMBL:'
	elif ldb == 27:
		ldbName = 'RefSeq:'
	elif ldb == 133:
		ldbName = 'ENSEMBL:'
	elif ldb == 131:
		ldbName = 'VEGA:'
	else:
		ldbName = ''

	fp.write(ldbName + r['rnaID'] + TAB)
	fp.write(r['rnaID'] + TAB)
	fp.write(TAB)
	fp.write(TAB)
	fp.write(DBTYPE_RNA + TAB)
	fp.write(SPECIES + TAB)
	fp.write(r['markerID'] + TAB)
	fp.write(TAB)
	fp.write(CRT)

# end RNAs
#

reportlib.finish_nonps(fp)
db.useOneConnection(0)

