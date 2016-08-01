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
# version 1.2
#
#   1.  DB                     required  1             1             MGI
#   2.  DB_Object_ID           required  1             2/17          MGI:87870
#   3.  DB_Object_Symbol       required  1             3             Acat1
#   4.  DB_Object_Name         optional  0 or greater  10            acetyl-Coenzyme A acetyltransferase 1
#   5.  DB_Object_Synonym(s)   optional  0 or greater  11            Acat|6330585C21Rik
#   6.  DB_Object_Type         required  1             12            gene
#   7.  Taxon                  required  1             13            taxon:10090
#   8.  Parent_Object_ID       optional  0 or 1        -             if DB_Object_ID = Isoform, then MGI:id of the Isoform
#   9.  DB_Xref(s)             optional  0 or greater  -             if DB_Object_ID = Isoform, then UniProtKB:id of the Isoform
#   10. Properties             optional  0 or greater  -             blank
#
#   DB_Object_Type = 'gene', DB = 'MGI', DB_Object_ID = 'MGI:xxxx'
#   DB_Object_Type = 'protein', DB = 'PR', DB_Object_ID = 'xxxx', Parent_Object_ID = 'MGI:MGI:xxxx', DB_Xref = 'UniProtDB:xxx'
#   DB_Object_Type = 'transcript', DB = 'EMBL', 'RefSeq', 'ENSEMBL', 'VEGA', Parent_Object_ID = 'MGI:MGI:xxxx'
#
# History:
#
# 06/27/2016	lec
#	- TR12349/12345/GPAD/GPI
#
'''

import sys
import os
import gzip
import mgi_utils
import reportlib
import db

db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE(False)

DB_PREFIX = 'MGI'
SPECIES = 'taxon:10090'
DBTYPE_MARKER = 'gene'
DBTYPE_ISOFORM = 'protein'
DBTYPE_RNA = 'transcript'

TAB = reportlib.TAB
CRT = reportlib.CRT

db.useOneConnection(1)

fp = reportlib.init('mgi_association', fileExt = '.gpi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('!gpi-version: 1.2\n')
fp.write('!date: %s $\n' % (mgi_utils.date("%m/%d/%Y")))
fp.write('!\n')
fp.write('! from Mouse Genome Database (MGD) & Gene Expression Database (GXD)\n')
fp.write('!\n')
fp.write('!  DB                     required  1             1             MGI\n')
fp.write('!  DB_Object_ID           required  1             2/17          MGI:87870\n')
fp.write('!  DB_Object_Symbol       required  1             3             Acat1\n')
fp.write('!  DB_Object_Name         optional  0 or greater  10            acetyl-Coenzyme A acetyltransferase 1\n')
fp.write('!  DB_Object_Synonym(s)   optional  0 or greater  11            Acat|6330585C21Rik\n')
fp.write('!  DB_Object_Type         required  1             12            gene\n')
fp.write('!  Taxon                  required  1             13            taxon:10090\n')
fp.write('!  Parent_Object_ID       optional  0 or 1        -             if DB_Object_ID = Isoform, then MGI:id of the Isoform\n')
fp.write('!  DB_Xref(s)             optional  0 or greater  -             if DB_Object_ID = Isoform, then UniProtKB:id of the Isoform\n')
fp.write('!  Properties             optional  0 or greater  -             blank\n\n')
fp.write('!  DB_Object_Type = "gene", DB = "MGI", DB_Object_ID = "MGI:xxxx"\n')
fp.write('!  DB_Object_Type = "protein", DB = "PR", DB_Object_ID = "xxxx", Parent_Object_ID = "MGI:MGI:xxxx", DB_Xref = "UniProtDB:xxx"\n')
fp.write('!  DB_Object_Type = "transcript", DB = "EMBL", "RefSeq", "ENSEMBL", "VEGA", Parent_Object_ID = "MGI:MGI:xxxx"\n')
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
	select a1.accID as markerID, a2.accID as prID
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
	select a.accID as accID, s.synonym
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
# markerUniProtKB primary
#
# hard-coded (for now)...
# read: ${DATADOWNLOADS}/ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/goa_mouse.gpi.gz
# field 1 = UniProtKB 
# field 2 = xxxxx
# find marker _object_key where accID = 'xxxxx'
#	and _logicaldb_key in (13, 14) and _mgitype_key = 2
# add UniProtDB:xxxx to gpi field 8
#

uniprotGPI = []
gpiFile = gzip.open(os.environ['DATADOWNLOADS'] + '/ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/goa_mouse.gpi.gz', 'r')
for line in gpiFile.readlines():
        if line[0] == '!':
            continue
        tokens = line[:-1].split('\t')
	uniprotGPI.append(tokens[1])
gpiFile.close()

markerUniProtKB = {}
results = db.sql('''
        select sm.accID as uniprotID, a.accID as markerID, sm._LogicalDB_key
        from SEQ_Marker_Cache sm, ACC_Accession a
        where sm._LogicalDB_key in (13,41)
        and sm._Marker_Type_key = 1 
        and sm._Organism_key = 1 
        and sm._Marker_key = a._Object_key
        and a._MGIType_key = 2 
        and a._LogicalDB_key = 1 
        and a.preferred = 1 
        order by sm._LogicalDB_key
	''', 'auto')
for r in results:
	key = r['markerID']
	value = r['uniprotID']
	if value in uniprotGPI:
	    if key not in markerUniProtKB:
		    markerUniProtKB[key] = []
	    markerUniProtKB[key].append('UniProtKB:' + value)

#
# markers:terms
#

results = db.sql('''
	select a.accID as accID, m.symbol, m.name
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

	fp.write(DB_PREFIX + TAB)
	fp.write(r['accID'] + TAB)
	fp.write(r['symbol'] + TAB)
	fp.write(r['name'] + TAB)

	if marker in markerSynonyms:
		fp.write("|".join(markerSynonyms[marker]))
	fp.write(TAB)

	fp.write(DBTYPE_MARKER + TAB)
	fp.write(SPECIES + TAB)

	if marker in markerUniProtKB:
		fp.write("|".join(markerUniProtKB[marker]))
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
	prefixPart, accID = r['accID'].split(':')
	isoform = r['accID']

	fp.write(prefixPart + TAB)
	fp.write(accID + TAB)
	fp.write(r['symbol'] + TAB)
	fp.write(r['name'] + TAB)

	if isoform in isoformSynonyms:
		fp.write("|".join(isoformSynonyms[isoform]))
	fp.write(TAB)

	fp.write(DBTYPE_ISOFORM + TAB)
	fp.write(SPECIES + TAB)

	# Parent_Object_ID/MGI id
	if isoform in markerByIsoform:
		fp.write('MGI:' + markerByIsoform[isoform][0])
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
	select distinct sm._Sequence_key, sm.accID as rnaID, sm._Marker_key, sm._LogicalDB_key
        	from SEQ_Marker_Cache sm, ACC_Accession a
        	where sm._SequenceType_key = 316346
        	and sm._Marker_Type_key = 1 
        	and sm._Organism_key = 1 
		and sm._LogicalDB_key in (9,27,131,133)
		and sm._Sequence_key = a._Object_key
		and a._MGIType_key = 19
	        and a.prefixPart not in ('XR_', 'XM_')
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
		ldbName = 'EMBL'
	elif ldb == 27:
		ldbName = 'RefSeq'
	elif ldb == 133:
		ldbName = 'ENSEMBL'
	elif ldb == 131:
		ldbName = 'VEGA'

	fp.write(ldbName + TAB)
	fp.write(r['rnaID'] + TAB)
	fp.write(r['rnaID'] + TAB)
	fp.write(TAB)
	fp.write(TAB)
	fp.write(DBTYPE_RNA + TAB)
	fp.write(SPECIES + TAB)
	fp.write('MGI:' + r['markerID'] + TAB)
	fp.write(TAB)
	fp.write(CRT)

# end RNAs

reportlib.finish_nonps(fp)
db.useOneConnection(0)

