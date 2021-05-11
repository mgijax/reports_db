'''
#
# GO_gpi.py
#
# cat mgi.gpi | grep -v '^!' | awk 'BEGIN {FS="\t"} {print NF}' | sort | uniq -c
#
# Report:
#       see: wiki/mediawiki/index.php/sw:GOload 
#       contains link to GO/GPI format
#
# 1  DB_Object_ID
# 2  DB_Object_Symbol
# 3  DB_Object_Name
# 4  DB_Object_Synonyms 
# 5  DB_Object_Type
# 6  DB_Object_Taxon
# 7  Encoded_By 
# 8  Parent_Protein 
# 9  Protein_Containing_Complex_Members 
# 10 DB_Xrefs 
# 11 Gene_Product_Properties 
#
# History:
#
# lec   04/01/2020 python 3 upgrade
#       - TR13272/GPI 2.0
#
# lec   04/01/2020 python 3 upgrade
#
# sc    03/21/2020 python 3 upgrade
#
# 06/27/2016    lec
#       - TR12349/12345/GPAD/GPI
#
'''

import sys
import os
import gzip
import mgi_utils
import reportlib
import db

db.setTrace()

SPECIES = 'NCBITaxon:10090'

TAB = reportlib.TAB
CRT = reportlib.CRT

db.useOneConnection(1)

fp = reportlib.init('mgi', fileExt = '.gpi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('!gpi-version: 2.0\n')
fp.write('!namespace: MGI\n')
fp.write('!generated-by: MGI\n')
fp.write('!date-generated: %s\n' % (mgi_utils.date("%Y-%m-%d")))

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
        and exists (select 1 from MRK_Marker m 
                where v._Object_key = m._Marker_key
                and m._Marker_Status_key = 1
                )
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
        and m._Marker_Status_key = 1
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
#       and _logicaldb_key in (13, 14) and _mgitype_key = 2
# add UniProtDB:xxxx to gpi field 8
#

uniprotGPI = []
gpiFile = gzip.open(os.environ['DATADOWNLOADS'] + '/ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/goa_mouse.gpi.gz', 'rt')
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
        and sm._Organism_key = 1 
        and sm._Marker_key = a._Object_key
        and a._MGIType_key = 2 
        and a._LogicalDB_key = 1 
        and a.preferred = 1 
        and exists (select 1 from MRK_Marker m 
                where a._Object_key = m._Marker_key
                and m._Marker_Status_key = 1
                )
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
# marker with feature type
#

results = db.sql('''
        select a.accID as accID, m.symbol, m.name, m._Marker_Type_key, t.term
        from ACC_Accession a, MRK_Marker m, VOC_Annot v, VOC_Term t
        where a._MGIType_key = 2
        and a._LogicalDB_key = 1
        and a.prefixPart = 'MGI:'
        and a.preferred = 1
        and a._Object_key = m._Marker_key
        and m._Marker_Status_key = 1
        and m._Marker_Type_key in (1,7,10)
        and m._Organism_key = 1
        and m._Marker_key = v._Object_key
        and v._AnnotType_key = 1011
        and v._Term_key = t._Term_key
        order by m.symbol
        ''', 'auto')

for r in results:
        marker = r['accID']

        fp.write('MGI:' + r['accID'] + TAB)
        fp.write(r['symbol'] + TAB)
        fp.write(r['name'] + TAB)

        if marker in markerSynonyms:
                fp.write("|".join(markerSynonyms[marker]))
        fp.write(TAB)

        # 5. DB_Object_Type
        if r['_Marker_Type_key'] == 1 and r['term'] == 'unclassified gene':
                fp.write('SO:0000704' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'protein coding gene':
                fp.write('SO:0001217' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'non-coding RNA gene':
                fp.write('SO:0001263' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'heritable phenotypic marker':
                fp.write('SO:0001645' + TAB)
        elif r['_Marker_Type_key'] == 1 and \
                        (r['term'] == 'lncRNA gene' \
                        or r['term'] == 'antisense lncRNA gene' \
                        or r['term'] == 'lincRNA gene' \
                        or r['term'] == 'sense intronic lncRNA gene' \
                        or r['term'] == 'sense overlapping lncRNA gene' \
                        or r['term'] == 'bidirectional promoter lncRNA gene' \
                        or r['term'] == 'rRNA gene' \
                        or r['term'] == 'tRNA gene' \
                        or r['term'] == 'snRNA gene' \
                        or r['term'] == 'snoRNA gene' \
                        or r['term'] == 'miRNA gene' \
                        or r['term'] == 'scRNA gene' \
                        or r['term'] == 'SRP RNA gene' \
                        or r['term'] == 'RNase P RNA gene' \
                        or r['term'] == 'RNase MRP RNA gene' \
                        or r['term'] == 'telomerase RNA gene' \
                        or r['term'] == 'unclassified non-coding RNA gene' \
                        or r['term'] == 'ribozyme gene'):
                fp.write('SO:0001263' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'gene segment':
                fp.write('SO:3000000' + TAB)
        elif r['_Marker_Type_key'] == 7:
                fp.write('SO:0000336' + TAB)
        elif r['_Marker_Type_key'] == 10:
                fp.write('SO:0001411' + TAB)
        else:
                fp.write('col 5 bad with feature type =  ' + r['term'] + TAB)

        fp.write(SPECIES + TAB)
        fp.write(TAB)
        fp.write(TAB)
        fp.write(TAB)

        if marker in markerUniProtKB:
                fp.write("|".join(markerUniProtKB[marker]))
        fp.write(TAB)

        fp.write(CRT)

#
# marker without feature type; use marker type
#

results = db.sql('''
	select a.accID as accID, m.symbol, m.name, m._Marker_Type_key
	from ACC_Accession a, MRK_Marker m
	where a._MGIType_key = 2
	and a._LogicalDB_key = 1
	and a.prefixPart = 'MGI:'
	and a.preferred = 1
	and a._Object_key = m._Marker_key
	and m._Marker_Status_key = 1
        and m._Marker_Type_key in (1,7,10)
	and m._Organism_key = 1
        and not exists (select v.* from VOC_Annot v
                where m._Marker_key = v._Object_key
                and v._AnnotType_key = 1011
                )
        order by m.symbol
   	''', 'auto')

for r in results:
	marker = r['accID']

	fp.write('MGI:' + r['accID'] + TAB)
	fp.write(r['symbol'] + TAB)
	fp.write(r['name'] + TAB)

	if marker in markerSynonyms:
		fp.write("|".join(markerSynonyms[marker]))
	fp.write(TAB)

	if r['_Marker_Type_key'] == 7:
                fp.write('SO:0000336' + TAB)
	elif r['_Marker_Type_key'] == 10:
                fp.write('SO:0001411' + TAB)
	else:
                fp.write('col 5 bad no feature type/marker type =  ' + str(r['_Marker_Type_key']) + TAB)

	fp.write(SPECIES + TAB)
	fp.write(TAB)
	fp.write(TAB)
	fp.write(TAB)

	if marker in markerUniProtKB:
		fp.write("|".join(markerUniProtKB[marker]))
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
        and exists (select 1 from MRK_Marker m 
                where v._Object_key = m._Marker_key
                and m._Marker_Status_key = 1
                )
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
        select a.accID, t.term as symbol, t.note as name
        from ACC_Accession a, VOC_Term t
        where t._Vocab_key = 112
        and t._term_key = a._Object_key 
        and a._LogicalDB_key = 183
        order by symbol, a.accID
        ''', 'auto')

for r in results:
        isoform = r['accID']

        if isoform not in markerByIsoform:
            continue

        fp.write(isoform + TAB)
        fp.write(r['symbol'] + TAB)
        fp.write(r['name'] + TAB)

        if isoform in isoformSynonyms:
                fp.write("|".join(isoformSynonyms[isoform]))
        fp.write(TAB)

        fp.write('PR:000000001' + TAB)
        fp.write(SPECIES + TAB)

        # 7. Encoded_By                                 ::= [ID] ('|' ID)*
        fp.write('MGI:' + markerByIsoform[isoform][0] + TAB)

        # 8 Parent_Protein                             ::= [ID] ('|' ID)*
        fp.write(TAB)

        # 9 Protein_Containing_Complex_Members         ::= [ID] ('|' ID)*
        fp.write(TAB)

        # 10 DB_Xrefs                                  ::= [ID] ('|' ID)*
        if isoform in isoformUniProt:
                fp.write(isoformUniProt[isoform][0])
        fp.write(TAB)

        # 11 Gene_Product_Properties  
        fp.write(CRT)

# attach PR/protein_complex generated from proisoformload 
#try:
gpiFileName = os.environ['DATALOADSOUTPUT'] + '/pro/proisoformload/output/gpi2.txt'
gpiFile = open(gpiFileName, 'r')
for line in gpiFile.readlines():
        fp.write(line)
#except:
#        exit(1, 'Could not open file %s\n' % gpiFileName)

# end isoforms
#

#
# RNAs that are associated with at most 1 marker
#

for ldbsearch in (9, 27, 133):

    results = db.sql('''
	    WITH rna AS
	    (
	    select distinct sm._Sequence_key, sm.accID as rnaID, sm._Marker_key, sm._LogicalDB_key, t.term
        	    from SEQ_Marker_Cache sm, ACC_Accession a, VOC_Annot v, VOC_Term t
        	    where sm._SequenceType_key = 316346
        	    and sm._Marker_Type_key = 1 
        	    and sm._Organism_key = 1 
		    and sm._LogicalDB_key = %s
		    and sm._Sequence_key = a._Object_key
		    and a._MGIType_key = 19
	            and a.prefixPart not in ('XR_', 'XM_')
                    and sm._Marker_key = v._Object_key
                    and v._AnnotType_key = 1011
                    and v._Term_key = t._Term_key
	        )
	    select rna.rnaID, a.accID as markerID, rna._LogicalDB_key, rna.term
	    from rna, ACC_Accession a
	    where rnaID in (select rnaID from rna group by rnaID having count(*) = 1)
            and rna._Marker_key = a._Object_key
            and a._MGIType_key = 2
            and a._LogicalDB_key = 1
            and a.preferred = 1
	    order by rna._LogicalDB_key
   	    ''' % (ldbsearch), 'auto')

    for r in results:
    
	    ldb = r['_LogicalDB_key']
	    if ldb == 9:
		    ldbName = 'EMBL'
	    elif ldb == 27:
		    ldbName = 'RefSeq'
	    elif ldb == 133:
		    ldbName = 'ENSEMBL'
    
	    # if protein coding gene
	    if r['term'] == 'protein coding gene':
                DBTYPE_RNA = 'SO:0000234'
	    else:
                DBTYPE_RNA = 'SO:0000655'
                
	    fp.write(ldbName + ':' + r['rnaID'] + TAB)
	    fp.write(r['rnaID'] + TAB)
	    fp.write(TAB)
	    fp.write(TAB)
	    fp.write(DBTYPE_RNA + TAB)
	    fp.write(SPECIES + TAB)
	    fp.write('MGI:' + r['markerID'] + TAB)
	    fp.write(TAB)
	    fp.write(TAB)
	    fp.write(TAB)
	    fp.write(CRT)

# end RNAs

reportlib.finish_nonps(fp)
db.useOneConnection(0)
