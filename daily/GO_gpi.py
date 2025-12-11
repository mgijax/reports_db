'''
#
# GO_gpi.py
#
# cat mgi.gpi | grep -v '^!' | awk 'BEGIN {FS="\t"} {print NF}' | sort | uniq -c
#
# 4 marker sets:
#   A: markers with feature types; contains col10 = RNAcentral
#   B: markers without feature types; contains col10 = RNAcentral
#   C: isoforms PR (183) : 1,7,10 only
#   D: RNAs that are associated with at most 1 marker : 1,7,10 only; contains col10 = RNAcentral
#
# Report:
#       see: wiki/mediawiki/index.php/sw:GOload 
#       contains link to GO/GPI format
#
# 1  DB:Object ID
# 2  Object Symbol
# 3  Object Name
# 4  Object Synonyms 
# 5  Object Type
# 6  Object Taxon
#
# 7  Encoded by
# 8  Canonical object ID
# if column 5 (Object Type) in (PR:000000001, SO:0000234, SO:0000655)
#   column 7, column 8 = MGI:xxx
# else:
#   column 8 = column 1
# 
# 9  Protein Complex Members
# 10 Cross-reference(s) DB_Xrefs 
# include the corresponding RNA_Central ID for each MGI feature
#
# 11 Gene Product Properties 
#
# History:
#
# lec   12/11/2025 python 3 upgrade
#       - wts2-1771/sprt-153/MGI GPI missing GCRPs as cross reference
#
'''

import sys
import os
import gzip
import mgi_utils
import reportlib
import db

#db.setTrace()

SPECIES = 'NCBITaxon:10090'

TAB = reportlib.TAB
CRT = reportlib.CRT

rnaTag = 'RNAcentral:%s_10090'
singleMgiToRnaCentral = {}
singleRnaCentralToMGI = {}

def rnaCentral():
    global singleMgiToRnaCentral
    global singleRnaCentralToMGI

    #
    # create set of single MGI -> RNACentral associations
    # create set of single RNACentral -> MGI associations
    # The single RNACentral in the single MGI -> RNACentral set
    # must also exist in the single RNACentral -> MGI set
    #

    allMgiToRnaCentral = {}
    allRnaCentralToMGI = {}

    inFile = open(os.environ['DATADOWNLOADS'] + '/ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/mgi.tsv', 'r')
    for line in inFile.readlines():
        tokens = line[:-1].split('\t')
        rnaId = tokens[0]
        mgiId = tokens[2]

        if mgiId not in allMgiToRnaCentral:
            allMgiToRnaCentral[mgiId] = []
        allMgiToRnaCentral[mgiId].append(rnaId)

        if rnaId not in allRnaCentralToMGI:
            allRnaCentralToMGI[rnaId] = []
        allRnaCentralToMGI[rnaId].append(mgiId)
    inFile.close()

    # create set of single MGI -> RNACentral associations
    for rnaId in allRnaCentralToMGI:
        if len(allRnaCentralToMGI[rnaId]) == 1:
            singleRnaCentralToMGI[rnaId] = allRnaCentralToMGI[rnaId]

    # create set of single RNACentral -> MGI associations
    # The single RNACentral in the single MGI -> RNACentral set
    # must also exist in the single RNACentral -> MGI set
    for mgiId in allMgiToRnaCentral:
        if len(allMgiToRnaCentral[mgiId]) == 1:
            rnaId = allMgiToRnaCentral[mgiId][0]
            if rnaId in singleRnaCentralToMGI:
                singleMgiToRnaCentral[mgiId] = allMgiToRnaCentral[mgiId]

    #print(allRnaCentralToMGI)
    #print(allMgiToRnaCentral)
    #print(singleMgiToRnaCentral)
    #print(singleRnaCentralToMGI)

#
# main
#

rnaCentral()

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
# hard-coded (for now)...
# read: ${DATADOWNLOADS}/ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/goa_mouse.gpi.gz
# field 1 = UniProtKB 
# add UniProtDB:xxxx to gpi field 8
#

uniprotGPI = {}
gpiFile = gzip.open(os.environ['DATADOWNLOADS'] + '/ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/goa_mouse.gpi.gz', 'rt')
for line in gpiFile.readlines():

        if line[0] == '!':
            continue

        tokens = line[:-1].split('\t')
        id = tokens[1]
        fullid = tokens[0] + ':' + tokens[1]

        # search uniprot id -> marker relationships (uniprotload)
        results = db.sql('''
            select m.symbol
            from acc_accession a1, mrk_marker m
            where a1.accid = '%s'
            and a1._mgitype_key = 2 
            and a1._logicaldb_key in (13,41)
            and a1._object_key = m._marker_key
            ''' % id, 'auto')

        # only interested in 1:1 uniprotload relationships
        if len(results) == 1:
            for r in results:
                symbol = r['symbol']
                if symbol not in uniprotGPI:
                    uniprotGPI[symbol] = []
                uniprotGPI[symbol].append(fullid)

        # if 1:N, if db_subset=Swiss-Prot, then use gpiFile/symbol
        elif len(results) > 1:
            dbSubSet = tokens[9]
            if dbSubSet == 'db_subset=Swiss-Prot':
                symbol = tokens[2]
                if symbol not in uniprotGPI:
                    uniprotGPI[symbol] = []
                uniprotGPI[symbol].append(fullid)
            else:
                print(id, results)
        else:
            print(id, results)

gpiFile.close()

#
# markers:terms
# A: marker with feature type
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
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'antisense lncRNA gene':      
                fp.write('SO:0002182' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'bidirectional promoter lncRNA gene':  
                fp.write('SO:0002185' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'lincRNA gene':        
                fp.write('SO:0001641' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'lncRNA gene': 
                fp.write('SO:0002127' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'miRNA gene':  
                fp.write('SO:0001265' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'non-coding RNA gene': 
                fp.write('SO:0001263' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'ribozyme gene':       
                fp.write('SO:0002181' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'RNase MRP RNA gene':  
                fp.write('SO:0001640' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'RNase P RNA gene':    
                fp.write('SO:0001639' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'rRNA gene':   
                fp.write('SO:0001637' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'scRNA gene':  
                fp.write('SO:0001266' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'sense intronic lncRNA gene':  
                fp.write('SO:0002184' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'sense overlapping lncRNA gene':       
                fp.write('SO:0002183' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'snoRNA gene': 
                fp.write('SO:0001267' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'snRNA gene':  
                fp.write('SO:0001268' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'SRP RNA gene':        
                fp.write('SO:0001269' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'telomerase RNA gene': 
                fp.write('SO:0001643' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'tRNA gene':   
                fp.write('SO:0001272' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'unclassified non-coding RNA gene':    
                fp.write('SO:0001263' + TAB)
        elif r['_Marker_Type_key'] == 1 and r['term'] == 'gene segment':
                fp.write('SO:3000000' + TAB)
        elif r['_Marker_Type_key'] == 7:
                fp.write('SO:0000336' + TAB)
        elif r['_Marker_Type_key'] == 10:
                fp.write('SO:0001411' + TAB)
        else:
                fp.write('col 5 bad with feature type = ' + r['term'] + TAB)

        fp.write(SPECIES + TAB)
        fp.write(TAB)
        fp.write('MGI:' + r['accID'] + TAB)
        fp.write(TAB)

        addPipe = ""
        symbol = r['symbol']
        if symbol in uniprotGPI:
                fp.write("|".join(uniprotGPI[symbol]))
                addPipe = "|"
        if marker in singleMgiToRnaCentral:
                fp.write(addPipe + rnaTag % (singleMgiToRnaCentral[marker][0]))
        fp.write(TAB)

        fp.write(CRT)

#
# B: marker without feature type; use marker type
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
	fp.write('MGI:' + r['accID'] + TAB) 
	fp.write(TAB)

	addPipe = ""
	symbol = r['symbol']
	if symbol in uniprotGPI:
		fp.write("|".join(uniprotGPI[symbol]))
		addPipe = "|"
	if marker in singleMgiToRnaCentral:
		fp.write(addPipe + rnaTag % (singleMgiToRnaCentral[marker][0]))
	fp.write(TAB)

	fp.write(CRT)

#
# end markers
#

#
# isoforms:
# 183 | Protein Isoform Ontology
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
        and exists (select 1 from VOC_Annot v where v._annottype_key = 1019 and t._term_key = v._term_key)
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
                and m._Marker_Type_key in (1,7,10)
                )
        ''', 'auto')
for r in results:
        key = r['prID']
        value = r['value']
        if key not in isoformUniProt:
                isoformUniProt[key] = []
        isoformUniProt[key].append(value)

#
# C: isoforms:terms
# 183 | Protein Isoform Ontology
#

results = db.sql('''
        select a.accID, t.term as symbol, t.note as name
        from ACC_Accession a, VOC_Term t
        where t._Vocab_key = 112
        and t._term_key = a._Object_key 
        and a._LogicalDB_key = 183
        and exists (select 1 from VOC_Annot v, MRK_Marker m 
                where v._annottype_key = 1019 
                and t._term_key = v._term_key
                and v._Object_key = m._Marker_key
                and m._Marker_Status_key = 1
                and m._Marker_Type_key in (1,7,10)
                )
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

        # 8  Canonical object ID
        fp.write('MGI:' + markerByIsoform[isoform][0] + TAB)

        # 9 Protein_Containing_Complex_Members         ::= [ID] ('|' ID)*
        fp.write(TAB)

        # 10 DB_Xrefs                                  ::= [ID] ('|' ID)*
        if isoform in isoformUniProt:
                fp.write(isoformUniProt[isoform][0])
        fp.write(TAB)

        # 11 Gene_Product_Properties  
        fp.write(CRT)

reportlib.finish_nonps(fp)
