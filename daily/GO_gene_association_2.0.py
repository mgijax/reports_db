'''
#
# GO_gene_association_2.0.py
#
# Generates:
#       
#	mgi2.gpad
#	gene_association.mgi2 (GAF)
#	gene_association_pro.mgi2 (GAF)
#
#
# GPAD
#
#!1. DB_Object_ID ::= ID       MGI or PR
#!2. Negation ::= 'NOT'
#!3. Relation ::= OBO_ID
#!4. Ontology_Class_ID ::= OBO_ID
#!5. Reference ::= [ID] ('|' ID)*
#!6. Evidence_type ::= OBO_ID
#!7. With_or_From ::= [ID] ('|' | ‘,’ ID)*
#!8. Interacting_taxon_ID ::= NCBITaxon:[Taxon_ID]
#!9. Date ::= YYYY-MM-DD
#!10. Assigned_by ::= Prefix
#!11. Annotation_Extensions ::= [Extension_Conj] ('|' Extension_Conj)*
#!12. Annotation_Properties ::= [Property_Value_Pair] ('|' Property_Value_Pair)*
#
# GAF
#
#!1.  DB                              required        1               UniProtKB\n')
#!2.  DB Object ID                    required        1               P12345\n')
#!3.  DB Object Symbol                required        1               PHO3\n')
#!4.  Qualifier                       optional        0 or greater    NOT\n')
#!5.  GO ID                           required        1               GO:0003993\n')
#!6.  DB:Reference (|DB:Reference)    required        1 or greater    PMID:2676709\n')
#!7.  Evidence Code                   required        1               IMP\n')
#!8.  With (or) From                  optional        0 or greater    GO:0000346\n')
#!9.  Aspect                          required        1               F\n')
#!10. DB Object Name                  optional        0 or 1          Toll-like receptor 4\n')
#!11. DB Object Synonym (|Synonym)    optional        0 or greater    hToll   Tollboot\n')
#!12. DB Object Type                  required        1               protein\n')
#!13. Taxon(|taxon)                   required        1 or 2          taxon:9606\n')
#!14. Date                            required        1               20090118\n')
#!15. Assigned By                     required        1               SGD\n')
#!16. Annotation Extension            optional        0 or greater    part_of(CL:0000576)\n')
#!17. Gene Product Form ID            optional        0 or 1          UniProtKB:P12345-2\n')
# 
# lec   08/25/2020
#       - TR13272/converting to GPI 2.0
#       mgi2.gpad : dph reviewing 09/09/2020
#       gene_association.mgi2 : no changes yet
#       gene_association_pro.mgi2 : no changes yet
#
'''

import sys
import os
import mgi_utils
import reportlib
import go_annot_extensions
import go_isoforms
import db

goloadpath = os.environ['GOLOAD'] + '/lib'
sys.path.insert(0, goloadpath)
import ecolib

db.setTrace()

MGIPREFIX = 'MGI'
SPECIES = 'NCBITaxon:10090'

#
# if in list 1, then use 'UniProt'
# if in list 2, then use itself
# else it will use MGIPREFIX
#
assignedByList1 = ['uniprotload']
assignedByList2 = ['GOC', 'GO_Central']

TAB = reportlib.TAB
CRT = reportlib.CRT

# see doSetup()
dag = {}
syns = {}
pubMed = {}

# mapping of load reference keys to GO_REF IDs/per TR11772
# see _LogicalDB_key = 185 (GO_REF)
goRefDict = {}

# gaf lookups
#
# see doGAFCol16() : list of column 16 object/evidence/printable format
gafCol16Lookup = {}

# see doIsoform() : isoformsProtein = {}
forPROC = {}

# see doProtein() : protein identifiers
proteins = {}

#
# gpad lookups
#
# translate dag to qualifier (col 3)
dagQualifier = {'C':'BFO:0000050', 'P':'RO:0002264', 'F':'RO:0002327'}
taxonLookup = {}
ecoLookupByEco = {}
ecoLookupByEvidence = {}
evidenceLookup = {}
goPropertyLookup = {}
goQualifierLookup = {}
gpadCol11Lookup = {}
gpadCol12Lookup = {}

#
# begin doSetup()
#
def doSetup():
    global dag
    global syns
    global pubMed
    global evidenceLookup
    global taxonLookup
    global ecoLookupByEco, ecoLookupByEvidence
    global goPropertyLookup
    global goQualifierLookup
    global gpadCol11Lookup
    global gpadCol12Lookup
    global goRefDict

    #
    # retrieve all dag abbrevations for each term
    #
    results = db.sql('select distinct _Object_key, rtrim(dagAbbrev) as dagAbbrev from DAG_Node_View where _Vocab_key = 4', 'auto')
    for r in results:
        dag[r['_Object_key']] = r['dagAbbrev']

    #
    # retrieve data set to process
    #
    #   and m.symbol = 'Mbd2'
    #   and m.symbol = 'Adipoq'
    #   and m.symbol = 'Birc3'
    #	and m.symbol = 'Hk1'
    #
    db.sql('''select distinct a._Term_key, t.term, ta.accID as termID, q.term as qualifier, a._Object_key, 
            e._AnnotEvidence_key, e.inferredFrom, e._EvidenceTerm_key, 
            e._Refs_key, e._CreatedBy_key, e._ModifiedBy_key, e.creation_date, e.modification_date,
            m.symbol, m.name, lower(mt.name) as markerType
        into temporary table gomarker1 
        from VOC_Annot a, 
             ACC_Accession ta, 
             VOC_Term t, 
             VOC_Evidence e, 
             MRK_Marker m, 
             MRK_Types mt, 
             VOC_Term q
        where a._AnnotType_key = 1000 
        and a._Annot_key = e._Annot_key 
        and a._Object_key = m._Marker_key 
        and m._Marker_Type_key = 1 
        and m._Marker_Status_key = 1
        and a._Term_key = t._Term_key 
        and a._Term_key = ta._Object_key 
        and ta._MGIType_key = 13 
        and ta.preferred = 1 
        and m._Marker_Type_key = mt._Marker_Type_key 
        and a._Qualifier_key = q._Term_key 
        ''', None)
    db.sql('create index gomarker1_idx1 on gomarker1(_Object_key)', None)
    db.sql('create index gomarker1_idx2 on gomarker1(_EvidenceTerm_key)', None)
    db.sql('create index gomarker1_idx3 on gomarker1(_Refs_key)', None)
    db.sql('create index gomarker1_idx4 on gomarker1(_ModifiedBy_key)', None)
    db.sql('create index gomarker1_idx5 on gomarker1(_AnnotEvidence_key)', None)

    #
    # retrieve synonyms for markers in data set
    #
    results = db.sql('''select distinct g._Object_key, s.synonym 
        from gomarker1 g, MGI_Synonym s, MGI_SynonymType st 
        where g._Object_key = s._Object_key 
        and s._MGIType_key = 2 
        and s._SynonymType_key = st._SynonymType_key 
        and st.synonymType = 'exact'
        order by g._Object_key''', 'auto')
    for r in results:
        key = r['_Object_key']
        value = r['synonym']
        if key not in syns:
            syns[key] = []
        syns[key].append(value)

    #
    # resolve foreign keys and store in "results" table
    #
    db.sql('''select distinct g._Refs_key, g._Term_key, g.termID, g.qualifier, g.inferredFrom, 
            g._Object_key, g._AnnotEvidence_key, g._EvidenceTerm_key, g.symbol, g.name, g.markerType, 
            to_char(g.creation_date, 'YYYY-MM-DD') as cDate,
            to_char(g.modification_date, 'YYYY-MM-DD') as mDate,
            g._CreatedBy_key,
            g._ModifiedBy_key,
            ma.accID as markerID, 
            b.accID as refID, 
            rtrim(t.abbreviation) as evidenceCode, 
            u.login as assignedBy
        into temporary table gomarker2 
        from gomarker1 g, ACC_Accession ma, ACC_Accession b, VOC_Term t, MGI_User u 
        where g._Object_key = ma._Object_key 
        and ma._MGIType_key = 2 
        and ma.prefixPart = 'MGI:' 
        and ma._LogicalDB_key = 1 
        and ma.preferred = 1 
        and g._Refs_key = b._Object_key 
        and b._MGIType_key = 1 
        and b.prefixPart = 'MGI:' 
        and b._LogicalDB_key = 1 
        and g._EvidenceTerm_key = t._Term_key 
        and g._ModifiedBy_key = u._User_key
        ''', None)
    db.sql('create index gomarker2_idx1 on gomarker2(symbol)', None)
    db.sql('create index gomarker2_idx2 on gomarker2(_Refs_key)', None)
    db.sql('create index gomarker2_idx3 on gomarker2(_AnnotEvidence_key)', None)
    db.sql('create index gomarker2_idx4 on gomarker2(_Object_key)', None)
    db.sql('create index gomarker2_idx5 on gomarker2(_EvidenceTerm_key)', None)
    
    #
    # resolve PubMed IDs for References
    #
    results = db.sql('''select distinct r._Refs_key, a.accID from gomarker2 r, ACC_Accession a 
        where r._Refs_key = a._Object_key 
        and a._MGIType_key = 1 
        and a._LogicalDB_key = 29''', 'auto')
    for r in results:
        key = r['_Refs_key']
        value = r['accID']
        pubMed[key] = value

    #
    # ecoLookupByEco
    # ecoLookupByEvidence
    # evidenceLookup : evidence
    #
    # use goload/ecolib to lookup eco using evidencevidenceCode
    ecoLookupByEco, ecoLookupByEvidence = ecolib.processECO()
    #print(ecoLookupByEvidence)

    results = db.sql('''select distinct a._AnnotEvidence_key, t.term, p.value
            from gomarker2 a,
                 VOC_Evidence_Property p,  
                 VOC_Term t
            where a._AnnotEvidence_key = p._AnnotEvidence_key
            and p._PropertyTerm_key = t._Term_key
            and t.term = 'evidence'
            ''', 'auto')
    for r in results:
        key = r['_AnnotEvidence_key']
        value = ''
        tokens = r['value'].split(' ')
        for t in tokens:
            if t.find('ECO:') >= 0:
                value = t
        if len(value) > 0 and value in ecoLookupByEco:
            if key not in evidenceLookup:
                evidenceLookup[key] = []
            evidenceLookup[key].append(value)
    #print(evidenceLookup)

    #
    # taxonLookup : dual-taxon ID
    #
    results = db.sql('''select distinct a._AnnotEvidence_key, t.term, p.value
            from gomarker2 a,  
                 VOC_Evidence_Property p,  
                 VOC_Term t
            where a._AnnotEvidence_key = p._AnnotEvidence_key
            and p._PropertyTerm_key = t._Term_key
            and t.term = 'dual-taxon ID'
            ''', 'auto')
    for r in results:
        key = r['_AnnotEvidence_key']
        value = r['value']
        if key not in taxonLookup:
            taxonLookup[key] = []
        taxonLookup[key].append(value)
    #print(taxonLookup)

    #
    # GO/Property : note (definition) contains RO/.etc
    #
    results = db.sql('''select t.term, t.note
            from VOC_Term t
            where t._vocab_key = 82
            and t.note is not null
            ''', 'auto')
    for r in results:
        key = r['term']
        value = r['note']
        if key not in goPropertyLookup:
            goPropertyLookup[key] = []
        goPropertyLookup[key].append(value)
    #print(goPropertyLookup)

    #
    # gpadCol3 : go_qualifier
    # go_qualifier -> value -> RO id
    #
    results = db.sql('''select distinct a._AnnotEvidence_key, t.term, p.value, t2.note
            from gomarker2 a,
                 VOC_Evidence_Property p,  
                 VOC_Term t,
                 VOC_Term t2
            where a._AnnotEvidence_key = p._AnnotEvidence_key
            and p._PropertyTerm_key = t._Term_key
            and t.term in ('go_qualifier')
            and p.value = t2.term
            and t2.note is not null
            ''', 'auto')
    for r in results:
        key = r['_AnnotEvidence_key']
        value = r['note']
        if key not in goQualifierLookup:
            goQualifierLookup[key] = []
        goQualifierLookup[key].append(value)
    #print(goQualifierLookup)

    #
    # gpadCol11 : convert properties to RO id
    # note that noctua-generated properties will *always* have one stanza
    #
    results = db.sql('''select distinct a._AnnotEvidence_key, t.note, p.value
            from gomarker2 a,
                 VOC_Evidence_Property p,  
                 VOC_Term t
            where a._AnnotEvidence_key = p._AnnotEvidence_key
            and p._PropertyTerm_key = t._Term_key
            and t.term not in (
                'evidence', 'anatomy', 'cell type', 'gene product', 'modification', 'target', 
                'external ref', 'text', 'dual-taxon ID',
                'noctua-model-id', 'contributor', 'individual', 'go_qualifier', 'model-state',
                'has_participant', 'regulates_o_has_participant'
                )
            and t.note is not null
            ''', 'auto')
    for r in results:
        key = r['_AnnotEvidence_key']

        value = r['value'].replace('MGI:', 'MGI:MGI:')

        # if value = "EMAPA:xxx TSxxx", then only use "EMAPA:xxx"
        if value.find('EMAPA') >= 0:
                 tokens = value.split(' ')
                 value = tokens[0]

        value = r['note'] + '(' + value + ')'

        if key not in gpadCol11Lookup:
            gpadCol11Lookup[key] = []
        gpadCol11Lookup[key].append(value)
    #print(gpadCol11Lookup)

    #
    # gpadCol12 
    # all annotations rows have creation/modification date
    # creation/modification dates
    #
    results = db.sql('''select distinct a._AnnotEvidence_key, a.cDate, a.mDate from gomarker2 a''', 'auto')
    for r in results:
        key = r['_AnnotEvidence_key']
        value = 'creation-date=' + r['cDate'] + '|' + 'modification-date=' + r['mDate']
        if key not in gpadCol12Lookup:
            gpadCol12Lookup[key] = []
        gpadCol12Lookup[key].append(value)
    #print(gpadCol12Lookup)

    #
    # gpadCol12 : use actual property values
    #
    # and (u.login like 'NOCTUA_%' or (u.orcid is not null and p._propertyterm_key = 18583062))
    #
    results = db.sql('''select distinct a._AnnotEvidence_key, t.term, p.value
            from gomarker2 a,
                 VOC_Evidence_Property p,  
                 VOC_Term t
            where a._AnnotEvidence_key = p._AnnotEvidence_key
            and p._PropertyTerm_key = t._Term_key
            and t.term in ('noctua-model-id', 'model-state', 
                'text', 'has_participant', 'regulates_o_has_participant'
                )
            ''', 'auto')
    for r in results:
        key = r['_AnnotEvidence_key']
        term = r['term']

        if term in ('noctua-model-id', 'model-state'):
                value = r['term'] + '=' + r['value']

        elif term in ('text'):
                value = 'comment=' + r['value']

        elif term in ('has_participant', 'regulates_o_has_participant'):
                value = 'comment=' + term + '(' + r['value'] + ')'

        if key not in gpadCol12Lookup:
            gpadCol12Lookup[key] = []
        gpadCol12Lookup[key].append(value)
    #print(gpadCol12Lookup)

    #
    # gpadCol12 
    # createdBy/orcid
    # modifiedBy/orcid
    #
    results = db.sql('''select distinct a._AnnotEvidence_key, u.orcid
            from gomarker2 a,
                 MGI_User u
            where a._CreatedBy_key = u._User_key
            and u.orcid is not null
            union
            select distinct a._AnnotEvidence_key, u.orcid
            from gomarker2 a,
                 MGI_User u
            where a._ModifiedBy_key = u._User_key
            and u.orcid is not null
            union
            select distinct a._AnnotEvidence_key, p.value
            from gomarker2 a,
                 VOC_Evidence_Property p,  
                 VOC_Term t
            where a._AnnotEvidence_key = p._AnnotEvidence_key
            and p._PropertyTerm_key = t._Term_key
            and t.term in ('contributor')
            ''', 'auto')
    for r in results:
        key = r['_AnnotEvidence_key']
        value = 'contributor=' + r['orcid']
        if key not in gpadCol12Lookup:
            gpadCol12Lookup[key] = []
        gpadCol12Lookup[key].append(value)
    #print(gpadCol12Lookup)

    results = db.sql('select _Object_key, accID from ACC_Accession where _LogicalDB_key = 185', 'auto')
    for r in results:
        key = r['_Object_key']
        value = r['accID']
        goRefDict[key] = value
    #print(goRefDict)

#
# end doSetup()
#

#
# begin doGAFCol16()
#
def doGAFCol16():

    global gafCol16Lookup

    #
    # list of object/evidence/stanza that are GO-properties
    #
    # objectKey = marker key:annotation/evidence key (must match objectKey in doGAFFinish())
    #
    # including properties that use go-sanctioned-property list
    # and we also include mgi-property ('cell type', 'anatomy', 'target')
    #
    # include: _Vocab_key = 82
    # include: go-sanctioned-properties (use _Term_key > 6481780 or sequenceNum > 9)
    # include: mgi-properties
    # exclude: annotations where evidence = ISO (3251466)
    #

    # Query the valid _term_keys for properties and evidence codes
    extensionProcessor = go_annot_extensions.Processor()
    sanctionedPropertyKeys = extensionProcessor.querySanctionedPropertyTermKeys()
    sanctionedEvidenceTermKeys = extensionProcessor.querySanctionedEvidenceTermKeys()

    propertyKeyClause = ",".join([str(k) for k in sanctionedPropertyKeys])
    evidenceKeyClause = ",".join([str(k) for k in sanctionedEvidenceTermKeys])

    cmd = '''
    select r.symbol, r._Object_key, r._AnnotEvidence_key, t.note as property, p.value, p.stanza
            from gomarker2 r, VOC_Evidence_Property p, VOC_Term t
            where r._EvidenceTerm_key in (%s)
            and r._AnnotEvidence_key = p._AnnotEvidence_key
            and p._PropertyTerm_key = t._Term_key
            and p._PropertyTerm_key in (%s)
            and t.note is not null
            order by r.symbol, 
                r._Object_key, 
                r._AnnotEvidence_key, 
                p.stanza, 
                p.sequenceNum, 
                property
    ''' % (evidenceKeyClause, propertyKeyClause)

    results = db.sql(cmd, 'auto')

    for r in results:
        objectKey = str(r['_Object_key']) + ':' + str(r['_AnnotEvidence_key'])
        value = r['value'].replace('MGI:', 'MGI:MGI:')

        # process out the comments, etc
        value = extensionProcessor.processValue(value)

        if objectKey not in gafCol16Lookup:
                gafCol16Lookup[objectKey] = []
                stanza = 0

        if stanza == 0:
                sep = ''
        elif r['stanza'] != stanza:
                sep = '|'
        else:
                sep = ','

        gafCol16Lookup[objectKey].append(sep + r['property'] + '(' + value + ')')
        stanza = r['stanza']

    #print9gafCol16Lookup)

## end doGAFCol16()

#
# begin doIsoform()
#
def doIsoform():
    global isoformsProtein
    global forPROC

    isoformsProtein = {}
    forPROC = {}

    # Query the valid _term_keys for properties
    isoformProcessor = go_isoforms.Processor()
    sanctionedPropertyKeys = isoformProcessor.querySanctionedPropertyTermKeys()

    propertyKeyClause = ",".join([str(k) for k in sanctionedPropertyKeys])

    results = db.sql('''select r._Object_key, r._AnnotEvidence_key, p.value
        from gomarker2 r, VOC_Evidence_Property p
        where r._AnnotEvidence_key = p._AnnotEvidence_key
        and p._PropertyTerm_key in (%s)
        order by r._Object_key, r._AnnotEvidence_key, p.stanza, p.sequenceNum
    ''' % propertyKeyClause, 'auto')

    for r in results:

        key = str(r['_Object_key']) + ':' + str(r['_AnnotEvidence_key'])
        value = r['value']

        isoformValues = isoformProcessor.processValue(value)

        for isoform in isoformValues:

            isoformsProtein.setdefault(key, []).append(isoform)

            #
            # TR11060
            # if UniProtKB:xxxx (any UniProtKB)
            # if PR:xxxx
            #
            if isoform.find('UniProtKB:') >= 0 or isoform.find('PR:') >= 0:
                forPROC.setdefault(key, []).append(isoform)

#
# end doIsoform()
#

#
# begin doProtein()
#
def doProtein():
    global proteins

    #
    # protein hash
    #
    # proteins hash
    #    representative proteins (615421) for marker type "gene" (1)
    #    13 = SwissProt
    #    41 = TrEMBL
    #    27 = RefSeq (XP, NP)
    #    134 = Ensembl Protein
    #
    # example of counts:
    #  3935 go-uniprot.txt   proteins
    #  4519 go-npxp.txt      proteins
    #

    results = db.sql('''
        select distinct r.symbol, mc._Marker_key, mc.accID as seqID, mc._LogicalDB_key, mc._Qualifier_key
        from gomarker2 r, SEQ_Marker_Cache mc 
        where r._Object_key = mc._Marker_key 
        and mc._Marker_Type_key = 1 
        and mc._Qualifier_key = 615421 
        union 
        select distinct r.symbol, mc._Marker_key, mc.accID as seqID, mc._LogicalDB_key, mc._Qualifier_key
        from gomarker2 r, SEQ_Marker_Cache mc, MRK_MCV_Cache mcv
        where r._Object_key = mc._Marker_key 
        and mc._Marker_Type_key = 1
        and mc._Qualifier_key = 615420
        and mc._Marker_key = mcv._Marker_key
        and mcv.term = 'miRNA Gene'
        ''', 'auto')
    
    proteins = {}
    
    for r in results:
        key = r['_Marker_key']
        symbol = r['symbol']
        seqID = r['seqID']
        logicalDB = r['_LogicalDB_key']    
        qualifier = r['_Qualifier_key']
    
        # UniProt
        if logicalDB in [13,41]:
            #print('ldb 13/41: ', str(symbol), str(seqID))
            proteins[key] = 'UniProtKB:' + seqID

        # RefSeq
        elif logicalDB in [27] and qualifier == 615421:
            #print('np/xp: ', str(symbol), str(seqID))
            proteins[key] = 'RefSeq:' + seqID

        # Ensembl
        elif logicalDB in [134]:
            proteins[key] = 'ENSEMBL:' + r['seqID']

#
# end doProtein()
#

#
# begin doGAFFinish()
#
def doGAFFinish():

    #
    # process results
    #
    results = db.sql('select * from gomarker2 order by symbol, termID', 'auto')

    for r in results:

        reportRow = ''    

        if r['_Term_key'] not in dag:
            continue

        if dag[r['_Term_key']] not in dagQualifier:
            continue

        objectKey = str(r['_Object_key']) + ':' + str(r['_AnnotEvidence_key'])

        # columns 1-5
        reportRow = MGIPREFIX + TAB
        reportRow = reportRow + str(r['markerID']) + TAB
        reportRow = reportRow + r['symbol'] + TAB

        if r['qualifier'] != None:
            qualifier = r['qualifier'].strip()
        else:
            qualifier = ''

        reportRow = reportRow + qualifier + TAB
        reportRow = reportRow + r['termID'] + TAB

        # column 6; reference
        references = []
        references.append(MGIPREFIX + ':' + r['refID'])
        if r['_Refs_key'] in pubMed:
            references.append('PMID:' + pubMed[r['_Refs_key']])
        else:
            if r['_Refs_key'] in goRefDict:
                references.append(goRefDict[r['_Refs_key']])
        reportRow = reportRow + '|'.join(references) + TAB

        # column 7
        reportRow = reportRow + r['evidenceCode'] + TAB

        # column 8
        inferredFrom = mgi_utils.prvalue(r['inferredFrom']).replace('MGI:', 'MGI:MGI:')
        reportRow = reportRow + inferredFrom + TAB

        # column 9-10
        reportRow = reportRow + dag[r['_Term_key']] + TAB
        reportRow = reportRow + r['name'] + TAB

        # column 11
        if r['_Object_key'] in syns:
            reportRow = reportRow + '|'.join(syns[r['_Object_key']]) + TAB
        else:
            reportRow = reportRow + TAB

        # column 12
        # if marker is associated with an isoform (via go/annotation)
        # or marker is associated with a protein (via marker/sequence cache)
        # 	print('protein')
        # else, print(marker type (ex. 'gene'))

        if objectKey in isoformsProtein or r['_Object_key'] in proteins:
            reportRow = reportRow + 'protein' + TAB
        else:
            reportRow = reportRow + r['markerType'] + TAB
                
        # column 13
        reportRow = reportRow + SPECIES + TAB

        # column 14
        reportRow = reportRow + str(r['mDate']) + TAB

        # column 15; assigned by

        # remove "GOA_"; for example:  "GOA_IntAct" ==> "IntAct"
        # remove "NOCTUA_"; for example:  "NOCTUA_MGI" ==> "MGI"
        if r['assignedBy'].find('NOCTUA_') >= 0:
            assignedBy = r['assignedBy'].replace('NOCTUA_', '')
            reportRow = reportRow + assignedBy + TAB

        elif r['assignedBy'].find('GOA_') >= 0:
            assignedBy = r['assignedBy'].replace('GOA_', '')
            reportRow = reportRow + assignedBy + TAB

        elif r['assignedBy'] in assignedByList1:
            reportRow = reportRow + 'UniProt' + TAB

        elif r['assignedBy'] in assignedByList2:
            reportRow = reportRow + r['assignedBy'] + TAB

        # else use default (MGIPREFIX)
        else:
            reportRow = reportRow + MGIPREFIX + TAB

        #
        # column 16
        # contains property/value information
        # see lib_py_report/go_annot_extensions.py for list of excluded properties
        properties = ''
        if objectKey in gafCol16Lookup:
            properties = ''.join(gafCol16Lookup[objectKey])
        reportRow = reportRow + properties + TAB

        # column 17
        # if isoformProtein = true
        #    then use isoformsProtein
        isoforms = ''
        if objectKey in isoformsProtein:
            isoforms = '|'.join(isoformsProtein[objectKey])
        reportRow = reportRow + isoforms + CRT

        fp.write(reportRow)

        #
        # TR11060
        # subset of UniProtKB:xxxx-?? only
        #
        if objectKey in forPROC:
            fp2.write(reportRow)
 
#
# end doGAFFinish()
#

# begin doGPADFinish()
#
def doGPADFinish():

    #!! 1. DB_Object_ID ::= ID       MGI or PR
    #!! 2. Negation ::= 'NOT'
    #!! 3. Relation ::= OBO_ID
    #!! 4. Ontology_Class_ID ::= OBO_ID
    #!! 5. Reference ::= [ID] ('|' ID)*
    #!! 6. Evidence_type ::= OBO_ID
    #!! 7. With_or_From ::= [ID] ('|' | ‘,’ ID)*
    #!! 8. Interacting_taxon_ID ::= NCBITaxon:[Taxon_ID]
    #!! 9. Date ::= YYYY-MM-DD
    #!! 10. Assigned_by ::= Prefix
    #!! 11. Annotation_Extensions ::= [Extension_Conj] ('|' Extension_Conj)*
    #!! 12. Annotation_Properties ::= [Property_Value_Pair] ('|' Property_Value_Pair)*

    #
    # process results
    #
    results = db.sql('select * from gomarker2 order by symbol, termID', 'auto')

    for r in results:

        reportRow = ''    

        if r['_Term_key'] not in dag:
            continue

        if dag[r['_Term_key']] not in dagQualifier:
            continue

        #
        #  1. DB_Object_ID ::= ID       MGI or PR
        #
        # if an Isoform (gene_product) exists, then create the annotation using:
        # 	isoformsProtein prefix (PR, RefSeq, UniProtDB, EMBL) ":" accession id  of isoformProtein object (Q92WPO-1)
        # else create the annotation using:
        # 	MGI:MGI:xxx
        #

        objectKey = str(r['_Object_key']) + ':' + str(r['_AnnotEvidence_key'])

        if objectKey in isoformsProtein:
            for i in isoformsProtein[objectKey]:
                try:
                        tokens = i.split(':')
                        prefixPart = tokens[0]
                        numericPart = tokens[1]
                        reportRow = tokens[0] + ':' + tokens[1] + TAB
                        reportRow = addGPADReportRow(reportRow, r)
                except:
                        reportRow = MGIPREFIX + TAB + str(r['markerID']) + TAB
                        reportRow = addGPADReportRow(reportRow, r)
        else:
            reportRow = MGIPREFIX + ':' + str(r['markerID']) + TAB
            reportRow = addGPADReportRow(reportRow, r)

        fp.write(reportRow)

#
# end doGPADFinish()
#

# begin addGPADReportRow(reportRow, r)
#
def addGPADReportRow(reportRow, r):

        objectKey = str(r['_Object_key']) + ':' + str(r['_AnnotEvidence_key'])
        key = r['_AnnotEvidence_key']

        #
        # for col 2 Negation and col 3 Relation
        #
        # MGI-Qualifier : VOC_Annot._qualifier_key (_vocab_key = 52)
        #
        # GO-Properties : VOC_EvidenceProperty (_vocab_key = 82) (goPropertyLookup)
        #
        # GO-Qualifier  : VOC_EvidenceProperty where value = 'go_qualifier' (goQualifierLookup)
        #
        # DAG-Qualifier : hard-coded RO of C, P, F (dagQualifier)
        #

        #
        # 2. Negation ::= "NOT"
        # if MGI-Qualifier contains "NOT", then attach the term "NOT" to col 2
        # if MGI-Qualifier = "NOT|xxxx"
        #       then find the RO id for "xxxx" and attach to col 3 (goPropertyLookup)
        #
        default_relation = ''
        if r['qualifier'] != None:
            tokens = r['qualifier'].split('|')
            try:
                qualifier = 'NOT'
                property = tokens[1]
                if property in goPropertyLookup:
                    default_relation = goPropertyLookup[property][0];
            except:
                if tokens[0] == 'NOT':
                    qualifier = 'NOT'
                else:
                    qualifier = ''
        else:
            qualifier = ''
        reportRow = reportRow + qualifier + TAB

        #
        # 3. Relation ::= OBO_ID
        # if col 3 was set by previous col 2 processing then done with col 3
        # else if GO Property exists in goQualifierLookup, then use its RO id
        # else if inferredFrom contains "InterPro", then use RO:0002331
        # else use dagQualifier (C, P, F)

        if default_relation == '':
            if key in goQualifierLookup:
                default_relation = '|'.join(goQualifierLookup[key])
            elif r['inferredFrom'] != None and r['inferredFrom'].find('InterPro:') >= 0 and dag[r['_Term_key']] == 'P':
                default_relation = 'RO:0002331'
            else:
                default_relation = dagQualifier[dag[r['_Term_key']]]
        reportRow = reportRow + default_relation + TAB

        # 4. Ontology_Class_ID ::= OBO_ID/GO ID
        reportRow = reportRow + r['termID'] + TAB

        # 5. Reference ::= [ID] ("|" ID)*
        references = []
        references.append(MGIPREFIX + ':' + r['refID'])
        if r['_Refs_key'] in pubMed:
            references.append('PMID:' + pubMed[r['_Refs_key']])
        else:
            if r['_Refs_key'] in goRefDict:
                references.append(goRefDict[r['_Refs_key']])
        reportRow = reportRow + '|'.join(references) + TAB

        # 6. Evidence_type ::= OBO_ID
        if key in evidenceLookup:
            reportRow = reportRow + evidenceLookup[key][0]
        elif r['evidenceCode'] in ecoLookupByEvidence:
            reportRow = reportRow + ecoLookupByEvidence[r['evidenceCode']]
        else:
            reportRow = reportRow + 'NOT FOUND'
        reportRow = reportRow + TAB

        # 7. With_or_From ::= [ID] ("|" | "," ID)*
        inferredFrom = mgi_utils.prvalue(r['inferredFrom']).replace('MGI:', 'MGI:MGI:')
        reportRow = reportRow + mgi_utils.prvalue(inferredFrom) + TAB

        # 8. Interacting_taxon_ID ::= NCBITaxon:[Taxon_ID]
        if key in taxonLookup:
            reportRow = reportRow + taxonLookup[key][0].replace('taxon', 'NCBITaxon')
        reportRow = reportRow + TAB

        # 9. Date ::= YYYY-MM-DD
        reportRow = reportRow + str(r['mDate']) + TAB

        # 10. Assigned_by ::= Prefix

        # remove "NOCTUA_"; for example:  "NOCTUA_MGI" ==> "MGI"
        if r['assignedBy'].find('NOCTUA_') >= 0:
            assignedBy = r['assignedBy'].replace('NOCTUA_', '')
            reportRow = reportRow + assignedBy + TAB

        # remove "GOA_"; for example:  "GOA_IntAct" ==> "IntAct"
        elif r['assignedBy'].find('GOA_') >= 0:
            assignedBy = r['assignedBy'].replace('GOA_', '')
            reportRow = reportRow + assignedBy + TAB

        elif r['assignedBy'] in assignedByList1:
            reportRow = reportRow + 'UniProt' + TAB

        elif r['assignedBy'] in assignedByList2:
            reportRow = reportRow + r['assignedBy'] + TAB

        # else use default (MGIPREFIX)
        else:
            reportRow = reportRow + MGIPREFIX + TAB

        # 11. Annotation_Extensions ::= [Extension_Conj] ("|" Extension_Conj)*
        properties = ''
        if key in gpadCol11Lookup:
            properties = ','.join(gpadCol11Lookup[key])
        elif objectKey in gafCol16Lookup:
            properties = ''.join(gafCol16Lookup[objectKey])
        reportRow = reportRow + properties + TAB

        # 12. Annotation_Properties ::= [Property_Value_Pair] ("|" Property_Value_Pair)*\n')
        #
        properties = ''
        if key in gpadCol12Lookup:
            properties = '|'.join(gpadCol12Lookup[key])
        reportRow = reportRow + properties + CRT

        return reportRow

#
# end addGPADReportRow()
#

#
# main
#

db.useOneConnection(1)

# querying the database/setting up lookups
doSetup()
doGAFCol16()
doIsoform()

#
# GAF 2.0
#

fp = reportlib.init('gene_association', fileExt = '.mgi2', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp2 = reportlib.init('gene_association_pro', fileExt = '.mgi2', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('!gaf-version: 2.1\n')
fp.write('!software version: $Revision$\n')
fp.write('!date: %s $\n' % (mgi_utils.date("%m/%d/%Y")))
fp.write('!\n')
fp.write('! from Mouse Genome Database (MGD) & Gene Expression Database (GXD)\n')
fp.write('!\n')
fp.write('!1.  DB                              required        1               UniProtKB\n')
fp.write('!2.  DB Object ID                    required        1               P12345\n')
fp.write('!3.  DB Object Symbol                required        1               PHO3\n')
fp.write('!4.  Qualifier                       optional        0 or greater    NOT\n')
fp.write('!5.  GO ID                           required        1               GO:0003993\n')
fp.write('!6.  DB:Reference (|DB:Reference)    required        1 or greater    PMID:2676709\n')
fp.write('!7.  Evidence Code                   required        1               IMP\n')
fp.write('!8.  With (or) From                  optional        0 or greater    GO:0000346\n')
fp.write('!9.  Aspect                          required        1               F\n')
fp.write('!10. DB Object Name                  optional        0 or 1          Toll-like receptor 4\n')
fp.write('!11. DB Object Synonym (|Synonym)    optional        0 or greater    hToll   Tollboot\n')
fp.write('!12. DB Object Type                  required        1               protein\n')
fp.write('!13. Taxon(|taxon)                   required        1 or 2          taxon:9606\n')
fp.write('!14. Date                            required        1               20090118\n')
fp.write('!15. Assigned By                     required        1               SGD\n')
fp.write('!16. Annotation Extension            optional        0 or greater    part_of(CL:0000576)\n')
fp.write('!17. Gene Product Form ID            optional        0 or 1          UniProtKB:P12345-2\n')
fp.write('!\n')

doProtein()
doGAFFinish()

# append GOA annotations, if exists : see goload/goamouse
try:
    goaFile = open(os.environ['GOAGAFMGI'], 'r')
    for line in goaFile.readlines():
        fp.write(line)
    goaFile.close()
except:
    pass

reportlib.finish_nonps(fp)
reportlib.finish_nonps(fp2)

#
# GPAD 2.0
#

fp = reportlib.init('mgi2', fileExt = '.gpad', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('!gpa-version: 2.0\n') 
fp.write('!\n')
fp.write('!date: %s $\n' % (mgi_utils.date("%m/%d/%Y")))
fp.write('!\n')
fp.write('! from Mouse Genome Database (MGD) & Gene Expression Database (GXD)\n')
fp.write('!\n')
fp.write('!1.  DB_Object_ID            ::= ID MGI or PR\n')
fp.write('!2.  Negation                ::= "NOT"\n')
fp.write('!3.  Relation                ::= OBO_ID\n')
fp.write('!4.  Ontology_Class_ID       ::= OBO_ID\n')
fp.write('!5.  Reference               ::= [ID] ("|" ID)*\n')
fp.write('!6.  Evidence_type           ::= OBO_ID\n')
fp.write('!7.  With_or_From            ::= [ID] ("|" | "," ID)*\n')
fp.write('!8.  Interacting_taxon_ID    ::= NCBITaxon:[Taxon_ID]\n')
fp.write('!9.  Date                    ::= YYYY-MM-DD\n')
fp.write('!10. Assigned_by             ::= Prefix\n')
fp.write('!11. Annotation_Extensions   ::= [Extension_Conj] ("|" Extension_Conj)*\n')
fp.write('!12. Annotation_Properties   ::= [Property_Value_Pair] ("|" Property_Value_Pair)*\n')
fp.write('!\n')

doGPADFinish()

# append GOA annotations, if exists : see goload/goamouse
# CONVERT
#try:
#    goafile = open(os.environ['GOAGPADMGI'], 'r')
#    for line in goafile.readlines():
#        fp.write(line)
#    goafile.close()
#except:
#    pass

reportlib.finish_nonps(fp)

db.useOneConnection(0)
