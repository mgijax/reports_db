'''
#
# GO_gene_association_2.0.py
#
# Generates:
#       
#	gene_association.mgi (GAF)
#	gene_association_pro (GAF)
#	mgi.gpa
#
#
# gpad col vs gaf col
#
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
# DB			1	1
# DB_Object_ID		2	2/17 (isoforms)
# Qualifier		3	4
# GO ID			4	5
# DB:Reference(s)	5	6
# Evidence code		6	7
# With (or) From	7	8
# Interacting taxon ID	8	13
# Date			9	14
# Assigned_by		10	15
# Annotation Extension	11	16
# Annotation Properties	12
# 
# lec   08/25/2020
#       - TR13272/converting to GPI 2.0
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
import rolib

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
dagQualifier = {'C':'part_of', 'P':'acts_upstream_of_or_within', 'F':'enables'}
taxonLookup = {}
ecoLookupByEco = {}
ecoLookupByEvidence = {}
evidenceLookup = {}
gpadCol3Lookup = {}
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
    global roLookupByRO, roLookupByName
    global gpadCol3Lookup
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
            e._AnnotEvidence_key, e.inferredFrom, e.modification_date, e._EvidenceTerm_key, 
            e._Refs_key, e._CreatedBy_key, e._ModifiedBy_key, 
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
            to_char(g.modification_date, 'YYYY-MM-DD') as mDate,
            g._CreatedBy_key,
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

    #
    # roLookupByRO
    # roLookupByName
    #
    # use goload/rolib to lookup RO using ...
    roLookupByRO, roLookupByName = rolib.processRO()
    #print(roLookupByName)

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
    # gpadCol3 : go_qualifier
    #
    results = db.sql('''select distinct a._AnnotEvidence_key, t.term, p.value
            from gomarker2 a,
                 VOC_Evidence_Property p,  
                 VOC_Term t
            where a._AnnotEvidence_key = p._AnnotEvidence_key
            and p._PropertyTerm_key = t._Term_key
            and t.term in ('go_qualifier')
            order by t.term, p.value
            ''', 'auto')
    for r in results:
        key = r['_AnnotEvidence_key']
        value = r['value']
        if key not in gpadCol3Lookup:
            gpadCol3Lookup[key] = []
        gpadCol3Lookup[key].append(value)
    #print(gpadCol3Lookup)

    #
    # gpadCol11 : (MGI_User.login like NOCTUA_%)
    #	exclude older terms (sequenceNum 1-9, 90,91,92, 93)
    #
    # note that noctua-generated properties will *always* have one stanza
    #
    results = db.sql('''select distinct a._AnnotEvidence_key, t.term, p.value
            from gomarker2 a,
                 VOC_Evidence_Property p,  
                 VOC_Term t,
                 MGI_User u
            where a._CreatedBy_key = u._User_key
            and u.login like 'NOCTUA_%'
            and a._AnnotEvidence_key = p._AnnotEvidence_key
            and p._PropertyTerm_key = t._Term_key
            and t.term not in (
                'evidence', 'anatomy', 'cell type', 'gene product', 'modification', 'target', 
                'external ref', 'text', 'dual-taxon ID',
                'noctua-model-id', 'contributor', 'individual', 'go_qualifier', 'model-state'
                )
            order by t.term, p.value
            ''', 'auto')
    for r in results:
        key = r['_AnnotEvidence_key']
        value = r['value'].replace('MGI:', 'MGI:MGI:')
        value = r['term'] + '(' + value + ')'
        if key not in gpadCol11Lookup:
            gpadCol11Lookup[key] = []
        gpadCol11Lookup[key].append(value)
    #print(gpadCol11Lookup)

    #
    # gpadCol12 : (MGI_User.login like NOCTUA_%)
    # exclude : occurs_in, part_of, go_qualifier, evidence
    #
    # TR13272
    # and (u.login like 'NOCTUA_%' or (u.orcid is not null and p._propertyterm_key = 18583062))
    #
    results = db.sql('''select distinct a._AnnotEvidence_key, t.term, p.value
            from gomarker2 a,
                 VOC_Evidence_Property p,  
                 VOC_Term t,
                 MGI_User u
            where a._CreatedBy_key = u._User_key
            and (u.login like 'NOCTUA_%' or (u.orcid is not null and p._propertyterm_key = 18583062))
            and a._AnnotEvidence_key = p._AnnotEvidence_key
            and p._PropertyTerm_key = t._Term_key
            and t.term not in ('occurs_in', 'part_of', 'go_qualifier', 'evidence')
            order by t.term, p.value
            ''', 'auto')
    for r in results:
        key = r['_AnnotEvidence_key']
        value = r['term'] + '=' + r['value']
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
    select r.symbol, r._Object_key, r._AnnotEvidence_key, t.term as property, p.value, p.stanza
            from gomarker2 r, VOC_Evidence_Property p, VOC_Term t
            where r._EvidenceTerm_key in (%s)
            and r._AnnotEvidence_key = p._AnnotEvidence_key
            and p._PropertyTerm_key = t._Term_key
            and p._PropertyTerm_key in (%s)
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

        # 2. Negation ::= "NOT"
        # qualifier from MGD annotations
        if r['qualifier'] != None:
            qualifier = r['qualifier'].strip()
        else:
            qualifier = ''
        reportRow = reportRow + qualifier + TAB

        # 3. Relation ::= OBO_ID
        # CONVERTED TO RO:xxx
        # use gadCol3 or DAG
        if key in gpadCol3Lookup:
            default_relation_for_aspect = '|'.join(gpadCol3Lookup[key])
        elif r['inferredFrom'] != None and r['inferredFrom'].find('InterPro:') >= 0 and dag[r['_Term_key']] == 'P':
            default_relation_for_aspect = 'involved_in'
        else:
            default_relation_for_aspect = dagQualifier[dag[r['_Term_key']]]

        # replace "_" with space
        default_relation_for_aspect = default_relation_for_aspect.replace('_', ' ')

        try:
                roId = roLookupByName[default_relation_for_aspect][0]
        except:
                roId = 'cannot find RO: ' + default_relation_for_aspect

        reportRow = reportRow + roId + TAB

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
        # CONVERTED TO : NCBITaxon:1280
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
        # CONVERT properties to RO
        properties = ''
        if key in gpadCol11Lookup:
            properties = ','.join(gpadCol11Lookup[key])
        elif objectKey in gafCol16Lookup:
            properties = ''.join(gafCol16Lookup[objectKey])
        reportRow = reportRow + properties + TAB

        # 12. Annotation_Properties ::= [Property_Value_Pair] ("|" Property_Value_Pair)*\n')
        # CONVERT
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
# GAF
#

fp = reportlib.init('gene_association', fileExt = '.mgi2', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp2 = reportlib.init('gene_association_pro', fileExt = '.mgi2', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('!gaf-version: 2.1\n')
fp.write('!software version: $Revision$\n')
fp.write('!date: %s $\n' % (mgi_utils.date("%m/%d/%Y")))
fp.write('!\n')
fp.write('! from Mouse Genome Database (MGD) & Gene Expression Database (GXD)\n')
fp.write('!\n')
fp.write('!! 1. DB_Object_ID              MGI or PR\n')
fp.write('!! 2. Negation                  NOT\n')
fp.write('!! 3. Relation                  OBO_ID\n')
fp.write('!! 4. Ontology_Class_ID         OBO_ID\n')
fp.write('!! 5. Reference                 ID|ID\n')
fp.write('!! 6. Evidence_type             OBO_ID\n')
fp.write('!! 7. With_or_From              [ID] ("|" | "," ID)\n')
fp.write('!! 8. Interacting_taxon_ID      NCBITaxon:[Taxon_ID]\n')
fp.write('!! 9. Date                      YYYY-MM-DD\n')
fp.write('!! 10. Assigned_by              Prefix\n')
fp.write('!! 11. Annotation_Extensions    [Extension_Conj] ("|" Extension_Conj)\n')
fp.write('!! 12. Annotation_Properties    [Property_Value_Pair] ("|" Property_Value_Pair)\n')
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
fp.write('! 1.  DB_Object_ID            ::= ID       MGI or PR\n')
fp.write('! 2.  Negation                ::= "NOT"\n')
fp.write('! 3.  Relation                ::= OBO_ID\n')
fp.write('! 4.  Ontology_Class_ID       ::= OBO_ID\n')
fp.write('! 5.  Reference               ::= [ID] ("|" ID)*\n')
fp.write('! 6.  Evidence_type           ::= OBO_ID\n')
fp.write('! 7.  With_or_From            ::= [ID] ("|" | "," ID)*\n')
fp.write('! 8.  Interacting_taxon_ID    ::= NCBITaxon:[Taxon_ID]\n')
fp.write('! 9.  Date                    ::= YYYY-MM-DD\n')
fp.write('! 10. Assigned_by             ::= Prefix\n')
fp.write('! 11. Annotation_Extensions   ::= [Extension_Conj] ("|" Extension_Conj)*\n')
fp.write('! 12. Annotation_Properties   ::= [Property_Value_Pair] ("|" Property_Value_Pair)*\n')
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
