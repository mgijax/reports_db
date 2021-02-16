'''
#
# GO_gene_association.py
#
# Generates:
# intentionally creating gaf & gpad using one script
# because they use much of the same logic/lookups/etc
#       
#	gene_association.mgi (GAF)
#	gene_association_pro.mgi (GAF)
#	mgi.gpad
#
# ALSO CHANGE:  weekly/GO_gene_association_nonmouse.py
# 09/11/2020 per David Hill, do nothing until we have Tues/09/15 meeting
#
# GAF 2.2 : see below
# fp.write('!gaf-version: 2.2\n')
#
# GPAD 1.1 : see below
# fp.write('!gpa-version: 1.1\n') 
#
# IMPORTANT THINGS TO KNOW:
#
#    gaf/col11 (annotation extension) and gpad/col16 should be equal.
#    except :  ISO annotations are excluded from gaf/col11 (blank)
#    see "lib_py_report/go_annot_extensions.py" for the list of
#	excluded Properties and excluded Evidence
#
#    gpad/col16 : will *never* contains > 1 stanza, and will always use the "," delimiter
#
# lec   08/25/2020
#       - TR13272/converting to GPI 2.0
#       mgi2.gpad : dph reviewing 09/10/2020
#       gene_association.mgi2 : no changes yet
#       gene_association_pro.mgi2 : no changes yet
#
# For GAF2.2
# If the annotation has the value GO:0008150 (biological_process) in column 5, 
#       then the value in column 4 should be involved_in.
# If the annotation has the value GO:0003674 (molecular_function) in column 5, 
#       then the value in column 4 should be enables.
# If the annotation has the value GO:0005575 (cellular_component) in column 5, 
#       then the value in column 4 should be is_active_in.
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

TAB = reportlib.TAB
CRT = reportlib.CRT

MGIPREFIX = 'MGI'

#
# if in list 1, then use 'UniProt'
# if in list 2, then use itself
# else it will use MGIPREFIX
#
assignedByList1 = ['uniprotload']
assignedByList2 = ['GOC', 'GO_Central']

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

# translate dag to qualifier
dagQualifierGAF = {'C':'located_in', 'P':'acts_upstream_of_or_within', 'F':'enables'}
goQualifierGAF = {}

# located_in, acts_upstream_of_or_within, enables
# GPAD 1.1
dagQualifierGPAD = {'C':'part_of', 'P':'acts_upstream_of_or_within', 'F':'enables'}
# GPAD 2.0
#dagQualifierGPAD = {'C':'RO:0001025', 'P':'RO:0002264', 'F':'RO:0002327'}
goQualifierGPAD = {}

#
# gpad lookups
#
taxonLookup = {}
ecoLookupByEco = {}
ecoLookupByEvidence = {}
evidenceLookup = {}
goPropertyLookup = {}
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
    global goQualifierGAF
    global goQualifierGPAD
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
    #   and m.symbol = 'Asap1'
    #   and m.symbol = 'Mbd2'
    #   and m.symbol = 'Adipoq'
    #   and m.symbol = 'Birc3'
    #	and m.symbol = 'Hk1'
    #
    # for Dustin, exclude:
    # _refs_key |    mgiid    |  jnumid  |                       short_citation
    # -----------+-------------+----------+-------------------------------------------------------------
    #     156949 | MGI:4417868 | J:155856 | Mouse Genome Informatics Scientific Curators,  2010 Jan;():
    #     165659 | MGI:4834177 | J:164563 | Mouse Genome Informatics Scientific Curators,  2010 Oct;():

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
             VOC_Term q,
             MGI_User u
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
        and e._ModifiedBy_key = u._User_key
        -- for Dustin/only MGI_curated
        --and u.orcid is not null
        and e._Refs_key not in (156949, 165659)
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

        -- for Dustin/only MGI_curated
        --and u.orcid is not null

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
        value = r['value'].replace('taxon', 'NCBITaxon')
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
    # goQualifierGAF : go_qualifier
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

        value = r['value']
        if key not in goQualifierGAF:
            goQualifierGAF[key] = []
        goQualifierGAF[key].append(value)

    #print(goQualifierGAF)

    #
    # goQualifierGPAD: go_qualifier
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
        if key not in goQualifierGPAD:
            goQualifierGPAD[key] = []
        goQualifierGPAD[key].append(value)
    #print goQualifierGPAD

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
    #print gpadCol11Lookup

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
    #print gpadCol12Lookup

    results = db.sql('select _Object_key, accID from ACC_Accession where _LogicalDB_key = 185', 'auto')
    for r in results:
        key = r['_Object_key']
        value = r['accID']
        goRefDict[key] = value
    #print goRefDict

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

    #print(gafCol16Lookup)

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

        if dag[r['_Term_key']] not in dagQualifierGAF:
            continue

        objectKey = str(r['_Object_key']) + ':' + str(r['_AnnotEvidence_key'])
        key = r['_AnnotEvidence_key']

        #!1  DB               
        #!2  DB Object ID    
        #!3  DB Object Symbol
        reportRow = MGIPREFIX + TAB
        reportRow = reportRow + str(r['markerID']) + TAB
        reportRow = reportRow + r['symbol'] + TAB

        #
        # GO-Qualifier : value where GO-Property = “go_qualifier”
        # For example : GO-Property = “go_qualifier”, value = “part_of”
        #
        # MGI-Qualifier :  MGI original qualifier list before GO started
        #       empty
        #       NOT
        #       colocalizes_with
        #       NOT|colocalizes_with
        #       contributes_to
        #       NOT|contributes_to
        #
        # Default Qualifier:
        # {'C':'located_in', 'P':'acts_upstream_of_or_within', 'F':'enables'
        #
        # If the annotation has the value GO:0008150 (biological_process) in column 5, 
        #       then the value in column 4 should be involved_in.
        # else if the annotation has the value GO:0003674 (molecular_function) in column 5, 
        #       then the value in column 4 should be enables.
        # else if the annotation has the value GO:0005575 (cellular_component) in column 5, 
        #       then the value in column 4 should be is_active_in.
        # else if the Annotation contains a GO-Qualifier, then use the “go_qualifier”/value
        #       If MGI-Qualifier = NOT, then attach NOT to front of GO-Qualifier
        #       For example :  NOT|enables
        # else if the MGI-Qualifier = "NOT, then "NOT:" + Default Qualifier
        # else if the MGI-Qualifier is not empty, then use the MGI-Qualifier value
        # else if the Annotation/Inferred From contains “InterPro:”, then use “involved_in”
        # else use Default Qualifier
        #

        #!4  Qualifier      
        qualifier = ""
        if r['termID'] == 'GO:0008150':
                qualifier = 'involved_in'
        elif r['termID'] == 'GO:0003674':
                qualifier = 'enables'
        elif r['termID'] == 'GO:0005575':
                qualifier = 'is_active_in'
        elif key in goQualifierGAF:
            if r['qualifier'] == 'NOT':
               qualifier = 'NOT|'
            qualifier = qualifier + '|'.join(goQualifierGAF[key])
        elif r['qualifier'] == 'NOT':
            qualifier = 'NOT|' + dagQualifierGAF[dag[r['_Term_key']]]
        elif r['qualifier'] != None:
            qualifier = r['qualifier'].strip()
        elif r['inferredFrom'] != None and r['inferredFrom'].find('InterPro:') >= 0 and dag[r['_Term_key']] == 'P':
            qualifier = 'involved_in'
        else:
            qualifier = dagQualifierGAF[dag[r['_Term_key']]]
        reportRow = reportRow + qualifier + TAB

        #!5  GO ID         
        reportRow = reportRow + r['termID'] + TAB

        #!6  DB:Reference (|DB:Reference) 
        references = []
        references.append(MGIPREFIX + ':' + r['refID'])
        if r['_Refs_key'] in pubMed:
            references.append('PMID:' + pubMed[r['_Refs_key']])
        else:
            if r['_Refs_key'] in goRefDict:
                references.append(goRefDict[r['_Refs_key']])
        reportRow = reportRow + '|'.join(references) + TAB

        #!7  Evidence Code               
        reportRow = reportRow + r['evidenceCode'] + TAB

        #!8  With (or) From             
        inferredFrom = mgi_utils.prvalue(r['inferredFrom']).replace('MGI:', 'MGI:MGI:')
        reportRow = reportRow + inferredFrom + TAB

        #!9  Aspect                     
        #!10 DB Object Name             
        reportRow = reportRow + dag[r['_Term_key']] + TAB
        reportRow = reportRow + r['name'] + TAB

        #!11 DB Object Synonym (|Synonym)
        if r['_Object_key'] in syns:
            reportRow = reportRow + '|'.join(syns[r['_Object_key']]) + TAB
        else:
            reportRow = reportRow + TAB

        #!12 DB Object Type             
        # if marker is associated with an isoform (via go/annotation)
        # or marker is associated with a protein (via marker/sequence cache)
        # 	print('protein')
        # else, print(marker type (ex. 'gene'))

        if objectKey in isoformsProtein or r['_Object_key'] in proteins:
            reportRow = reportRow + 'protein' + TAB
        else:
            reportRow = reportRow + r['markerType'] + TAB
                
        #!13 Taxon(|taxon)             
        reportRow = reportRow + 'taxon:10090' + TAB

        #!14 Date                     
        reportRow = reportRow + str(r['mDate']) + TAB

        #!15 Assigned By             
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

        #!16 Annotation Extension   
        # column 16
        # contains property/value information
        # see lib_py_report/go_annot_extensions.py for list of excluded properties
        properties = ''
        if objectKey in gafCol16Lookup:
            properties = ''.join(gafCol16Lookup[objectKey])
        reportRow = reportRow + properties + TAB

        #!17 Gene Product Form ID  
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

    # Output format:
    #
    #   1. DB                       MGI
    #   2. DB Object ID             MGI:xxxx
    #   3. Qualifier    
    #   4. GO ID                    GO:xxxx
    #   5. DB:Reference(s)          MGI:MGI:xxxx|PMID:xxxx
    #   6. Evidence Code            ECO:xxxx
    #   7. With (or)From            (optional)
    #   8. Interacting taxon ID     (optional)
    #   9. Date                     YYYYMMDD
    #   10. Assigned by
    #   11. Annotation Extension    (optional) same as GAF/col 16
    #   12. Annotation Properties   (optional) properties if creator like 'NOCTUA_%'

    #
    # process results
    #
    results = db.sql('select * from gomarker2 order by symbol, termID', 'auto')

    for r in results:

        reportRow = ''    

        if r['_Term_key'] not in dag:
            continue

        if dag[r['_Term_key']] not in dagQualifierGPAD:
            continue

        #
        #   1. DB
        #   2. DB Object ID
        #
        # if an Isoform (gene_product) exists, then create the annotation using:
        # 	col 1 =  isoformsProtein prefix (PR, RefSeq, UniProtDB, EMBL)
        # 	col 2 =  accession id  of isoformProtein object (Q92WPO-1)
        # else create the annotation using:
        # 	col 1 =  MGI
        # 	col 2 =  MGI:xxx
        #

        objectKey = str(r['_Object_key']) + ':' + str(r['_AnnotEvidence_key'])

        if objectKey in isoformsProtein:
            for i in isoformsProtein[objectKey]:
                try:
                        tokens = i.split(':')
                        prefixPart = tokens[0]
                        numericPart = tokens[1]
                        reportRow = tokens[0] + TAB + tokens[1] + TAB
                        reportRow = addGPADReportRow(reportRow, r)
                except:
                        reportRow = MGIPREFIX + TAB + str(r['markerID']) + TAB
                        reportRow = addGPADReportRow(reportRow, r)
        else:
            reportRow = MGIPREFIX + TAB + str(r['markerID']) + TAB
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

        #   3. Qualifier

        if key in goQualifierGPAD:
            default_relation_for_aspect = '|'.join(goQualifierGPAD[key])
        elif r['inferredFrom'] != None and r['inferredFrom'].find('InterPro:') >= 0 and dag[r['_Term_key']] == 'P':
            default_relation_for_aspect = 'involved_in'
        else:
            default_relation_for_aspect = dagQualifierGPAD[dag[r['_Term_key']]]

        # qualifier from MGD annotations
        if r['qualifier'] != None:
            qualifier = r['qualifier'].strip()
        else:
            qualifier = ''

        if qualifier == '':
            gap_qualifier = default_relation_for_aspect
        elif qualifier == 'NOT':
            gap_qualifier = qualifier + '|' + default_relation_for_aspect
        else:
            gap_qualifier = qualifier

        reportRow = reportRow + gap_qualifier + TAB

        #   4. GO ID
        reportRow = reportRow + r['termID'] + TAB

        #   5. DB:Reference(s)
        references = []
        references.append(MGIPREFIX + ':' + r['refID'])
        if r['_Refs_key'] in pubMed:
            references.append('PMID:' + pubMed[r['_Refs_key']])
        else:
            if r['_Refs_key'] in goRefDict:
                references.append(goRefDict[r['_Refs_key']])
        reportRow = reportRow + '|'.join(references) + TAB

        #   6. Evidence Code
        if key in evidenceLookup:
            reportRow = reportRow + evidenceLookup[key][0]
        elif r['evidenceCode'] in ecoLookupByEvidence:
            reportRow = reportRow + ecoLookupByEvidence[r['evidenceCode']]
        else:
            reportRow = reportRow + 'NOT FOUND'
        reportRow = reportRow + TAB

        #   7. With (or)From
        inferredFrom = mgi_utils.prvalue(r['inferredFrom']).replace('MGI:', 'MGI:MGI:')
        reportRow = reportRow + mgi_utils.prvalue(inferredFrom) + TAB

        #   8. Interacting taxon ID
        if key in taxonLookup:
            reportRow = reportRow + taxonLookup[key][0]
        reportRow = reportRow + TAB

        #   9. Date
        reportRow = reportRow + str(r['mDate']) + TAB

        #   10. Assigned by

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

        #   11. Annotation Extension
        properties = ''
        if key in gpadCol11Lookup:
            properties = ','.join(gpadCol11Lookup[key])
        elif objectKey in gafCol16Lookup:
            properties = ''.join(gafCol16Lookup[objectKey])
        reportRow = reportRow + properties + TAB

        #   12. Annotation Properties
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
# GAF 2.2
#
fp = reportlib.init('gene_association', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp.write('!gaf-version: 2.2\n')
fp.write('!generated-by: MGI\n')
fp.write('!date-generated: %s\n' % (mgi_utils.date("%Y-%m-%d")))
fp.write('!\n')
fp.write('!1  DB                              required        1       UniProtKB\n')
fp.write('!2  DB Object ID                    required        1       P12345\n')
fp.write('!3  DB Object Symbol                required        1       PHO3\n')
fp.write('!4  Qualifier                       required        1 or 2  NOT|involved_in\n')
fp.write('!5  GO ID                           required        1       GO:0003993\n')
fp.write('!6  DB:Reference (|DB:Reference)    required        1 or greater    PMID:2676709\n')
fp.write('!7  Evidence Code                   required        1       IMP\n')
fp.write('!8  With (or) From                  optional        0 or greater    GO:0000346\n')
fp.write('!9  Aspect                          required        1       F\n')
fp.write('!10 DB Object Name                  optional        0 or 1  Toll-like receptor 4\n')
fp.write('!11 DB Object Synonym (|Synonym)    optional        0 or greater    hToll\n')
fp.write('!12 DB Object Type                  required        1       protein\n')
fp.write('!13 Taxon(|taxon)                   required        1 or 2  taxon:9606\n')
fp.write('!14 Date                            required        1       20090118\n')
fp.write('!15 Assigned By                     required        1       SGD\n')
fp.write('!16 Annotation Extension            optional        0 or greater    part_of(CL:0000576)\n')
fp.write('!17 Gene Product Form ID            optional        0 or 1  UniProtKB:P12345-2\n')
fp.write('!\n')

fp2 = reportlib.init('gene_association_pro', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp2.write('!gaf-version: 2.2\n')

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
# GPAD 1.1
#

fp = reportlib.init('mgi', fileExt = '.gpad', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('!gpa-version: 1.1\n') 
fp.write('!generated-by: MGI\n')
fp.write('!date: %s $\n' % (mgi_utils.date("%m/%d/%Y")))
fp.write('!\n')
fp.write('!1.  DB                       MGI or PR\n')
fp.write('!2.  DB Object ID             MGI:xxxx or xxxxx\n')
fp.write('!3.  Qualifier                enables, causally_upstream_of_or_within, part_of\n')
fp.write('!4.  GO ID                    GO:xxxx\n')
fp.write('!5.  DB:Reference(s)          MGI:MGI:xxxx|PMID:xxxx\n')
fp.write('!6.  Evidence Code            ECO:xxxx\n')
fp.write('!7.  With (or)From            optional\n')
fp.write('!8.  Interacting taxon ID     optional\n')
fp.write('!9.  Date                     YYYYMMDD\n')
fp.write('!10. Assigned By\n')
fp.write('!11. Annotation Extension     optional same as GAF/col 16\n')
fp.write('!12. Annotation Properties    optional\n')
fp.write('!\n')

doGPADFinish()

# append GOA annotations, if exists : see goload/goamouse
try:
    goafile = open(os.environ['GOAGPADMGI'], 'r')
    for line in goafile.readlines():
        fp.write(line)
    goafile.close()
except:
    pass

doGPADFinish()

# append GOA annotations, if exists : see goload/goamouse
try:
        goafile = open(os.environ['GOAGPADMGI'], 'r')
        for line in goafile.readlines():
                fp.write(line)
        goafile.close()
except:
        pass

reportlib.finish_nonps(fp)

db.useOneConnection(0)
