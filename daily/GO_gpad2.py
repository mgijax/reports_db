'''
#
# GO_gpad2.py
#
# Generates:
#       
#	mgi2.gpad
#
# GPAD 2.0 : see below
# fp.write('!gpa-version: 2.0\n') 
#
# lec   08/25/2020
#       - TR13272/converting to GPI 2.0
#       mgi2.gpad : dph reviewing 09/10/2020
#
'''

import sys
import os
import mgi_utils
import reportlib
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
pubMed = {}

# mapping of load reference keys to GO_REF IDs/per TR11772
# see _LogicalDB_key = 185 (GO_REF)
goRefLookup = {}

# see doIsoform() : isoformsProtein = {}
forPROC = {}

# located_in, acts_upstream_of_or_within, enables
dagQualifier = {'C':'RO:0001025', 'P':'RO:0002264', 'F':'RO:0002327'}

#
# gpad lookups
#
taxonLookup = {}
ecoLookupByEco = {}
ecoLookupByEvidence = {}
evidenceLookup = {}
goPropertyLookup = {}
goQualifierLookup = {}
col11Lookup = {}
col12Lookup = {}

#
# begin doSetup()
#
def doSetup():
    global dag
    global pubMed
    global evidenceLookup
    global taxonLookup
    global ecoLookupByEco, ecoLookupByEvidence
    global goPropertyLookup
    global goQualifierLookup
    global col11Lookup
    global col12Lookup
    global goRefLookup

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
    # col3 : go_qualifier
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
    # col11 : convert properties to RO id
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

        # MGI: -> MGI:MGI:
        value = r['value'].replace('MGI:', 'MGI:MGI:')

        # xxxx ; zzzz -> zzzz
        tokens = value.split(';')
        value = tokens[-1]
        value = value.strip()

        # EMAPA:xxx TS:xxx -> EMAPA:xxx
        if value.startswith('EMAPA:'):
                value = value.split(' ')[0]

        value = r['note'] + '(' + value + ')'

        if key not in col11Lookup:
            col11Lookup[key] = []
        col11Lookup[key].append(value)
    #print(col11Lookup)

    #
    # col12 
    # all annotations rows have creation/modification date
    # creation/modification dates
    #
    results = db.sql('''select distinct a._AnnotEvidence_key, a.cDate, a.mDate from gomarker2 a''', 'auto')
    for r in results:
        key = r['_AnnotEvidence_key']
        value = 'creation-date=' + r['cDate'] + '|' + 'modification-date=' + r['mDate']
        if key not in col12Lookup:
            col12Lookup[key] = []
        col12Lookup[key].append(value)
    #print(col12Lookup)

    #
    # col12 : use actual property values
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
                'has_participant', 'regulates_o_has_participant',
                'text'
                )
            ''', 'auto')
    for r in results:
        key = r['_AnnotEvidence_key']
        term = r['term']

        # to remove any non-ascii from "text" values
        value = ''.join(c for c in r['value'] if ord(c) >= 32)

        # text with "|" -> space
        # text with "," -> space
        value = value.replace('|', ' ')
        value = value.replace(',', ' ')

        if term in ('noctua-model-id', 'model-state'):
                value = term + '=' + value

        # "comment" is part of Dustin's initial noctua load; then it can be removed
        elif term in ('text'):
                value = 'comment=' + value
        elif term in ('has_participant', 'regulates_o_has_participant'):
                value = 'comment=' + term + '(' + value + ')'

        if key not in col12Lookup:
            col12Lookup[key] = []
        col12Lookup[key].append(value)
    #print(col12Lookup)

    #
    # col12 
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
        value = 'contributor-id=' + r['orcid']
        if key not in col12Lookup:
            col12Lookup[key] = []
        col12Lookup[key].append(value)
    #print(col12Lookup)

    results = db.sql('select _Object_key, accID from ACC_Accession where _LogicalDB_key = 185', 'auto')
    for r in results:
        key = r['_Object_key']
        value = r['accID']
        goRefLookup[key] = value
    #print(goRefLookup)

#
# end doSetup()
#

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

# begin doFinish()
#
def doFinish():

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

        #! 1  DB_Object_ID
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
                        reportRow = addReportRow(reportRow, r)
                except:
                        reportRow = MGIPREFIX + TAB + str(r['markerID']) + TAB
                        reportRow = addReportRow(reportRow, r)
        else:
            reportRow = MGIPREFIX + ':' + str(r['markerID']) + TAB
            reportRow = addReportRow(reportRow, r)

        fp.write(reportRow)

#
# end doFinish()
#

# begin addReportRow(reportRow, r)
#
def addReportRow(reportRow, r):

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
        #! 2  Negation
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
        #
        #! 3  Relation
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

        #! 4  Ontology_Class_ID
        reportRow = reportRow + r['termID'] + TAB

        #! 5  Reference
        references = []
        references.append(MGIPREFIX + ':' + r['refID'])
        if r['_Refs_key'] in pubMed:
            references.append('PMID:' + pubMed[r['_Refs_key']])
        else:
            if r['_Refs_key'] in goRefLookup:
                references.append(goRefLookup[r['_Refs_key']])
        reportRow = reportRow + '|'.join(references) + TAB

        #! 6  Evidence_type
        if key in evidenceLookup:
            reportRow = reportRow + evidenceLookup[key][0]
        elif r['evidenceCode'] in ecoLookupByEvidence:
            reportRow = reportRow + ecoLookupByEvidence[r['evidenceCode']]
        else:
            reportRow = reportRow + 'NOT FOUND'
        reportRow = reportRow + TAB

        #! 7  With_or_From
        inferredFrom = mgi_utils.prvalue(r['inferredFrom']).replace('MGI:', 'MGI:MGI:')
        reportRow = reportRow + mgi_utils.prvalue(inferredFrom) + TAB

        #! 8  Interacting_taxon_ID
        if key in taxonLookup:
            reportRow = reportRow + taxonLookup[key][0]
        reportRow = reportRow + TAB

        #! 9  Date
        reportRow = reportRow + str(r['mDate']) + TAB

        #! 10 Assigned_by

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

        #! 11 Annotation_Extensions
        properties = ''
        if key in col11Lookup:
            properties = ','.join(col11Lookup[key])
        reportRow = reportRow + properties + TAB

        #! 12 Annotation_Properties
        properties = ''
        if key in col12Lookup:
            properties = '|'.join(col12Lookup[key])
        reportRow = reportRow + properties + CRT

        return reportRow

#
# end addReportRow()
#

#
# main
#

db.useOneConnection(1)
doSetup()
doIsoform()

fp = reportlib.init('mgi2', fileExt = '.gpad', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp.write('!gpa-version: 2.0\n') 
fp.write('!generated-by: MGI\n')
fp.write('!date-generated: %s\n' % (mgi_utils.date("%Y-%m-%d")))
fp.write('!\n')
fp.write('!1  DB_Object_ID            ::= ID MGI or PR\n')
fp.write('!2  Negation                ::= "NOT"\n')
fp.write('!3  Relation                ::= OBO_ID\n')
fp.write('!4  Ontology_Class_ID       ::= OBO_ID\n')
fp.write('!5  Reference               ::= [ID] ("|" ID)*\n')
fp.write('!6  Evidence_type           ::= OBO_ID\n')
fp.write('!7  With_or_From            ::= [ID] ("|" | "," ID)*\n')
fp.write('!8  Interacting_taxon_ID    ::= NCBITaxon:[Taxon_ID]\n')
fp.write('!9  Date                    ::= YYYY-MM-DD\n')
fp.write('!10 Assigned_by             ::= Prefix\n')
fp.write('!11 Annotation_Extensions   ::= [Extension_Conj] ("|" Extension_Conj)*\n')
fp.write('!12 Annotation_Properties   ::= [Property_Value_Pair] ("|" Property_Value_Pair)*\n')
fp.write('!\n')

doFinish()

# append GOA annotations, if exists : see goload/goamouse
try:
        goafile = open(os.environ['GOAGPAD2MGI'], 'r')
        for line in goafile.readlines():
                fp.write(line)
        goafile.close()
except:
        pass

reportlib.finish_nonps(fp)

db.useOneConnection(0)
