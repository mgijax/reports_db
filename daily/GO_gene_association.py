#!/usr/local/bin/python

'''
#
# GO_gene_association.py
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
# IMPORTANT THINGS TO KNOW:
#
#    gaf/col11 (annotation extension) and gpad/col16 should be equal.
#    except :  ISO annotations are excluded from gaf/col11 (blank)
#    see "lib_py_report/go_annot_extensions.py" for the list of
#	excluded Properties and excluded Evidence
#
#    gpad/col16 : will *never* contains > 1 stanza, and will always use the "," delimiter
#
# History:
#
# lec	07/14/2016
#	- TR12349/12345/GO_Noctua/GPAD added
#
# kstone 09/14/2015
#	- TR12070 Refactored col16 and col17 logic into 'go_annot_extensions' and 'go_isoforms' modules
#		Removed extra protein values for col17. Use only isoform 'gene product' GO properties now
#
# lec	04/22/2015
#	- TR11932/added "GO_" (for GO_Central) to set proper assignedBy value
#
# sc	11/06/2014
#	- TR11772/Modification of GAF file references
#
# lec   10/24/2014
#       - TR11750/postres complient
#
# lec	07/02/2014
#	- TR11693/only include official/interum markers
#
# lec	07/01/2014
#	- TR11710/doGAFCol16
#
# lec	01/13/2014
#	- TR11570/fix GO qualifier; use VOC_Term
#
# lec	01/02/2014
#	- TR11518/add "new" doGAFCol16()
#	- all mgi-properties shoudl now be coverted to go-properties
#
# lec	08/27/2013
#	- TR11459/GOANNOT_RELATIONSHIP/fix
#
# lec	05/28/2013
#	- TR11060/add all UniProtKB: and PR:
#
# lec	01/15/2013
#	- TR11112/use GOANNOTRELATIONSHIP to generate column 16
#
# lec	05/08/2012
#	- TR11060/add secondary file/subset/
#	  only contains annotations where column 17 (gene product) has UniProtKB:xxxx-??
#	  gene_association_pro.mgi
#
# lec	03/13/2012
# lec	03/08/2012
#	pubMed = {} was being called/selected twice
#	pubMed = change query to 'distinct'
#
# lec	06/20/2011
#   - TR10044/MGI_Notes --> VOC_Evidence_Property
#	this affects col16LookupByEvidence, isoformsProtein, column 12, 16, 17
#
# lec	03/30/2011
#   - TR10652/change 'NCBI:' to 'RefSeq:'
#
# lec	03/21/2011
#   - TR9962/add 'RefGenome' to column 15
#
# lec	03/15/2011
#   - TR10633/allow PRO ids (PR:) in column 17
#
# lec
#   - TR6839/marker types
#   - marker type 11 (microRNA) moved to marker type 1 (gene)
#   - changes to proteins/proteinsGene hash/column 17
#
# lec   06/17/2010
#   - check logicalDB for proteins and proteinsGene hash
#     and fix prefix name for Vega, Ensembl
#
# lec	06/10/2010
#   - cleanup up cell ontology, isoform protein and protein hashes
#   - added TR9901/date history (see below)
#
# lec	06/03/2010
#   - TAB not being written in between column 15/column 16
#   - column 16/17: replace MRK_Marker with go/marker temp table
#
# lec	04/29/2010
#   - TR9777/"swissload" login name changed to "uniprotload"
#
# lec	03/02/2010
#   - TR10035; added RGD check for column 15
#
# mhall 02/02/2010
#   - TR 9901; column 12, 16, 17
#
# lec   07/24/2008
#   - TR 9134; change UniProt to UniProtKB
#
# lec   05/07/2008
#   - TR 8997; lowercase the marker types
#
# lec   01/25/2007
#   - TR 8122; don't convert inferredFrom delimiters; use as is in the database
#
# lec   09/14/2006
#   - TR 7904; append GOA annotations   
#
# lec   10/19/2005
#   - added PMID, TR 7173
#
# lec   10/04/2005
#   - TR 5188; GO Qualifier
#
# lec   10/03/2005
#   - replace SWALL with UniProt
#
# lec   03/10/2004
#       - only include Markers of type Gene
#
# lec   02/11/2003
#   - TR 4511; add new column "assigned by"
#
# lec   04/02/2002
#   - TR 3527; some terms may not have DAGs; this shouldn't happen
#   but we don't want the report to bomb if it does.
#
# lec   01/28/2002
#   - new - revision 2.0
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
db.setAutoTranslate(False)
db.setAutoTranslateBE(False)

MGIPREFIX = 'MGI'
SPECIES = 'taxon:10090'
UNIPROTKB = 'uniprotload'
assignedByList = ['RGD', 'GOC', 'GO_Central', 'GO_Noctua']

TAB = reportlib.TAB
CRT = reportlib.CRT

# see doSetup()
dag = {}
syns = {}
pubMed = {}

# mapping of load reference keys to GO_REF IDs/per TR11772
# JNUM          _Refs_key      Desired addition
# J:155856	156949         Rat2mouseISO     GO_REF:0000096 ???
# J:161428	162524         GO_REF:0000033
# J:164563	165659         human2mouseISO	 GO_REF:0000096 ???
# J:60000	61933          GO_REF:0000004
# J:72245	73197          GO_REF:0000003
# J:72247	73199          GO_REF:0000002
# J:73065	74017          GO_REF:0000008
# J:73796	74750          GO_REF:0000015
goRefDict = {156949:'GO_REF:0000096',
             162524:'GO_REF:0000033',
             165659:'GO_REF:0000096',
             61933:'GO_REF:0000004',
             73197:'GO_REF:0000003',
             73199:'GO_REF:0000002',
             74017:'GO_REF:0000008',
             74750:'GO_REF:0000015',
             }

#
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
dagQualifier = {'C':'part_of', 'P':'involved_in', 'F':'enables'}
taxonLookup = {}
ecoLookupByEco = {}
ecoLookupByEvidence = {}
evidenceLookup = {}
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
    global gpadCol11Lookup
    global gpadCol12Lookup

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
        if not syns.has_key(key):
            syns[key] = []
        syns[key].append(value)

    #
    # resolve foreign keys and store in "results" table
    #
    db.sql('''select distinct g._Refs_key, g._Term_key, g.termID, g.qualifier, g.inferredFrom, 
    	    g._Object_key, g._AnnotEvidence_key, g._EvidenceTerm_key, g.symbol, g.name, g.markerType, 
	    to_char(g.modification_date, 'YYYYMMDD') as mDate,
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
    #print ecoLookupByEvidence

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
    #print evidenceLookup

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
    #print taxonLookup 

    #
    # gpadCol11 : noctua (_CreatedBy_key = 1559)
    #	exclude older terms (sequenceNum 1-9, 90,91,92, 93)
    #
    # note that noctua-generated properties will *always* have one stanza
    #
    results = db.sql('''select distinct a._AnnotEvidence_key, t.term, p.value
            from gomarker2 a,
                 VOC_Evidence_Property p,  
                 VOC_Term t
	    where a._CreatedBy_key = 1559
	    and a._AnnotEvidence_key = p._AnnotEvidence_key
	    and p._PropertyTerm_key = t._Term_key
	    and t.term not in (
    		'evidence', 'anatomy', 'cell type', 'gene product', 'modification', 'target', 
    		'external ref', 'text', 'dual-taxon ID',
    		'lego-model-id', 'contributor', 'individual', 'go_qualifier'
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
    # gpadCol12 : noctua (_CreatedBy_key = 1559) 
    # exclude : occurs_in, part_of, go_qualifier, evidence
    #
    results = db.sql('''select distinct a._AnnotEvidence_key, t.term, p.value
            from gomarker2 a,
                 VOC_Evidence_Property p,  
                 VOC_Term t
	    where a._CreatedBy_key = 1559
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

	if not gafCol16Lookup.has_key(objectKey):
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

    #print gafCol16Lookup

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
    #    132 = VEGA Protein
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
	    #print 'ldb 13/41: ', str(symbol), str(seqID)
            proteins[key] = 'UniProtKB:' + seqID

        # RefSeq
        elif logicalDB in [27] and qualifier == 615421:
	    #print 'np/xp: ', str(symbol), str(seqID)
            proteins[key] = 'RefSeq:' + seqID

        # Vega
        elif logicalDB in [132]:
            proteins[key] = 'VEGA:' + r['seqID']

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
    # Output format:
    #
    # The GO format has the following columns:
    #
    #   1.  Database designation (MGI)
    #   2.  MGI Marker ID (MGI:xxxx)
    #   3.  Symbol
    #   4.  Qualifier
    #   5.  GO id
    #   6.  MGI ID of Reference (MGI:MGI:xxxx|PMID:xxxx)
    #   7.  Evidence abbreviation
    #   8.  Inferred From
    #   9.  GO DAG Abbreviation (F, P, C)
    #   10. Gene name
    #   11. Gene synonym(s) - list of |-delimited synonyms
    #   12. Marker Type or Protein (gene)
    #   13. Species (taxon:10090)
    #   14. Modification Date (YYYYMMDD)
    #   15. Assigned By
    #   16. Properites/Values (occurs_in, part_of, etc.)
    #   17. Isorform
    #

    #
    # process results
    #
    results = db.sql('select * from gomarker2 order by symbol, termID', 'auto')

    for r in results:

        reportRow = ''    

        if r['_Term_key'] not in dag:
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
        if pubMed.has_key(r['_Refs_key']):
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
        if syns.has_key(r['_Object_key']):
            syn_string = ','.join(syns[r['_Object_key']])
            reportRow = reportRow + syn_string + TAB
        else:
            reportRow = reportRow + TAB

        # column 12
	# if marker is associated with an isoform (via go/annotation)
	# or marker is associated with a protein (via marker/sequence cache)
        # 	print 'protein' 
        # else, print marker type (ex. 'gene')

        if isoformsProtein.has_key(objectKey) or proteins.has_key(r['_Object_key']):
            reportRow = reportRow + 'protein' + TAB
        else:
            reportRow = reportRow + r['markerType'] + TAB
                
	# column 13
        reportRow = reportRow + SPECIES + TAB

	# column 14
        reportRow = reportRow + str(r['mDate']) + TAB

	# column 15; assigned by
        if r['assignedBy'] == UNIPROTKB:
            reportRow = reportRow + 'UniProtKB' + TAB

	# remove "GOA_"; for example:  "GOA_IntAct" ==> "IntAct"
	elif r['assignedBy'].find('GOA_') >= 0:
            assignedBy = r['assignedBy'].replace('GOA_', '')
            reportRow = reportRow + assignedBy + TAB

	# check list
        elif r['assignedBy'] in assignedByList:
            reportRow = reportRow + r['assignedBy'] + TAB

	# else use default
        else:
            reportRow = reportRow + MGIPREFIX + TAB

	#
	# column 16
	# contains property/value information
	# see lib_py_report/go_annot_extensions.py for list of excluded properties
	properties = ''
        if gafCol16Lookup.has_key(objectKey):
	    properties = ''.join(gafCol16Lookup[objectKey])
        reportRow = reportRow + properties + TAB

        # column 17
	# if isoformProtein = true
	#    then use isoformsProtein
	isoforms = ''
        if isoformsProtein.has_key(objectKey):
	    isoforms = '|'.join(isoformsProtein[objectKey])
        reportRow = reportRow + isoforms + CRT

        fp.write(reportRow)

	#
	# TR11060
	# subset of UniProtKB:xxxx-?? only
	#
	if forPROC.has_key(objectKey):
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
    #   12. Annotation Properties   (optional) properties if creator = 'GO_Noctua' (1559)

    #
    # process results
    #
    results = db.sql('select * from gomarker2 order by symbol, termID', 'auto')

    for r in results:

        reportRow = ''    

        if r['_Term_key'] not in dag:
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

        if isoformsProtein.has_key(objectKey):
	    for i in isoformsProtein[objectKey]:
	        prefixPart, numericPart = i.split(':')
                reportRow = prefixPart + TAB + numericPart + TAB
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
	#   always add new C/P/F qualifier
	#   if any additional qualifiers (contributes_to, colocalizes_with, etc.) exists, attach them too
	qualifier = dagQualifier[dag[r['_Term_key']]]
	if r['qualifier'] != None:
	    qualifier = qualifier + '|' + r['qualifier'].strip()
        reportRow = reportRow + qualifier + TAB

	#   4. GO ID
        reportRow = reportRow + r['termID'] + TAB

	#   5. DB:Reference(s)
	references = []
	references.append(MGIPREFIX + ':' + r['refID'])
        if pubMed.has_key(r['_Refs_key']):
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
        if r['assignedBy'] == UNIPROTKB:
            reportRow = reportRow + 'UniProtKB' + TAB

	# remove "GOA_"; for example:  "GOA_IntAct" ==> "IntAct"
	elif r['assignedBy'].find('GOA_') >= 0:
            assignedBy = r['assignedBy'].replace('GOA_', '')
            reportRow = reportRow + assignedBy + TAB

	# check list
        elif r['assignedBy'] in assignedByList:
            reportRow = reportRow + r['assignedBy'] + TAB

	# else use default
        else:
            reportRow = reportRow + MGIPREFIX + TAB

	#   11. Annotation Extension
	#   if GO_Nocuta, then use gpadCol11Lookup else use gafCol16Lookup
	properties = ''
	if key in gpadCol11Lookup:
            properties = ','.join(gpadCol11Lookup[key])
	elif gafCol16Lookup.has_key(objectKey):
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
# GAF
#

fp = reportlib.init('gene_association', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp2 = reportlib.init('gene_association_pro', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('!gaf-version: 2.1\n')
fp.write('!software version: $Revision$\n')
fp.write('!date: %s $\n' % (mgi_utils.date("%m/%d/%Y")))
fp.write('!\n')
fp.write('! from Mouse Genome Database (MGD) & Gene Expression Database (GXD)\n')
fp.write('!\n')
fp.write('!	1.  DB                       MGI\n')
fp.write('!	2.  DB Object ID             MGI:xxxx\n')
fp.write('!	3.  DB Object Symbol\n')
fp.write('!	4.  Qualifier\n')
fp.write('!	5.  GO ID                    GO:xxxx\n')
fp.write('!	6.  DB:Reference(s)          MGI:MGI:xxxx|PMID:xxxx\n')
fp.write('!	7.  Evidence Code            3-digit (not ECO:xxxx)\n')
fp.write('!	8.  With (or)From            optional\n')
fp.write('!	9.  Aspect (GO DAG Abbreviation (F, P, C))\n')
fp.write('!	10. DB Object Name           optional\n')
fp.write('!	11. DB Object Synonym(s)     optional\n')
fp.write('!	12. DB Object Type\n')
fp.write('!	13. Taxon                    taxon:10090\n')
fp.write('!	14. Date                     YYYYMMDD\n')
fp.write('!	15. Assigned By\n')
fp.write('!	16. Annotation Extension     same as GPAD/col 11\n')
fp.write('!	17. Gene Product Form ID     Isorform\n')
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
# GPAD
#

fp = reportlib.init('mgi', fileExt = '.gpa', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('!gpa-version: 1.1\n') 
fp.write('!\n')
fp.write('!date: %s $\n' % (mgi_utils.date("%m/%d/%Y")))
fp.write('!\n')
fp.write('! from Mouse Genome Database (MGD) & Gene Expression Database (GXD)\n')
fp.write('!\n')
fp.write('!	1.  DB                       MGI or PR\n')
fp.write('!	2.  DB Object ID             MGI:xxxx or xxxxx\n')
fp.write('!	3.  Qualifier                enables, involved_in, part_of\n')
fp.write('!	4.  GO ID                    GO:xxxx\n')
fp.write('!	5.  DB:Reference(s)          MGI:MGI:xxxx|PMID:xxxx\n')
fp.write('!	6.  Evidence Code            ECO:xxxx\n')
fp.write('!	7.  With (or)From            optional\n')
fp.write('!	8.  Interacting taxon ID     optional\n')
fp.write('!	9.  Date                     YYYYMMDD\n')
fp.write('!	10. Assigned By\n')
fp.write('!	11. Annotation Extension     optional same as GAF/col 16\n')
fp.write('!	12. Annotation Properties    optional\n')
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

reportlib.finish_nonps(fp)

db.useOneConnection(0)

