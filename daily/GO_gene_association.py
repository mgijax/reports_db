#!/usr/local/bin/python

'''
#
# GO_gene_association.py 01/28/2002
#
# Report:
#       Tab-delimited file of all MGI GO/Marker associations for Stanford
#
# Usage:
#       GO_gene_association.py
#
# Used by:
#   Stanford - central GO site
#   also made publically available on MGI FTP site
#
# Output format:
#
# The GO format has the following columns:
#
#   1.  Database designation (MGI)
#   2.  MGI Marker ID
#   3.  Symbol
#   4.  NOT
#   5.  GO id
#   6.  MGI ID of Reference (in format MGI:MGI:####) (double MGI: necessary)
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
#   17. Isorform/Protein Data/MGI ID (Depending on data)
#
# History:
#
# ANY ERRORS WITH THE ANNOTATION RELATIONSHIP FILE:
#	MAKE SURE THE MAC EXCEL SPREADSHEET IS SAVED AS "WINDOWS FORMATED TEXT"
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
import string
import re
import mgi_utils
import reportlib

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
	db.setAutoTranslate(False)
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db

DBABBREV = 'MGI'
SPECIES = 'taxon:10090'
UNIPROTKB = 'uniprotload'
assignedByList = ['RGD', 'GOC', 'RefGenome']
annotRelationshipFileName = os.environ['GOANNOT_RELATIONSHIP']
annotRelationshipFileError = os.environ['GOANNOT_ERROR']

TAB = reportlib.TAB
CRT = reportlib.CRT

# see doSetup()
dag = {}
syns = {}
pubMed = {}

# see doIsoform()
isoformsProtein = {}
isoformsProtein2 = {}

# see doProtein()
proteins = {}
proteinsGene = {}

# see doCol16
# list of column 16 object/evidence/printable format
col16PrintLookup = {}
col16Error = ''

#
# begin doSetup()
#
def doSetup():
    global dag
    global syns
    global pubMed

    #
    # retrieve all dag abbrevations for each term
    #
    results = db.sql('select distinct _Object_key, rtrim(dagAbbrev) as dagAbbrev from DAG_Node_View where _Vocab_key = 4', 'auto')
    for r in results:
        dag[r['_Object_key']] = r['dagAbbrev']

    #
    # retrieve data set to process
    #
    #   and m.symbol = 'Slc14a2'
    #   and ta.accID in ('GO:0045944', 'GO:0006355', 'GO:0000122')
    #   and m.symbol = 'Slc14a2'
    #	and ta.accID in ('GO:0015204','GO:0002687','GO:0008594')
    #   and m.symbol = 'Mbd2'
    #
    db.sql('''select a._Term_key, t.term, ta.accID as termID, q.synonym as qualifier, a._Object_key, 
    	    e._AnnotEvidence_key, e.inferredFrom, e.modification_date, e._EvidenceTerm_key, 
    	    e._Refs_key, e._ModifiedBy_key, 
    	    m.symbol, m.name, lower(mt.name) as markerType
        into #gomarker1 
        from VOC_Annot a, 
	     ACC_Accession ta, 
	     VOC_Term t, 
	     VOC_Evidence e, 
	     MRK_Marker m, 
	     MRK_Types mt, 
	     MGI_Synonym q 
        where a._AnnotType_key = 1000 
        and a._Annot_key = e._Annot_key 
        and a._Object_key = m._Marker_key 
        and m._Marker_Type_key = 1 
        and a._Term_key = t._Term_key 
        and a._Term_key = ta._Object_key 
        and ta._MGIType_key = 13 
        and ta.preferred = 1 
        and m._Marker_Type_key = mt._Marker_Type_key 
        and a._Qualifier_key = q._Object_key 
        and q._SynonymType_key = 1023
        ''', None)
    db.sql('create index gomarker1_idx1 on #gomarker1(_Object_key)', None)
    db.sql('create index gomarker1_idx2 on #gomarker1(_EvidenceTerm_key)', None)
    db.sql('create index gomarker1_idx3 on #gomarker1(_Refs_key)', None)
    db.sql('create index gomarker1_idx4 on #gomarker1(_ModifiedBy_key)', None)
    db.sql('create index gomarker1_idx5 on #gomarker1(_AnnotEvidence_key)', None)

    #
    # retrieve synonyms for markers in data set
    #
    results = db.sql('''select distinct g._Object_key, s.synonym 
        from #gomarker1 g, MGI_Synonym s, MGI_SynonymType st 
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
    db.sql('''select g._Refs_key, g._Term_key, g.termID, g.qualifier, g.inferredFrom, 
    	    g._Object_key, g._AnnotEvidence_key, g._EvidenceTerm_key, g.symbol, g.name, g.markerType, 
    	    convert(varchar(10), g.modification_date, 112) as mDate, 
    	    ma.accID as markerID, 
    	    b.accID as refID, 
    	    rtrim(t.abbreviation) as eCode, 
    	    u.login as assignedBy
        into #gomarker2 
        from #gomarker1 g, ACC_Accession ma, ACC_Accession b, VOC_Term t, MGI_User u 
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
    db.sql('create index gomarker2_idx1 on #gomarker2(symbol)', None)
    db.sql('create index gomarker2_idx2 on #gomarker2(_Refs_key)', None)
    db.sql('create index gomarker2_idx3 on #gomarker2(_AnnotEvidence_key)', None)
    db.sql('create index gomarker2_idx4 on #gomarker2(_Object_key)', None)
    db.sql('create index gomarker2_idx5 on #gomarker2(_EvidenceTerm_key)', None)
    
    #
    # resolve PubMed IDs for References
    #
    results = db.sql('''select distinct r._Refs_key, a.accID from #gomarker2 r, ACC_Accession a 
        where r._Refs_key = a._Object_key 
        and a._MGIType_key = 1 
        and a._LogicalDB_key = 29''', 'auto')
    for r in results:
        key = r['_Refs_key']
        value = r['accID']
        pubMed[key] = value
#
# end doSetup()
#

#
# begin doCol16()
#
def doCol16():

    global col16PrintLookup
    global col16Error

    #
    # create error file
    #
    try:
	col16Error = open(annotRelationshipFileError, 'w')
    except:
	print 'cannot open error file %s : no errors can be reported\n' % (annotRelationshipFileError)

    #
    # TR11112/use hand-coded annotation relationship file to determine column 16
    #
    annotRelationshipLookup = {}
    annotRelationshipFile = open(annotRelationshipFileName, 'r')
    for line in annotRelationshipFile.readlines():
        tokens = line[:-1].split(TAB)
        key = str(tokens[0]) + ':' + tokens[2]
        value = tokens[3].strip()
        if not annotRelationshipLookup.has_key(key):
	    annotRelationshipLookup[key] = []
        annotRelationshipLookup[key].append(value)
    annotRelationshipFile.close()
    #print annotRelationshipLookup

    cellLineSelect = re.compile(r'(CL:[0-9]+)', re.I)
    targetSelect = re.compile(r'(MGI:[0-9]+|UniProtKB:+[0-9])', re.I)
    anatomySelect = re.compile(r'(MA:[0-9]+|EMAP:[0-9]+)', re.I)

    # list of object/evidence/stanza rows
    col16EvidenceLookup = {}
    # list of property/value for one object/evidence row
    col16PropertyLookup = {}

    #
    # list of object/evidence/stanza with 'anatomy', 'cell type', 'target' properties
    # key = marker key:annotation/evidence key
    # values = all results (includes  property, value)
    #
    # exclude annotations where evidence = ISO (3251466)
    #
    results = db.sql('''select r.symbol, r._Object_key, r._AnnotEvidence_key, r._Term_key, r.termID, t.term as property, p.value, p.stanza
	    from #gomarker2 r, VOC_Evidence_Property p, VOC_Term t
	    where r._EvidenceTerm_key not in (3251466)
            and r._AnnotEvidence_key = p._AnnotEvidence_key
	    and p._PropertyTerm_key in (6481774, 6481773, 6481777)
	    and p._PropertyTerm_key = t._Term_key
	    order by r.symbol, r._Object_key, r._AnnotEvidence_key, p.stanza, p.sequenceNum, property''', 'auto')
    for r in results:
        key = str(r['_Object_key']) + ':' + str(r['_AnnotEvidence_key']) + ':' + str(r['stanza'])
        value = r
        if not col16EvidenceLookup.has_key(key):
	    col16EvidenceLookup[key] = []
        col16EvidenceLookup[key].append(value)
    #print col16EvidenceLookup

    #
    # iterate thru list of each object/evidence
    #
    for r in col16EvidenceLookup:
        col16PropertyLookup.clear()

        #
        # iterate thru list of property/value for a specific object/evidence/stanza
        # create the "property relationship(value)" text
        # for example:
        #   go id: GO:0000122
        #	property = target
        #   value = Ihh ; MGI:96533
        # text = "has_regulation_target(MGI:MGI:96533)"
        #
        for v in col16EvidenceLookup[r]:
       
	   isError = 0

           # determine property relationship
           aKey = v['termID'] + ':' + v['property']
           if annotRelationshipLookup.has_key(aKey):
               startPart = annotRelationshipLookup[aKey][0]

	   # cellular component default
	   elif dag[v['_Term_key']] == 'C':
               startPart = 'part_of'

	   # molecular function and property in ('cell type', 'anatomy')
	   elif dag[v['_Term_key']] == 'F' and v['property'] in ['cell type', 'anatomy']:
               startPart = 'occurs_in'

	   # molecular function and property not in ('cell type', 'anatomy')
	   # biological process
           else:
	       isError = 1
	       try:
		   col16Error.write(v['termID'] + '\t' + v['symbol'] + '\t' + v['property'] + '\t' + v['value'] + '\n')
               except:
		   pass

	   # if errors were found...skip this property/value
	   if isError == 1:
	       continue

           p = v['property']

           if not col16PropertyLookup.has_key(p):
	       col16PropertyLookup[p] = []

           # extract anatomy information
           if p == 'anatomy':
               for pv in anatomySelect.findall(v['value']):
	          col16PropertyLookup[p].append(startPart + '(' + pv + ')')

           # extract cell line information
           if p == 'cell type':
               for pv in cellLineSelect.findall(v['value']):
	           if pv not in ('CL:0000001', 'CL:0000002'):
	               col16PropertyLookup[p].append(startPart + '(' + pv + ')')

           # extract target information
           if p == 'target':
               for pv in targetSelect.findall(v['value']):
	          pv = pv.replace('MGI', 'MGI:MGI')
	          col16PropertyLookup[p].append(startPart + '(' + pv + ')')

        #
        # all properties collected for given object/evidence
        # construct column 16 print value
        #

        # make sure there is at least one value per property
        if not col16PropertyLookup.has_key('anatomy'):
	    col16PropertyLookup['anatomy'] = []
	    col16PropertyLookup['anatomy'].append('')
        if not col16PropertyLookup.has_key('cell type'):
	    col16PropertyLookup['cell type'] = []
	    col16PropertyLookup['cell type'].append('')
        if not col16PropertyLookup.has_key('target'):
	    col16PropertyLookup['target'] = []
	    col16PropertyLookup['target'].append('')

        # prepare col16PrintLookup
	# remove 'stanza' and put all stanza's 
	# for given object/evidence into one print statement
	# this assumes that stanza's will not exceed '9'
	ckey = r[:-2]
        if not col16PrintLookup.has_key(ckey):
            col16PrintLookup[ckey] = []

        #
        # iterate thru each property
        # concatenate one value per property with "," separator
	#
        # for example:
        #  anatomy 1,cell type 1,target 1
        #  anatomy 1,cell type 2,target 2
        #  anatomy 1,cell type 2,target 1
        #  anatomy 1,cell type 2,target 2
        #
        for ap in col16PropertyLookup['anatomy']:
            for cp in col16PropertyLookup['cell type']:
                for tp in col16PropertyLookup['target']:
		    temp = []
		    if len(ap) > 0:
		        temp.append(ap)
		    if len(cp) > 0:
		        temp.append(cp)
		    if len(tp) > 0:
		        temp.append(tp)
		    col16PrintLookup[ckey].append(string.join(temp, ','))
    #print col16PrintLookup

    col16Error.close()

## end doCol16()

#
# begin doIsoform()
#
def doIsoform():
    global isoformsProtein
    global isoformsProtein2

    #
    # isoformsProtein hash
    # select all "gene_product:" properties
    # key = marker key:annotation/evidence key ()
    # values = [UniProt:XXXXX,UniProt:XXXXX]
    #

    r1 = re.compile(r'([^\s\\\n]*)', re.I)
    isoformPattern1 = re.compile(r'UniProtKB:', re.I)
    isoformPattern2 = re.compile(r'protein_id', re.I)
    isoformPattern3 = re.compile(r'NCBI:NP_', re.I)
    isoformPattern4 = re.compile(r'NCBI:XP_', re.I)
    isoformPattern5 = re.compile(r'PR:', re.I)

    isoformsProtein = {}
    isoformsProtein2 = {}
    results = db.sql('''select r._Object_key, r._AnnotEvidence_key, p.value
	from #gomarker2 r, VOC_Evidence_Property p
	where r._AnnotEvidence_key = p._AnnotEvidence_key
	and p._PropertyTerm_key = 6481775
	order by r._Object_key, r._AnnotEvidence_key, p.stanza, p.sequenceNum''', 'auto')

    for r in results:

        key = str(r['_Object_key']) + ':' + str(r['_AnnotEvidence_key'])
        value = r['value']

        for a in r1.findall(value):
	    for b in a.split('|'):
	        b = b.strip()

                # Only certain patterns actually count, they are listed above.
                if isoformPattern1.match(b) != None or \
	           isoformPattern2.match(b) != None or \
                   isoformPattern3.match(b) != None or \
	           isoformPattern4.match(b) != None or \
	           isoformPattern5.match(b) != None:

	           # TR10652
	           # convert 'NCBI:' to 'RefSeq'
	           # once TR10044 is implemented and all NCBI are migrated to RefSeq,
	           # this will no longer be necessary

	           b = b.replace('NCBI:', 'RefSeq:')

                   if not isoformsProtein.has_key(key):
	               isoformsProtein[key] = []
                   isoformsProtein[key].append(b)

	           #
	           # TR11060
	           # if UniProtKB and contains "-"
	           #
                   if isoformPattern1.match(b) != None and string.find(b, '-') >= 0:
                       if not isoformsProtein2.has_key(key):
	                   isoformsProtein2[key] = []
                       isoformsProtein2[key].append(b)
#
# end doIsoform()
#

#
# begin doProtein()
#
def doProtein():
    global proteins
    global proteinsGene

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
    # proteinsGene
    #    representative transcript (615420) for marker type "gene" (1)
    #    marker symbol like 'mir%' (microRNAs)
    #    9 = GenBank
    #    27 = RefSeq
    #
    # example of counts:
    #  3935 go-uniprot.txt   proteins
    #  4519 go-npxp.txt      proteins
    #     0 go-genbank.txt   proteinsGene
    # 16946 go-all.txt       proteinsGene
    #

    results = db.sql('''
        select distinct r.symbol, mc._Marker_key, mc.accID as seqID, mc._LogicalDB_key, mc._Qualifier_key
        from #gomarker2 r, SEQ_Marker_Cache mc 
        where r._Object_key = mc._Marker_key 
        and mc._Marker_Type_key = 1 
        and mc._Qualifier_key = 615421 
        union 
        select distinct r.symbol, mc._Marker_key, mc.accID as seqID, mc._LogicalDB_key, mc._Qualifier_key
        from #gomarker2 r, SEQ_Marker_Cache mc, MRK_MCV_Cache mcv
        where r._Object_key = mc._Marker_key 
        and mc._Marker_Type_key = 1
        and mc._Qualifier_key = 615420
        and mc._Marker_key = mcv._Marker_key
        and mcv.term = 'miRNA Gene'
        ''', 'auto')
    
    proteins = {}
    proteinsGene = {}
    
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

        # GenBank
        elif logicalDB in [9]:
	    #print 'ldb 9: ', str(symbol), str(seqID)
            proteinsGene[key] = 'EMBL:' + seqID

        else:
	    #print 'all else: ', str(symbol), str(seqID)
            proteinsGene[key] = 'RefSeq:' + seqID
#
# end doProtein()
#

#
# begin doFinish()
#
def doFinish():

    #
    # process results
    #
    results = db.sql('select * from #gomarker2 order by symbol, termID', 'auto')
    for r in results:

        reportRow = ''    

        if dag.has_key(r['_Term_key']):

	    objectKey = str(r['_Object_key']) + ':' + str(r['_AnnotEvidence_key'])

	    # columns 1-5
            reportRow = DBABBREV + TAB
            reportRow = reportRow + str(r['markerID']) + TAB
            reportRow = reportRow + r['symbol'] + TAB
            reportRow = reportRow + string.strip(r['qualifier']) + TAB
            reportRow = reportRow + r['termID'] + TAB

            # column 6; reference
            referenceID = DBABBREV + ':' + r['refID']

            if pubMed.has_key(r['_Refs_key']):
                referenceID = referenceID + '|PMID:' + pubMed[r['_Refs_key']]

            reportRow = reportRow + referenceID + TAB

	    # column 7
            reportRow = reportRow + r['eCode'] + TAB

	    # column 8
            inferredFrom = re.sub('MGI:','MGI:MGI:',mgi_utils.prvalue(r['inferredFrom']))
            reportRow = reportRow + inferredFrom + TAB

	    # column 9-10
            reportRow = reportRow + dag[r['_Term_key']] + TAB
            reportRow = reportRow + r['name'] + TAB

	    # column 11
            if syns.has_key(r['_Object_key']):
                syn_string = string.join(syns[r['_Object_key']], '|')
                reportRow = reportRow + syn_string + TAB
            else:
                reportRow = reportRow + TAB

            # column 12 is populated in a special way
            # If there is a isoform or a protein, print out protein
            # otherwise print out the marker type.

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
            elif string.find(r['assignedBy'], 'GOA_') >= 0:
                assignedBy = re.sub('GOA_', '', r['assignedBy'])
                reportRow = reportRow + assignedBy + TAB

	    # check list
            elif r['assignedBy'] in assignedByList:
                reportRow = reportRow + r['assignedBy'] + TAB

	    # else use default
            else:
                reportRow = reportRow + DBABBREV + TAB

	    #
	    # column 16
	    # contains property/value information
	    #

            if col16PrintLookup.has_key(objectKey):
                reportRow = reportRow + string.join(col16PrintLookup[objectKey], '|') + TAB
            else:                    
                reportRow = reportRow + TAB

            # column 17
	    #
	    # if isoformProtein = true
	    #    then use isoformsProtein
	    #
	    # else if isoformProtein = false, proteins = true
	    #    then use proteins
	    #
	    # else if isoformProtein = false, proteins = false, proteinsGene = true
	    #    then use proteinsGene
	    #
	    # else if isoformProtein = false, proteins = false, proteinsGene = false
	    #    then leave blank
	    #

            if isoformsProtein.has_key(objectKey):
	        #print 'isoform:  ', str(r['symbol']), string.join(isoformsProtein[objectKey])
                reportRow = reportRow + string.join(isoformsProtein[objectKey], '|') + CRT
            else:
	        row = ''
                if proteins.has_key(r['_Object_key']):
	            #print 'protein:  ', str(r['symbol']), str(symbol), str(proteins[r['_Object_key']])
                    row = str(proteins[r['_Object_key']])
    
                elif proteinsGene.has_key(r['_Object_key']):
	            #print 'protein/gene:  ', str(r['symbol']), str(proteinsGene[r['_Object_key']])
                    row = str(proteinsGene[r['_Object_key']])

	        # intentionally commented out
                #else:
	        #    print 'blank:  ', str(r['symbol'])

	        reportRow = reportRow + row + CRT

            fp.write(reportRow)

	    #
	    # TR11060
	    # subset of UniProtKB:xxxx-?? only
	    #
	    if isoformsProtein2.has_key(objectKey):
                fp2.write(reportRow)
 
#
# end doFinish()
#

#
# main
#

db.useOneConnection(1)

fp = reportlib.init('gene_association', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp2 = reportlib.init('gene_association_pro', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# Header information
# This file gets submitted to CVS so don't muck with this header
# Update the Revision number every time you modify this code 
# (like adding a new column)
#

fp.write('!gaf-version: 2.0\n')
fp.write('!software version: $Revision$\n')
fp.write('!date: %s $\n' % (mgi_utils.date("%m/%d/%Y")))
fp.write('!\n')
fp.write('! from Mouse Genome Database (MGD) & Gene Expression Database (GXD)\n')
fp.write('!\n')

doSetup()
doCol16()
doIsoform()
doProtein()
doFinish()

#
# append GOA annotations, if they exist
#

try:

    goaFile = open(os.environ['GOAMGI'], 'r')

    for line in goaFile.readlines():
        fp.write(line.rstrip()+TAB+TAB+CRT)
    goaFile.close()

except:

    pass

reportlib.finish_nonps(fp)
db.useOneConnection(0)
