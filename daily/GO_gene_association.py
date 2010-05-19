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
#   16. Cell Ontology (occurs_in or part_of)
#   17. Isorform/Protein Data/MGI ID (Depending on data)
#
# History:
#
# lec	04/29/2010
#   - TR9777/"swissload" login name changed to "uniprotload"
#
# lec	03/02/2010
#   - TR10035; added RGD check for column 15
#   note that the isoform edits have no TR # attached
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
import db
import reportlib
import mgi_utils


DBABBREV = 'MGI'
SPECIES = 'taxon:10090'
UNIPROTKB = 'uniprotload'
assignedByList = ['RGD', 'GOC']

TAB = reportlib.TAB
CRT = reportlib.CRT

db.useOneConnection(1)

fp = reportlib.init('gene_association', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

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

#
# retrieve all dag abbrevations for each term
#
results = db.sql('select distinct _Object_key, dagAbbrev = rtrim(dagAbbrev) from DAG_Node_View where _Vocab_key = 4', 'auto')
dag = {}
for r in results:
    dag[r['_Object_key']] = r['dagAbbrev']

#
# retrieve data set to process
#
db.sql('select a._Term_key, t.term, termID = ta.accID, qualifier = q.synonym, a._Object_key, ' + \
    'e.inferredFrom, e.modification_date, e._EvidenceTerm_key, e._Refs_key, e._ModifiedBy_key, ' + \
    'm.symbol, m.name, markerType = lower(mt.name) ' + \
    'into #gomarker ' + \
    'from VOC_Annot a, ACC_Accession ta, VOC_Term t, VOC_Evidence e, MRK_Marker m, MRK_Types mt, MGI_Synonym q ' + \
    'where a._AnnotType_key = 1000 ' + \
    'and a._Annot_key = e._Annot_key ' + \
    'and a._Object_key = m._Marker_key ' + \
    'and m._Marker_Type_key = 1 ' + \
    'and a._Term_key = t._Term_key ' + \
    'and a._Term_key = ta._Object_key ' + \
    'and ta._MGIType_key = 13 ' + \
    'and ta.preferred = 1 ' + \
    'and m._Marker_Type_key = mt._Marker_Type_key ' + \
    'and a._Qualifier_key = q._Object_key ' + \
    'and q._SynonymType_key = 1023', None)
db.sql('create index idx1 on #gomarker(_Object_key)', None)
db.sql('create index idx2 on #gomarker(_EvidenceTerm_key)', None)
db.sql('create index idx3 on #gomarker(_Refs_key)', None)
db.sql('create index idx4 on #gomarker(_ModifiedBy_key)', None)

#
# retrieve synonyms for markers in data set
#
results = db.sql('select distinct g._Object_key, s.synonym ' + \
    'from #gomarker g, MGI_Synonym s, MGI_SynonymType st ' + \
    'where g._Object_key = s._Object_key ' + \
    'and s._MGIType_key = 2 ' + \
    'and s._SynonymType_key = st._SynonymType_key ' + \
    'and st.synonymType = "exact" ' + \
    'order by g._Object_key', 'auto')
syns = {}
for r in results:
    key = r['_Object_key']
    value = r['synonym']
    if not syns.has_key(key):
        syns[key] = []
    syns[key].append(value)

#
# resolve foreign keys
#
db.sql('select g._Refs_key, g._Term_key, g.termID, g.qualifier, g.inferredFrom, ' + \
    'g._Object_key, g.symbol, g.name, g.markerType, ' + \
    'mDate = convert(varchar(10), g.modification_date, 112), ' + \
    'markerID = ma.accID, ' + \
    'refID = b.accID, ' + \
    'eCode = rtrim(t.abbreviation), ' + \
    'assignedBy = u.login ' + \
    'into #results ' + \
    'from #gomarker g, ACC_Accession ma, ACC_Accession b, VOC_Term t, MGI_User u ' + \
    'where g._Object_key = ma._Object_key ' + \
    'and ma._MGIType_key = 2 ' + \
    'and ma.prefixPart = "MGI:" ' + \
    'and ma._LogicalDB_key = 1 ' + \
    'and ma.preferred = 1 ' + \
    'and g._Refs_key = b._Object_key ' + \
    'and b._MGIType_key = 1 ' + \
    'and b.prefixPart = "MGI:" ' + \
    'and b._LogicalDB_key = 1 ' + \
    'and g._EvidenceTerm_key = t._Term_key ' + \
    'and g._ModifiedBy_key = u._User_key', None)
db.sql('create index idx1 on #results(symbol)', None)
db.sql('create index idx2 on #results(_Refs_key)', None)

#
# resolve PubMed IDs for References
#
pubMed = {}
results = db.sql('select r._Refs_key, a.accID from #results r, ACC_Accession a ' + \
    'where r._Refs_key = a._Object_key ' + \
    'and a._MGIType_key = 1 ' + \
    'and a._LogicalDB_key = 29 ', 'auto')
for r in results:
    key = r['_Refs_key']
    value = r['accID']
    pubMed[key] = value

#
# resolve PubMed IDs for References
#
pubMed = {}
results = db.sql('select r._Refs_key, a.accID from #results r, ACC_Accession a ' + \
        'where r._Refs_key = a._Object_key ' + \
        'and a._MGIType_key = 1 ' + \
        'and a._LogicalDB_key = 29 ', 'auto')
for r in results:
    key = r['_Refs_key']
    value = r['accID']
    pubMed[key] = value

#
# Setup the cell ontology hash
#
#

clHash = {}

results = db.sql('''select distinct  mm._Marker_key as 'key', nc.note as 'gpe',
	nc.sequenceNum, vt._Term_key as 'tk', ve._Refs_key
	from MGI_NoteChunk nc, MGI_Note n, VOC_Evidence ve, 
	VOC_Annot va, MRK_Marker mm, VOC_Term vt 
	where nc.note like '%cell type%CL:%'   and nc._Note_key = n._Note_key 
	and n._Object_key = ve._AnnotEvidence_key 
	and ve._Annot_key = va._Annot_key 
	and va._Object_key = mm._Marker_key 
	and va._AnnotType_key = 1000 and va._Term_key = vt._Term_key 
	order by mm._Marker_key, sequenceNum, vt._Term_key''', 'auto')

#r1 = re.compile(r'cell.type:([^\s\\\n]*)', re.I)
r1 = re.compile(r'(CL:[0-9]{7}?)', re.I)

clPattern1 = re.compile(r'CL:', re.I)

tempString = ''
currentKey = ''
first = 1
workString = ''

for r in results:

    # The key to this hash is compound in nature, if any of the marker key, term key or reference key changes
    # we are on a new note.

    newKey = 'mk:' +str(r['key'])+'tk:'+str(r['tk'])+'rk:'+str(r['_Refs_key'])

    if currentKey != newKey:
        # Is this the first pass through the loop?
        if first != 1:
            temp = r1.findall(tempString)
            for word in temp:
                workString += word + '#'
            if workString.strip() != '':
                clHash[currentKey] = workString.strip()
        else:
            first = 0
        workString = ''
        tempString = r['gpe']
        # Construct the key for the current row (The row we are about to leave behind)
        currentKey = 'mk:' +str(r['key'])+'tk:'+str(r['tk'])+'rk:'+str(r['_Refs_key'])
    else:        
        tempString = tempString + r['gpe']

#
# Setup the isoformsProtein hash.  
# The items in this hash are seperated by a #, and are broken up later on into multiple rows.
#

isoformsProtein = {}

results = db.sql('select distinct  mm._Marker_key as "key", nc.note as "gpe", ' + \
    'nc.sequenceNum, vt._Term_key as "tk", ve._Refs_key ' + \
    'from MGI_NoteChunk nc, MGI_Note n, VOC_Evidence ve, ' + \
    'VOC_Annot va, MRK_Marker mm, VOC_Term vt ' + \
    'where nc.note like "%%gene_product:%%" and nc._Note_key = n._Note_key ' + \
    'and n._Object_key = ve._AnnotEvidence_key ' + \
    'and ve._Annot_key = va._Annot_key ' + \
    'and va._Object_key = mm._Marker_key ' + \
    'and va._AnnotType_key = 1000 and va._Term_key = vt._Term_key ' + \
    'order by mm._Marker_key, sequenceNum, vt._Term_key', 'auto')



r1 = re.compile(r'gene.product:([^\s\\\n]*)', re.I)

isoformPattern1 = re.compile(r'UniProtKB:', re.I)
isoformPattern2 = re.compile(r'protein_id', re.I)
isoformPattern3 = re.compile(r'NCBI:NP_', re.I)
isoformPattern4 = re.compile(r'NCBI:XP_', re.I)

tempString = ''
currentKey = ''
first = 1
workString = ''

for r in results:

    # The key to this hash is compound in nature, if any of the marker key, term key or reference key changes
    # we are on a new note.

    newKey = 'mk:' +str(r['key'])+'tk:'+str(r['tk'])+'rk:'+str(r['_Refs_key'])

    if currentKey != newKey:
        # Is this the first pass through the loop?
        if first != 1:
            temp = r1.findall(tempString)
            for word in temp:
                # Some of the entries are seperated by pipes, so we break on them as well.
                for word2 in word.split('|'):
                    word2 = word2.strip()
                    # Only certain patterns actually count, they are listed above.
                    if isoformPattern1.match(word2) != None or isoformPattern2.match(word2) != None \
                    or isoformPattern3.match(word2) != None or isoformPattern4.match(word2) != None:
                        workString += word2 + '#'
            # If the resulting string isn't blank add it to the hash.                        
            if workString.strip() != '':
                isoformsProtein[currentKey] = workString.strip()
        else:
            first = 0
        workString = ''
        tempString = r['gpe']
        # Construct the key for the current row (The row we are about to leave behind)
        currentKey = 'mk:' +str(r['key'])+'tk:'+str(r['tk'])+'rk:'+str(r['_Refs_key'])
    else:        
        tempString = tempString + r['gpe']

#
# Setup the protein hash
# This is as marker key <- protein text 
#

results = db.sql('select distinct mm._Marker_key, seqID=mc.accID, mc._LogicalDB_key ' + \
    'from SEQ_Marker_Cache mc, MRK_Mouse_View mm ' + \
    'where mc._Marker_key = mm._Marker_key ' + \
    'and mm._Marker_Type_key = 1 ' + \
    'and mc._Qualifier_key = 615421 ' + \
    'union ' + \
    'select distinct mm._Marker_key, seqID=mc.accID, mc._LogicalDB_key ' + \
    'from SEQ_Marker_Cache mc, MRK_Mouse_View mm ' + \
    'where mc._Marker_key = mm._Marker_key ' + \
    'and mm._Marker_Type_key = 11 ' + \
    'and mc._Qualifier_key = 615420', 'auto')


proteins = {}
proteinsGene = {}

proteinPattern1 = re.compile(r'NP_', re.I)
proteinPattern2 = re.compile(r'XP_', re.I)

for r in results:
    key = r['_Marker_key']
    
    logicalDB = r['_LogicalDB_key']    
    
    if logicalDB in [13,41]:
        proteins[key] = 'UniProtKB:' + r['seqID']
    elif logicalDB in [9]:
        proteinsGene[key] = 'EMBL:' + r['seqID'] 
    elif proteinPattern1.match(r['seqID']) != None or proteinPattern2.match(r['seqID']) != None:
        proteins[key] = 'NCBI:' + r['seqID']   
    else:
        proteinsGene[key] = 'NCBI:' + r['seqID']   

#
# process results
#
results = db.sql('select * from #results order by symbol', 'auto')
for r in results:

    reportRow = ''    

    if dag.has_key(r['_Term_key']):

	# columns 1-5
        reportRow = DBABBREV + TAB
        reportRow = reportRow + r['markerID'] + TAB
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

        # Calculate the key for this pass

        isoformKey = 'mk:'+str(r['_Object_key'])+'tk:'+str(r['_Term_key'])+'rk:'+str(r['_Refs_key'])

        # column 12 is populated in a special way
        # If there is a isoform or a protein, print out protein
        # otherwise print out the marker type.

        if isoformsProtein.has_key(isoformKey) or proteins.has_key(r['_Object_key']):
            reportRow = reportRow + 'protein' + TAB
        else:
            reportRow = reportRow + r['markerType'] + TAB
                
	# column 13
        reportRow = reportRow + SPECIES + TAB

	# column 14
        reportRow = reportRow + r['mDate'] + TAB

	# column 15; assigned by
        if r['assignedBy'] == UNIPROTKB:
            reportRow = reportRow + 'UniProtKB' + TAB

	# remove "GOA_"; for example:  "GOA_IntAct" ==> "IntAct"
        elif string.find(r['assignedBy'], 'GOA_') >= 0:
            assignedBy = re.sub('GOA_', '', r['assignedBy'])
            reportRow = reportRow + assignedBy + TAB

	# check list
        elif r['assignedBy'] in assignedByList:
            reportRow = reportRow + r['assignedBy']

	# else use default
        else:
            reportRow = reportRow + DBABBREV + TAB

	# column 16
        # cell ontology: 'occurs_in' or 'part_of'

        if clHash.has_key(isoformKey):
            startPart = 'occurs_in('
            if dag[r['_Term_key']] == 'C':
                startPart = 'part_of('
            row = ''
            for word3 in clHash[isoformKey].split('#'):
                word3 = word3.replace(';', '')
                if word3 != '' and row == '':
                    row = startPart + str(word3) + ')'
                elif word3 != '': 
                    row = row + '|' + startPart + str(word3) + ')'
            reportRow = reportRow + row + TAB
        else:                    
            reportRow = reportRow + '' + TAB

        # column 17 is populated in a special way.  
        # If there is an isoform, use that, parsing through it to check for multiple annotations.
        # If not, if there is a protein use that.
        # If not, leave it blank for now.

        if isoformsProtein.has_key(isoformKey):
            for word3 in isoformsProtein[isoformKey].split('#'):
                if word3 != '':
                    row = str(word3)
                    fp.write(reportRow + row + CRT)
        else:
            if proteins.has_key(r['_Object_key']):
                row = str(proteins[r['_Object_key']])
                #row = 'Protein: ' + str(proteins[r['_Object_key']])
            elif proteinsGene.has_key(r['_Object_key']):
                row = str(proteinsGene[r['_Object_key']])
                #row = 'Protein gene: ' +str(proteinsGene[r['_Object_key']])
            else:
                row = ''
            fp.write(reportRow + row + CRT)
    
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
