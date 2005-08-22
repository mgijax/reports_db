#!/usr/local/bin/python

'''
#
# GO_gene_association_nonmouse.py
#
# Report:
#       Tab-delimited file of all ISS MGI GO/Marker associations
#
# Usage:
#       GO_gene_association_nonmouse.py
#
# Used by:
#	GOA
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
#   12. Marker Type (gene)
#   13. Species (taxon:10090)
#   14. Modification Date (YYYYMMDD)
#   15. Assigned By
#
# History:
#
# lec	08/22/2005
#       - converted to standard GO format per request (David Hill)
#
'''

import sys
import os
import string
import regsub
import db
import reportlib
import mgi_utils

DBABBREV = 'MGI'
SPECIES = 'taxon:10090'

# mapping between MGI Marker Type and what gets printed in the gene association file
markerTypes = {1: 'gene'}

TAB = reportlib.TAB
CRT = reportlib.CRT

def writeRecord(r, refID):

	# if we can't find the DAG for the Term, skip it

	if not dag.has_key(r['_Term_key']):
		return

	fp.write(DBABBREV + TAB)
	fp.write(r['markerID'] + TAB)
	fp.write(r['symbol'] + TAB)

	if r['isNot'] == 1:
		fp.write('NOT')

	fp.write(TAB)
	fp.write(r['termID'] + TAB)
	fp.write(refID + TAB)
	fp.write(r['eCode'] + TAB)

	# substitute | for ", " in inferredFrom

	if r['inferredFrom'] != None:
		inferredFrom = regsub.gsub(',', '|', r['inferredFrom'])
		inferredFrom = regsub.gsub(';', '|', inferredFrom)
		inferredFrom = regsub.gsub(' ', '', inferredFrom)
	else:
		inferredFrom = r['inferredFrom']

	fp.write(mgi_utils.prvalue(inferredFrom) + TAB)

	fp.write(dag[r['_Term_key']] + TAB)
	fp.write(r['name'] + TAB)

	if syns.has_key(r['_Object_key']):
		fp.write(string.join(syns[r['_Object_key']], '|'))

	fp.write(TAB)

	if markerTypes.has_key(r['_Marker_Type_key']):
		fp.write(markerTypes[r['_Marker_Type_key']] + TAB)

	fp.write(SPECIES + TAB)

	fp.write(r['mDate'] + TAB)

	if r['modifiedBy'] == 'swissload':
		fp.write('SWALL')
	else:
		fp.write(DBABBREV)

	fp.write(CRT)

#
# Main
#

db.useOneConnection(1)
fp = reportlib.init('gene_association_nonmouse', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

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
db.sql('select a._Term_key, t.term, termID = ta.accID, a.isNot, a._Object_key, ' + \
	'e._AnnotEvidence_key, e.inferredFrom, e.modification_date, e._EvidenceTerm_key, e._Refs_key, e._ModifiedBy_key, ' + \
	'm._Marker_Type_key, m.symbol, m.name ' + \
	'into #gomarker ' + \
	'from VOC_Annot a, ACC_Accession ta, VOC_Term t, VOC_Evidence e, VOC_Term et,  MRK_Marker m ' + \
	'where a._AnnotType_key = 1000 ' + \
	'and a._Annot_key = e._Annot_key ' + \
	'and a._Object_key = m._Marker_key ' + \
	'and m._Marker_Type_key = 1 ' + \
	'and a._Term_key = t._Term_key ' + \
	'and a._Term_key = ta._Object_key ' + \
	'and ta._MGIType_key = 13 ' + \
	'and ta.preferred = 1 ' + \
	'and e._EvidenceTerm_key = et._Term_key ' + \
	'and et.abbreviation = "ISS" ' + \
	'and e.inferredFrom like "%UniProt%" ' + \
	'and e._Refs_key not in (80961, 89196)' , None)
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
db.sql('select g._AnnotEvidence_key, g._Term_key, g.termID, g.isNot, g.inferredFrom, ' + \
	'g._Object_key, g._Marker_Type_key, g.symbol, g.name, ' + \
	'mDate = convert(varchar(10), g.modification_date, 112), ' + \
	'markerID = ma.accID, ' + \
	'refID = b.accID, ' + \
	'eCode = rtrim(t.abbreviation), ' + \
	'modifiedBy = u.login ' + \
	'into #results ' + \
	'from #gomarker g, ACC_Accession ma, ACC_Accession b, VOC_Term t, MGI_User u ' + \
	'where g._Object_key = ma._Object_key ' + \
	'and ma._MGIType_key = 2 ' + \
	'and ma.prefixPart = "MGI:" ' + \
	'and ma._LogicalDB_key = 1 ' + \
	'and ma.preferred = 1 ' + \
	'and g._Refs_key *= b._Object_key ' + \
	'and b._MGIType_key = 1 ' + \
	'and b._LogicalDB_key = 29 ' + \
	'and g._EvidenceTerm_key = t._Term_key ' + \
	'and g._ModifiedBy_key = u._User_key', None)
db.sql('create index idx1 on #results(symbol)', None)
db.sql('create index idx2 on #results(_AnnotEvidence_key)', None)

#
# notes
#
results = db.sql('select g._AnnotEvidence_key, c.note, c.sequenceNum ' + \
	'from #results g, MGI_Note n, MGI_NoteChunk c ' + \
	'where g._AnnotEvidence_key = n._Object_key ' + \
	'and n._MGIType_key = 25 ' + \
	'and n._Note_key = c._Note_key ' + \
	'order by g._AnnotEvidence_key, c.sequenceNum', 'auto')
allnotes = {}
for r in results:
    key = r['_AnnotEvidence_key']
    value = r['note']
    value = regsub.gsub('PMID: ', 'PMID:', value)
    value = regsub.gsub('\n', '', value)
    if not allnotes.has_key(key):
	     allnotes[key] = []
    allnotes[key].append(value)

notes = {}
for n in allnotes.keys():
    value = string.join(allnotes[n])
    pmids = []

    i = string.find(value, 'PMID:')
    while i >= 0:
	t = value[i:]
	j = 5
	while j < len(t) and t[j] in string.digits:
	    j = j + 1
	pmid = t[:j]
	value = t[j:]
        i = string.find(value, 'PMID:')
	if pmid not in pmids:
	    pmids.append(pmid)

    if len(pmids) > 0:
        if not notes.has_key(n):
	     notes[n] = []
        notes[n].append(string.join(pmids, ';'))

#
# process results
#

results = db.sql('select * from #results order by symbol', 'auto')

for r in results:

	# if the notes have a pub med id, use it
	if notes.has_key(r['_AnnotEvidence_key']):
	    tokens = string.split(string.join(notes[r['_AnnotEvidence_key']], ''), ';')
	    for t in tokens:
	        writeRecord(r, t)

	# else use the annotation reference
	else:
	    writeRecord(r, 'PMID:' + mgi_utils.prvalue(r['refID']))
	
reportlib.finish_nonps(fp)
db.useOneConnection(0)

