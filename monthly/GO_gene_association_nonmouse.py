#!/usr/local/bin/python

'''
#
# GO_gene_association_nonmouse.py
#
# Report:
#       Tab-delimited file of all ISS/UniProt MGI GO/Marker associations
#
# Usage:
#       GO_gene_association_nonmouse.py
#
# Used by:
#	GOA
#
# Output format:
#
# Special format for this report:
#
#   1.  Database designation (UniProt)
#   2.  Inferred From (id minus the UniProt prefix)
#   3.  NOT
#   4.  null
#   5.  GO id
#   6.  PubMed ID of Reference from GO notes (external ref:)
#   7.  Evidence abbreviation (external ref:)
#   8.  DBXRef (external ref:)
#   9.  GO DAG Abbreviation (F, P, C)
#   10. null
#   11. null
#   12. protein
#   13. null
#   14. Modification Date (YYYYMMDD)
#   15. Assigned By (MGI)
#
# exclude J:88213 (olfactory load)
#
# History:
#
# 01/30/2006	lec
#	- TR 7424; modified to new format
#
# 08/22/2005	lec
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

FIELD1 = 'UniProt'
FIELD12 = 'protein'
FIELD15 = 'MGI'

TAB = reportlib.TAB
CRT = reportlib.CRT

def writeRecord(i, r, e):

	# if we can't find the DAG for the Term, skip it

	if not dag.has_key(r['_Term_key']):
		return

	# field 1
	fp.write(FIELD1 + TAB)

	# field 2
	fp.write(i + TAB)

	# field 3
	if r['isNot'] == 1:
		fp.write('NOT')
	fp.write(TAB)

	# field 4
	fp.write(TAB)

	# field 5
	fp.write(r['termID'] + TAB)

	# field 6
	fp.write(mgi_utils.prvalue(e[0]) + TAB)

	# field 7
	fp.write(mgi_utils.prvalue(e[1]) + TAB)

	# field 8 
	fp.write(mgi_utils.prvalue(e[2]) + TAB)

	# field 9
	fp.write(dag[r['_Term_key']] + TAB)

	# field 10
	fp.write(TAB)

	# field 11
	fp.write(TAB)

	# field 12
	fp.write(FIELD12 + TAB)

	# field 13
	fp.write(TAB)

	# field 14
	fp.write(r['mDate'] + TAB)

	# field 15
	fp.write(FIELD15)
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
# retrieve all ISS annotations that have a "with" value that begins "UniProt"
#
db.sql('select a._Term_key, termID = ta.accID, a.isNot, a._Object_key, ' + \
	'e._AnnotEvidence_key, uniprotIDs = e.inferredFrom, e.modification_date, e._Refs_key, e._ModifiedBy_key ' + \
	'into #gomarker ' + \
	'from VOC_Annot a, ACC_Accession ta, VOC_Term t, VOC_Evidence e, VOC_Term et ' + \
	'where a._AnnotType_key = 1000 ' + \
	'and a._Annot_key = e._Annot_key ' + \
	'and a._Term_key = t._Term_key ' + \
	'and a._Term_key = ta._Object_key ' + \
	'and ta._MGIType_key = 13 ' + \
	'and ta.preferred = 1 ' + \
	'and e._EvidenceTerm_key = et._Term_key ' + \
	'and et.abbreviation = "ISS" ' + \
	'and e.inferredFrom like "UniProt:%" ' + \
	'and e._Refs_key not in (89196)' , None)
db.sql('create index idx1 on #gomarker(_Object_key)', None)
db.sql('create index idx2 on #gomarker(_Refs_key)', None)

#
# resolve foreign keys (reference ID)
#
db.sql('select g._AnnotEvidence_key, g._Term_key, g.termID, g.isNot, g.uniprotIDs, g._ModifiedBy_key, ' + \
	'mDate = convert(varchar(10), g.modification_date, 112), ' + \
	'refID = b.accID ' + \
	'into #results ' + \
	'from #gomarker g, ACC_Accession b ' + \
	'where g._Refs_key *= b._Object_key ' + \
	'and b._MGIType_key = 1 ' + \
	'and b._LogicalDB_key = 29', None)
db.sql('create index idx1 on #results(_AnnotEvidence_key)', None)

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
    value = string.strip(regsub.gsub('\n', '', r['note']))
    if not allnotes.has_key(key):
	     allnotes[key] = []
    allnotes[key].append(value)

evidence = {}
for n in allnotes.keys():
    value = string.join(allnotes[n], '')

    # grab all text between "external ref:" and "text:"

    i = string.find(value, 'external ref:')
    j = string.find(value, 'text:')

    # parse it for pmid, evidence code, and cross-reference

    if j > i and i >= 0 and len(value) > 0:

	s1 = value[i + 13:j]

	# split by 'external ref:' (may be more than one)

	s1tokens = string.split(s1, 'external ref:')

	for s in s1tokens:
	    s2tokens = string.split(s, '|')

            pmid = ''
            ecode = ''
            dbxref = ''

	    pmid = s2tokens[0]
	    if len(s2tokens) > 1:
	        ecode = s2tokens[1]
	    if len(s2tokens) > 2:
	        dbxref = s2tokens[2]

	    dictvalue = (pmid, ecode, dbxref)

	    if not evidence.has_key(n):
	        evidence[n] = []
	    evidence[n].append(dictvalue)

#
# process results
#

results = db.sql('select * from #results order by uniprotIDs', 'auto')

for r in results:

    ids = string.split(r['uniprotIDs'], 'UniProt:')

    for i in ids:

	if len(i) == 0:
	    continue

        eKey = r['_AnnotEvidence_key']

	# if no evidence (no pub med it) and "tbreddy", skip it
        if not evidence.has_key(eKey) and r['_ModifiedBy_key'] == 1095:
	    continue

        # make up a bogus evidence record if there isn't one
        if not evidence.has_key(eKey):
	    evidence[eKey] = []
	    if r['refID'] == None:
	        value = ('', '', '')
	    else:
	        value = ('PMID:' + r['refID'], '', '')
	    evidence[eKey].append(value)

        for e in evidence[eKey]:
            writeRecord(i, r, e)

reportlib.finish_nonps(fp)
db.useOneConnection(0)

