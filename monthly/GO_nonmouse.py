#!/usr/local/bin/python

'''
#
# GO_nonmouse.py
#
# Report:
#       Tab-delimited file of all ISS MGI GO/Marker associations
#
# Usage:
#       GO_nonmouse.py
#
# Used by:
#	GO_nonmouse.py
#
# Output format:
#
#   Symbol
#
#	MGI Gene Symbol
#	GO ID
#	GO Term
#	Evidence Code
#	PubMed ID
#	Inferred From
#	Notes
#
# History:
#
# lec	10/27/2004
#	- new
#
'''

import sys
import os
import string
import regsub
import db
import reportlib
import mgi_utils

TAB = reportlib.TAB
CRT = reportlib.CRT

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []
cmds.append('select a._Term_key, a.term, termID = a.accID, ' + \
	'm.symbol, m._Marker_key, markerID = ma.accID, ' + \
	'e._AnnotEvidence_key, e.inferredFrom, eCode = rtrim(t.abbreviation), ' + \
	'refID = b.accID ' + \
	'into #go ' + \
	'from VOC_Annot_View a, MRK_Marker m, ACC_Accession ma, ' + \
	'VOC_Evidence e, VOC_Term t, ACC_Accession b ' + \
	'where a._AnnotType_key = 1000 ' + \
	'and a._Object_key = m._Marker_key ' + \
	'and m._Marker_key = ma._Object_key ' + \
	'and ma._MGIType_key = 2 ' + \
	'and ma.prefixPart = "MGI:" ' + \
	'and ma._LogicalDB_key = 1 ' + \
	'and ma.preferred = 1 ' + \
	'and a._Annot_key = e._Annot_key ' + \
	'and e._EvidenceTerm_key = t._Term_key ' + \
	'and t.abbreviation = "ISS" ' + \
	'and e.inferredFrom like "%UniProt%" ' + \
	'and e._Refs_key not in (80961, 89196) ' + \
	'and e._Refs_key *= b._Object_key ' + \
	'and b._MGIType_key = 1 ' + \
	'and b._LogicalDB_key = 29 ')
cmds.append('create index idx1 on #go(_AnnotEvidence_key)')
cmds.append('create index idx2 on #go(symbol)')
db.sql(cmds, None)

results = db.sql('select g._AnnotEvidence_key, c.note, c.sequenceNum ' + \
	'from #go g, MGI_Note n, MGI_NoteChunk c ' + \
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

results = db.sql('select * from #go order by symbol', 'auto')
for r in results:
	fp.write(r['symbol'] + TAB)
	fp.write(r['termID'] + TAB)
	fp.write(r['term'] + TAB)
	fp.write(r['eCode'] + TAB)
	fp.write(mgi_utils.prvalue(r['refID']) + TAB)
	fp.write(mgi_utils.prvalue(r['inferredFrom']) + TAB)

	if notes.has_key(r['_AnnotEvidence_key']):
	    fp.write(string.join(notes[r['_AnnotEvidence_key']], '') + TAB)

	fp.write(CRT)
	
reportlib.finish_nonps(fp)

