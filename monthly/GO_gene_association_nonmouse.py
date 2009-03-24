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
#   1.  Database designation (UniProt) required
#   2.  DB_Object_ID (Inferred From (id minus the UniProt prefix))
#   3.  DB_Object_Symbol (null) optional
#   4.  Qualifier (NOT) optional
#   5.  GO id required
#   6.  DB:Reference (|DB:Reference) (PubMed ID of Reference from GO notes) required
#   7.  Evidence code required
#   8.  With (or) From	optional
#   9.  Aspect required (GO DAG Abbreviation (F, P, C))
#   10. DB_Object_Name optional (null)
#   11. DB_Object_Synonym optional (null)
#   12. DB_Object_Type required (protein)
#   13. taxon (|taxon) required (null)
#   14. Modification Date required (YYYYMMDD)
#   15. Assigned By required (MGI)
#
# exclude J:88213 (olfactory load)
#
# History:
#
# 03/24/2009	lec
#	- TR 9569; fix columns 3,4,8,13
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
import re
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
	fp.write(TAB)

	# field 4
	fp.write(string.strip(r['qualifier']) + TAB)

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
	modifiedBy = r['login']
	if modifiedBy[0:4] == 'GOA_':
		fp.write(modifiedBy[4:])
	else:
		fp.write(FIELD15)
	fp.write(CRT)

#
# Main
#

db.useOneConnection(1)
fp = reportlib.init('gene_association', fileExt = '.mgi_nonmouse', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# retrieve all dag abbrevations for each term
#
results = db.sql('select distinct _Object_key, dagAbbrev = rtrim(dagAbbrev) from DAG_Node_View where _Vocab_key = 4', 'auto')
dag = {}
for r in results:
	dag[r['_Object_key']] = r['dagAbbrev']

#
# retrieve all ISS annotations that have a "with" value that begins "UniProtKB"
#
db.sql('select a._Term_key, termID = ta.accID, qualifier = q.synonym, a._Object_key, ' + \
	'e._AnnotEvidence_key, uniprotIDs = e.inferredFrom, e.modification_date, e._Refs_key, e._ModifiedBy_key, u.login ' + \
	'into #gomarker ' + \
	'from VOC_Annot a, ACC_Accession ta, VOC_Term t, VOC_Evidence e, VOC_Term et, MGI_Synonym q, MGI_User u ' + \
	'where a._AnnotType_key = 1000 ' + \
	'and a._Annot_key = e._Annot_key ' + \
	'and a._Term_key = t._Term_key ' + \
	'and a._Term_key = ta._Object_key ' + \
	'and ta._MGIType_key = 13 ' + \
	'and ta.preferred = 1 ' + \
	'and e._EvidenceTerm_key = et._Term_key ' + \
	'and et.abbreviation in ("ISS","ISO","ISM","ISA") ' + \
	'and e.inferredFrom like "UniProtKB:%" ' + \
	'and e._Refs_key not in (89196) ' + \
        'and a._Qualifier_key = q._Object_key ' + \
	'and q._SynonymType_key = 1023 ' + \
	'and e._ModifiedBy_key = u._User_key', None)

db.sql('create index idx1 on #gomarker(_Object_key)', None)
db.sql('create index idx2 on #gomarker(_Refs_key)', None)

#
# resolve pub med id
#
db.sql('select g._AnnotEvidence_key, g._Term_key, g.termID, g.qualifier, g.uniprotIDs, g._ModifiedBy_key, g.login, ' + \
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
    value = string.strip(re.sub('\n', '', r['note']))
    if not allnotes.has_key(key):
	     allnotes[key] = []
    allnotes[key].append(value)

evidence = {}
for n in allnotes.keys():
    value = string.join(allnotes[n], '')

    # grab all text between "external ref:" and "text:"
    # there may be multiple instances of "external ref:"

    i = string.find(value, 'external ref:')
    j = string.find(value, 'text:')

    # parse it for pmid, evidence code, and cross-reference

    if j > i and i >= 0 and len(value) > 0:

	s1 = value[i + 13:j]

	# split by 'external ref:'

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

	    if len(pmid) > 0:
	        dictvalue = (pmid, ecode, dbxref)

	        if not evidence.has_key(n):
	            evidence[n] = []
	        evidence[n].append(dictvalue)

#
# process results
#

results = db.sql('select * from #results order by uniprotIDs', 'auto')

for r in results:

    # there may be multiple instances of UniProt ids
    # write out one record per UniProt id

    ids = string.split(r['uniprotIDs'], 'UniProtKB:')

    for i in ids:

	if len(i) == 0:
	    continue

	# get rid of any dangling delimiters

        i = i.replace('|', '')

        eKey = r['_AnnotEvidence_key']

	# if no evidence (no pub med id) and "tbreddy", skip it
	# these are Rat Genome (J:104715)
        if not evidence.has_key(eKey) and r['_ModifiedBy_key'] == 1095:
	    continue

        # make up a bogus evidence record if there isn't one

        if not evidence.has_key(eKey):
	    evidence[eKey] = []
	    if r['refID'] == None:
	        value = ('', '', '')
	    else:
	        value = ('PMID:' + r['refID'], 'IDA', '')
	    evidence[eKey].append(value)

	# write out one record per External Reference

        for e in evidence[eKey]:
            writeRecord(i, r, e)

reportlib.finish_nonps(fp)
db.useOneConnection(0)

