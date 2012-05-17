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
# exclude J:104715 (tbreddy's Rat)
#
# History:
#
# 12/28/2011	lec
#	- changed non-ansi-standard query to left outer join
#
# 06/27/2011	lec
#	- TR10044/replace notes with properties
#
# 10/13/2010	lec
#	- TR 10393; exclude UniProtKB (GOA/Human)
#
# 06/22/2010	lec
#	- TR 10260; multiple start/end notes after "external ref"
#	  exclude GOA, RGD, GOC
#	  J:88213 (olfactory load) is still exlucded (bobs)
#	  J:104715 (tbreddy's Rat) is still excluded (tbreddy)
#
# 04/13/2010	lec
#	- TR 10163; skip ISO/J:155856
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
	fp.write(str(r['mDate']) + TAB)

	# field 15
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
results = db.sql('select distinct _Object_key, rtrim(dagAbbrev) as dagAbbrev from DAG_Node_View where _Vocab_key = 4', 'auto')
dag = {}
for r in results:
	dag[r['_Object_key']] = r['dagAbbrev']

#
# retrieve all ISS,ISO, ISM, ISA annotations 
# that have a "with" value 
# that begins "UniProtKB"
#
db.sql('''
	select a._Term_key, ta.accID as termID, q.synonym as qualifier, a._Object_key, 
	         e._AnnotEvidence_key, e.inferredFrom as uniprotIDs, 
	         e.modification_date, e._Refs_key
	into #gomarker1 
	from VOC_Annot a, ACC_Accession ta, VOC_Term t, VOC_Evidence e, 
	     VOC_Term et, MGI_Synonym q, MGI_User u 
	where a._AnnotType_key = 1000 
	and a._Annot_key = e._Annot_key 
	and a._Term_key = t._Term_key 
	and a._Term_key = ta._Object_key 
	and ta._MGIType_key = 13 
	and ta.preferred = 1 
	and e._EvidenceTerm_key = et._Term_key 
	and et.abbreviation in ('ISS','ISO','ISM','ISA') 
	and e.inferredFrom like 'UniProtKB:%' 
        and a._Qualifier_key = q._Object_key 
	and q._SynonymType_key = 1023 
	and e._ModifiedBy_key = u._User_key
	and e._Refs_key not in (89196,105787) 
	and u.login not in ('GOC', 'RGD', 'UniProtKB')
	and u.login not like 'GOA%'
	''', None)

db.sql('create index gomarker1_idx1 on #gomarker1(_Object_key)', None)
db.sql('create index gomarker1_idx2 on #gomarker1(_Refs_key)', None)

#
# resolve pub med id
#
db.sql('''select g._AnnotEvidence_key, g._Term_key, g.termID, g.qualifier, g.uniprotIDs,
	         convert(varchar(10), g.modification_date, 112) as mDate,
	         b.accID as refID
	into #gomarker2
	from #gomarker1 g LEFT OUTER JOIN ACC_Accession b on
		(g._Refs_key = b._Object_key 
		and b._MGIType_key = 1 
		and b._LogicalDB_key = 29)
	''', None)
db.sql('create index gomarker2_idx1 on #gomarker2(_AnnotEvidence_key)', None)

#
# properties
# external ref = 6481778
#
results = db.sql('''select g._AnnotEvidence_key, p.value
        from #gomarker2 g, VOC_Evidence_Property p
        where g._AnnotEvidence_key = p._AnnotEvidence_key
        and p._PropertyTerm_key = 6481778
        order by g._AnnotEvidence_key, p.stanza, p.sequenceNum''', 'auto')

evidence = {}
for r in results:

    eKey = r['_AnnotEvidence_key']
    value = r['value']
    value = value.replace(';', '|')

    # parse it for pmid, evidence code, and cross-reference

    tokens = string.split(value, '|')

    pmid = tokens[0] 
    ecode = ''
    dbxref = ''

    if len(tokens) > 1:
	ecode = string.strip(tokens[1])

    if len(tokens) > 2:
	dbxref = string.strip(tokens[2])

    dictvalue = (pmid, ecode, dbxref)

    if not evidence.has_key(eKey):
        evidence[eKey] = []
    evidence[eKey].append(dictvalue)

#
# process results
#

results = db.sql('select * from #gomarker2 order by uniprotIDs', 'auto')

for r in results:

    # there may be multiple instances of UniProt ids
    # write out one record per UniProt id

    tokens = string.split(r['uniprotIDs'], '|')
    ids = []

    for t in tokens:
        if string.find(t, 'UniProtKB:') > -1:
            id = string.split(t, 'UniProtKB:')
	    ids.append(id[1])

    for i in ids:

	if len(i) == 0:
	    continue

	# get rid of any dangling delimiters

        i = i.replace('|', '')

        eKey = r['_AnnotEvidence_key']

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

