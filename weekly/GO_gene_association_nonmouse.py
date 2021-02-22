
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
# !1  DB                              required        1       UniProtKB\n')
# !2  DB Object ID                    required        1       P12345\n')
# !3  DB Object Symbol                required        1       PHO3\n')
# !4  Qualifier                       required        1 or 2  NOT|involved_in\n')
# !5  GO ID                           required        1       GO:0003993\n')
# !6  DB:Reference (|DB:Reference)    required        1 or greater    PMID:2676709\n')
# !7  Evidence Code                   required        1       IMP\n')
# !8  With (or) From                  optional        0 or greater    GO:0000346\n')
# !9  Aspect                          required        1       F\n')
# !10 DB Object Name                  optional        0 or 1  Toll-like receptor 4\n')
# !11 DB Object Synonym (|Synonym)    optional        0 or greater    hToll\n')
# !12 DB Object Type                  required        1       protein\n')
# !13 Taxon(|taxon)                   required        1 or 2  taxon:9606\n')
# !14 Date                            required        1       20090118\n')
# !15 Assigned By                     required        1       SGD\n')
# !16 Annotation Extension            optional        0 or greater    part_of(CL:0000576)\n')
# !17 Gene Product Form ID            optional        0 or 1  UniProtKB:P12345-2\n')
#
# exclude J:88213 (olfactory load)
# exclude J:104715 (tbreddy's Rat)
#
# History:
#
# lec	01/13/2014
#	- TR11570/fix GO qualifier; use VOC_Term
#
# 03/26/2013
#	- TR11330 (scrum-dog TR11248)/added "gaf" header
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
import db

db.setTrace()

FIELD1 = 'UniProt'
FIELD12 = 'protein'
FIELD15 = 'MGI'

TAB = reportlib.TAB
CRT = reportlib.CRT

dagQualifierGAF = {'C':'located_in', 'P':'acts_upstream_of_or_within', 'F':'enables'}

def writeRecord(i, r, e):

        # if we can't find the DAG for the Term, skip it

        if r['_Term_key'] not in dag:
                return

        # field 1
        fp.write(FIELD1 + TAB)

        # field 2
        fp.write(i + TAB)

        # field 3
        fp.write(TAB)

        # field 4
        qualifier = ""
        if r['termID'] == 'GO:0008150':
                qualifier = 'involved_in'
        elif r['termID'] == 'GO:0003674':
                qualifier = 'enables'
        elif r['termID'] == 'GO:0005575':
                qualifier = 'is_active_in'
        elif r['_AnnotEvidence_key'] in goQualifierGAF:
            if r['qualifier'] == 'NOT':
               qualifier = 'NOT|'
            qualifier = qualifier + '|'.join(goQualifierGAF[r['_AnnotEvidence_key']])
        elif r['qualifier'] == 'NOT':
            qualifier = 'NOT|' + dagQualifierGAF[dag[r['_Term_key']]]
        elif r['qualifier'] != None:
            qualifier = r['qualifier'].strip()
        elif r['uniprotIDs'] != None and r['uniprotIDs'].find('InterPro:') >= 0 and dag[r['_Term_key']] == 'P':
            qualifier = 'involved_in'
        else:
            qualifier = dagQualifierGAF[dag[r['_Term_key']]]
        fp.write(qualifier + TAB)

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
        fp.write(FIELD15 + TAB)

        fp.write(TAB)
        fp.write(CRT)

#
# Main
#

db.useOneConnection(1)
fp = reportlib.init('gene_association', fileExt = '.mgi_nonmouse', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp.write('!gpa-version: 2.2\n')
fp.write('!generated-by: MGI\n')
fp.write('!date-generated: %s\n' % (mgi_utils.date("%Y-%m-%d")))

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
        select a._Term_key, ta.accID as termID, q.term as qualifier, a._Object_key, 
                 e._AnnotEvidence_key, e.inferredFrom as uniprotIDs, 
                 e.modification_date, e._Refs_key
        into temporary table gomarker1 
        from VOC_Annot a, ACC_Accession ta, VOC_Term t, VOC_Evidence e, 
             VOC_Term et, VOC_Term q, MGI_User u , MRK_Marker m
        where a._AnnotType_key = 1000 
        and a._Annot_key = e._Annot_key 
        and a._Term_key = t._Term_key 
        and a._Term_key = ta._Object_key 
        and ta._MGIType_key = 13 
        and ta.preferred = 1 
        and e._EvidenceTerm_key = et._Term_key 
        and et.abbreviation in ('ISS','ISO','ISM','ISA') 
        and e.inferredFrom like 'UniProtKB:%' 
        and a._Qualifier_key = q._Term_key 
        and e._ModifiedBy_key = u._User_key
        and e._Refs_key not in (89196,105787) 
        and u.login not in ('GOC', 'RGD', 'UniProtKB')
        and u.login not like 'GOA%'
        and a._Object_key = m._Marker_key
        ''', None)

db.sql('create index gomarker1_idx1 on gomarker1(_Object_key)', None)
db.sql('create index gomarker1_idx2 on gomarker1(_Refs_key)', None)

#
# resolve pub med id
#
db.sql('''select g._AnnotEvidence_key, g._Term_key, g.termID, g.qualifier, g.uniprotIDs,
                 to_char(g.modification_date, 'YYYYMMDD') as mDate,
                 b.accID as refID
        into temporary table gomarker2
        from gomarker1 g LEFT OUTER JOIN ACC_Accession b on
                (g._Refs_key = b._Object_key 
                and b._MGIType_key = 1 
                and b._LogicalDB_key = 29)
        ''', None)
db.sql('create index gomarker2_idx1 on gomarker2(_AnnotEvidence_key)', None)

#
# properties
# external ref = 6481778
#
results = db.sql('''select g._AnnotEvidence_key, p.value
        from gomarker2 g, VOC_Evidence_Property p
        where g._AnnotEvidence_key = p._AnnotEvidence_key
        and p._PropertyTerm_key = 6481778
        order by g._AnnotEvidence_key, p.stanza, p.sequenceNum''', 'auto')

evidence = {}
for r in results:

    eKey = r['_AnnotEvidence_key']
    value = r['value']
    value = value.replace(';', '|')

    # parse it for pmid, evidence code, and cross-reference

    tokens = str.split(value, '|')

    pmid = tokens[0] 
    ecode = ''
    dbxref = ''

    if len(tokens) > 1:
        ecode = str.strip(tokens[1])

    if len(tokens) > 2:
        dbxref = str.strip(tokens[2])

    dictvalue = (pmid, ecode, dbxref)

    if eKey not in evidence:
        evidence[eKey] = []
    evidence[eKey].append(dictvalue)

#
# goQualifierGAF : go_qualifier
#
goQualifierGAF = {}
results = db.sql('''select distinct a._AnnotEvidence_key, t.term, p.value, t2.note
    from gomarker2 a,
        VOC_Evidence_Property p,  
        VOC_Term t,
        VOC_Term t2
    where a._AnnotEvidence_key = p._AnnotEvidence_key
    and p._PropertyTerm_key = t._Term_key
    and t.term in ('go_qualifier')
    and p.value = t2.term
    ''', 'auto')
for r in results:
    key = r['_AnnotEvidence_key']
    value = r['value']
    if key not in goQualifierGAF:
        goQualifierGAF[key] = []
    goQualifierGAF[key].append(value)

#print(goQualifierGAF)

#
# process results
#

results = db.sql('select * from gomarker2 order by uniprotIDs', 'auto')

for r in results:

    # there may be multiple instances of UniProt ids
    # write out one record per UniProt id

    tokens = str.split(r['uniprotIDs'], '|')
    ids = []

    for t in tokens:
        if str.find(t, 'UniProtKB:') > -1:
            id = str.split(t, 'UniProtKB:')
            ids.append(id[1])

    for i in ids:

        if len(i) == 0:
            continue

        # get rid of any dangling delimiters

        i = i.replace('|', '')

        eKey = r['_AnnotEvidence_key']

        # make up a bogus evidence record if there isn't one

        if eKey not in evidence:
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
