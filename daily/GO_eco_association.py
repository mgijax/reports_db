#!/usr/bin/env python

'''
#
# GO_eco_association.py
#
# Report:
#       
# Generate a report of GO annotations that have
# the ECO code int the notes. This report will be
# publicly available and used by the gaf_fprocessor programs
# to map the Evidence code to the ECO code
#
#
# Fields:
# 1) mgiid
# 2) symbol
# 3) qualifier (GO term) 
# 4) goid
# 5) gaf_qualifier
# 6) evcode,
# 7) ref (J#)
# 8) mgiref
# 9) pubmedid
# 10)author
# 11)value (eco note) 
# 12)inferredfrom
# 13)stanza
# 14)sequencenum
#
# Usage:
#      GO_eco_association.py 
#
# History:
#
# 04/16/2015	lnh
#	- TR 11992
#
'''
 
import sys 
import os
import string
import reportlib
from datetime import datetime

#
#Set up email contact
#
sender = "mgiadmin@lindon.informatics.jax.org"
receiver = "Mgi-go@jax.org"
#receiver = "lnh@jax.org"

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db


CRT = reportlib.CRT
TAB = reportlib.TAB

def printResults(cmd):
    
    results = db.sql(cmd, 'auto')
    fp.write('#Total: %s\n' % (len(results)))
    log.write('#Total: %s\n' % (len(results)))
    fp.write('mgiid' + TAB)
    fp.write('symbol' + TAB)
    fp.write('go_term' + TAB)
    fp.write('goid' + TAB)
    fp.write('gaf_qualifier' + TAB)
    fp.write('evcode' + TAB)
    fp.write('ref' + TAB)
    fp.write('mgiref' + TAB)
    fp.write('pubmedid' + TAB)
    fp.write('author' + TAB)
    fp.write('evidence_note' + TAB)
    fp.write('ecoCode' + TAB)
    fp.write('inferredfrom' + TAB)
    fp.write('stanza' + TAB)
    fp.write('sequencenum' + 2*CRT)
    for r in results:
        fp.write(r['mgiid'].strip() + TAB)
        fp.write(r['symbol'].strip() + TAB)
        fp.write(r['go_term'].strip() + TAB)
        fp.write(r['goid'].strip() + TAB)
        if r['gaf_qualifier']:
            fp.write(r['gaf_qualifier'].strip() + TAB)
        else: fp.write("" + TAB)
        fp.write(r['evCode'] + TAB)
        fp.write(r['ref'] + TAB)
        fp.write(r['mgiref'] + TAB)
        if r['pubmedid']:
            fp.write("PMID:"+r['pubmedid'] + TAB)
        else: fp.write("" + TAB)
        if r['author']:   
            fp.write(r['author'] + TAB)
        else: fp.write("" + TAB)
        ecolist=r['evidence_note'].split(";")
        eco=""
        for eco in ecolist:
            if "ECO:" in eco:break
        fp.write(r['evidence_note'] + TAB)
        fp.write(eco + TAB)
        if r['inferredfrom']:
            fp.write(r['inferredfrom'] + TAB)
        else: fp.write("" + TAB)
        fp.write("%d%s"%(r['stanza'],TAB))
        fp.write("%d%s"%(r['sequencenum'],CRT))

#
# Main
#
#Note: for now the generated report should be stored under the tr directory
#     When done testing, we will set outputdir=os.environ['REPORTOUTPUTDIR']
#
fp = reportlib.init(sys.argv[0], fileExt = '.rpt', outputdir = '/mgi/all/wts_projects/11900/11992/', printHeading = None)
log = reportlib.init(sys.argv[0], fileExt = '.log', outputdir = '/mgi/all/wts_projects/11900/11992/', printHeading = None)

i = datetime.now()
fp.write("Date:"+i.strftime('%Y/%m/%d %I:%M:%S \n\n'))
log.write("Date:"+i.strftime('%Y/%m/%d %I:%M:%S \n\n'))

#
#
#Index ECO-Evidence code with the associated annotation key and reference key
#
# GO GAF.evidenceCode: voc_term.term/abbreviation of voc_vocab._vocab_key=3 (GO Evidence Codes)
# GO GAF.annotation_properties: voc_term.term of voc_vocab._vocab_key=82 (GO Property)
# GO ECO code: voc_evidenceproperties.value of voc_term._term_key=6481772 
#              where voc_vocab._vocab_key=82 (EOCO code evidence)
#

db.sql('''select vp.value as ecoCodeNote ,vt2.abbreviation as evCode,
    ve._refs_key,ve._annot_key,ve.inferredfrom,vp.stanza, vp.sequencenum
    into #goEvProperties
    from voc_evidence_property  vp, voc_term vt2, voc_evidence ve,voc_term vt
    where vp._propertyterm_key=6481772
    and vp._propertyterm_key= vt._term_key
    and vt._vocab_key=82
    and vt.isObsolete=0
    and vp._annotevidence_key= ve._annotevidence_key
    and ve._evidenceterm_key=vt2._term_key
    and vt2._vocab_key=3
    and vt2.isObsolete=0
    ''', None)
db.sql('create index annot_idx1 on #goEvProperties(_annot_key)', None)

#
# Add GO term,GO id, and gaf.qualifier to the list
# GO GAF.Qualifier: voc_term.term of voc_vocab._vocab_key=52 (GO Qualifier)
# GO/Marker annotation type key _annottype_key= 1000 in voc_annot
#
db.sql('''
    select va._object_key as _marker_key,vp.ecoCodeNote ,
    vp.evCode,vp._refs_key,vp.inferredfrom, vp.stanza, vp.sequencenum,
    acc.accid as goid,vt2.term as go_term,vt.term as gaf_qualifier
    into #goEvQualifier
    from voc_annot va,#goEvProperties vp, voc_term vt2,voc_term vt,acc_accession acc
    where va._annot_key=vp._annot_key
    and va._annottype_key= 1000
    and va._qualifier_key=vt._term_key
    and vt._vocab_key=52
    and vt.isobsolete=0
    and va._term_key=vt2._term_key
    and vt2._vocab_key=4
    and vt2.isobsolete=0
    and va._term_key= acc._object_key
    and acc.preferred = 1
    and acc._logicaldb_key = 31
    and acc._mgitype_key = 13
    ''', None)
db.sql('create index marker_key_idx on #goEvQualifier(_marker_key)', None)
db.sql('create index ref_idx on #goEvQualifier(_refs_key)', None)

#
# Add GO Ref to the list
#
db.sql('''
    select g._marker_key,g.ecoCodeNote ,acc.accid as ref,acc2.accid as mgiref,
    b._primary as author,c.pubmedid,g.evCode,g._refs_key,g.inferredfrom, 
    g.stanza, g.sequencenum,g.goid,g.go_term,g.gaf_qualifier
    into #goEvRefs
    from #goEvQualifier g,acc_accession acc,acc_accession acc2,bib_refs b,bib_citation_cache c
    where g._refs_key=acc._object_key
    and acc._logicaldb_key = 1
    and acc._mgitype_key = 1
    and acc.prefixpart = 'J:'
    and g._refs_key=acc2._object_key
    and acc2._logicaldb_key = 1
    and acc2._mgitype_key = 1
    and acc2.prefixpart = 'MGI:'
    and  g._refs_key= b._refs_key
    and g._refs_key= c._refs_key
    ''', None)
db.sql('create index marker_idx1 on #goEvRefs(_marker_key)', None)

#
# Generate report for annotations with ECO code in the note
#
cmd ="""
   select acc.accid as mgiid, m.symbol,g.ecoCodeNote as evidence_note,
   g.ref,g.mgiref,g.author,g.pubmedid,g.evCode,g.inferredfrom,
   g.stanza,g.sequencenum,g.goid,g.go_term,g.gaf_qualifier
   from #goEvRefs g, acc_accession acc, mrk_marker m 
   where g.ecoCodeNote like '%ECO%' and g._marker_key=acc._object_key 
   and acc._logicaldb_key = 1
   and acc._mgitype_key = 2 and acc.prefixpart = 'MGI:' and acc.preferred=1
   and g._marker_key=m._marker_key and m._Marker_Status_key=1
"""

printResults(cmd)

i = datetime.now()
log.write("Date:"+i.strftime('%Y/%m/%d %I:%M:%S \n\n'))
print "\nProgram Complete\n"

reportlib.finish_nonps(fp)
reportlib.finish_nonps(log)

