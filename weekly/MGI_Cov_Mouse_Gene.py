'''
#
# MGI_Cov_Mouse_Gene.py
#
# Report:
#       TR13459
#       http://wts.informatics.jax.org/searches/tr.detail.cgi?TR_Nr=TR13459
#
# Output format:
#       1. gene symbol
#       2. gene name
#       3. feature type
#       4. mgi id
#
#       markers of type 'gene'
#       markers that have references with COV:xxx tags
#       markers that have alleles that have references with COV:xxx tags
#       exclude Gt(ROSA)26Sor/key = 37270
#
# History:
#
# lec   02/02/2021
#       - TR13459/Coronavirus Reports for Jax Web Resources
#
'''

import sys
import os
import reportlib
import db

db.setTrace()

TAB = reportlib.TAB
CRT = reportlib.CRT

db.useOneConnection(1)

fp = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])

fp.write('#\n')
fp.write('#1. gene symbol\n')
fp.write('#2. gene name\n')
fp.write('#3. feature type\n')
fp.write('#4. mgi id\n')
fp.write('#5. short citation\n')
fp.write('#6. jnum id\n')
fp.write('#7. covid tags\n')
fp.write('#\n')

cmd = '''
(
select distinct m._marker_key, m.symbol, m.name, a.accid, v.term, c._refs_key, c.jnumid, c.short_citation
into temp table covidrefs
from MRK_Marker m, ACC_Accession a, MRK_MCV_Cache v, 
        MGI_Reference_Assoc r, BIB_Citation_Cache c
where m._organism_key = 1
and m._marker_status_key in (1,3)
and m._marker_type_key = 1
and m._marker_key not in (37270)
and m._marker_key = a._object_key
and a._mgitype_key = 2
and a._logicaldb_key = 1
and a.prefixPart = 'MGI:'
and a.preferred = 1
and m._marker_key = v._marker_key
and v.qualifier = 'D'
and m._marker_key = r._object_key
and r._mgitype_key = 2
and r._refs_key = c._refs_key
and exists (select 1 from BIB_Workflow_Tag tg, VOC_Term tt
        where r._refs_key = tg._refs_key
        and tg._tag_key = tt._term_key
        and tt._vocab_key = 129
        and tt.term ilike 'COV:%'
        )
union
select distinct m._marker_key, m.symbol, m.name, a.accid, v.term, c._refs_key, c.jnumid, c.short_citation
from MRK_Marker m, ACC_Accession a, MRK_MCV_Cache v,
        ALL_Allele al, MGI_Reference_Assoc r, BIB_Citation_Cache c
where m._organism_key = 1
and m._marker_status_key in (1,3)
and m._marker_type_key = 1
and m._marker_key not in (37270)
and m._marker_key = a._object_key
and a._mgitype_key = 2
and a._logicaldb_key = 1
and a.prefixPart = 'MGI:'
and a.preferred = 1
and m._marker_key = v._marker_key
and v.qualifier = 'D'
and m._marker_key = al._marker_key
and al._allele_key = r._object_key
and r._mgitype_key = 11
and r._refs_key = c._refs_key
and exists (select 1 from BIB_Workflow_Tag tg, VOC_Term tt
        where r._refs_key = tg._refs_key
        and tg._tag_key = tt._term_key
        and tt._vocab_key = 129
        and tt.term ilike 'COV:%'
        )
)
'''
db.sql(cmd, None)
db.sql('create index covidKey on covidrefs(_refs_key)', None)

#
# covid tags lookup
#
covidtags = {}
cmd = '''
select distinct r._refs_key, tt.term as covidtag
from covidrefs r, bib_workflow_tag tg, voc_term tt
where r._refs_key = tg._refs_key
and tg._tag_key = tt._term_key
and tt._vocab_key = 129
and tt.term ilike 'COV:%'
'''
results = db.sql(cmd, 'auto')
for r in results:
        key = r['_refs_key']
        value = r['covidtag']
        if key not in covidtags:
                covidtags[key] = []
        covidtags[key].append(value)
#print(covidtags)

results = db.sql('select * from covidrefs order by symbol, short_citation', 'auto')

for r in results:
    refsKey = r['_refs_key']
    fp.write(r['symbol'] + TAB)
    fp.write(r['name'] + TAB)
    fp.write(r['term'] + TAB)
    fp.write(r['accid'] + TAB)
    fp.write(r['short_citation'] + TAB)
    fp.write(r['jnumid'] + TAB)
    fp.write('|'.join(covidtags[refsKey]) + CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
