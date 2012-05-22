#!/usr/local/bin/python

'''
#
# MGI_iphone_app.py
#
# Report:
#       Tab-delimited file for iphone app
#	counts and accession ids by mouse marker
#	see TR11046 for more details
#
# Used by:
#	Jill Recla (jillrecla@jax.org), Carol Bult (carol.bult@jax.org)
#
# Output format:
#
#	1: MGI Marker ID
#	2: Symbol
#	3: Name
#	4: Start Coordinate
#	5: End Coordinate
#	6: Strand
#	7: # of References
#	8: MGI Ref ID (MGI:xxx|MGI:xxx|...)
#	9: # of Alleles
#	10: MGI Allele ID (MGI:xxx|MGI:xxx|...)
#	11: # of GO annotations (C group)
#	12: GO ID (C group): (GO:xxxx|GO:xxxx|...)
#	13: # of GO annotations (F group)
#	14: GO ID (F group): (GO:xxxx|GO:xxxx|...)
#	15: # of GO annotations (P group)
#	16: GO ID (P group): (GO:xxxx|GO:xxxx|...)
#	17: # of MP annotations
#	18: MP ID: (MP:xxxx|MP:xxxx|...)
#	19: # of OMIM annotations (via orthology)
#	20: OMIM ID: (xxxx|xxxx|...)
#	21: # of Nomenclature Events (via marker history)
#	22: Nomenclature sequence number:event:  
#			1||assigned|2||rename|3||split|4||deletion
# History:
#
# lec	05/17/2012
#	- TR11046; new
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

TAB = reportlib.TAB
CRT = reportlib.CRT

#
# report lives in repoort directory (not the FTP directory)
# archive lives in the FTP directory (reports/archive/...)
#

reportDir = os.environ['REPORTOUTPUTDIR']
archiveDir = os.environ['IPHONEARCHIVE']
currentDate = mgi_utils.date('%m-%d-%Y')

reportName = 'iphone_app_snapshot'
reportWithDate = '%s-%s' % (reportName, currentDate)
currentReport = '%s-current.rpt' % (reportName)

# remove current report link if it exists
if os.path.isfile('%s/%s' % (reportDir, currentReport)):
        os.remove('%s/%s' % (reportDir, currentReport))

# move existing reports to the archive
if os.path.isfile('%s/%s.rpt' % (reportDir, reportWithDate)):
	os.system('mv -f %s/%s.rpt %s' % (reportDir, reportWithDate, archiveDir))

fp = reportlib.init(reportWithDate, outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# for testing
#	and m.chromosome = '5'
#	and m.symbol = 'Kit'
#
db.sql('''
	select m._marker_key, m.symbol, m.name, a.accid
	into #markers
	from MRK_Marker m, ACC_Accession a
	where m._organism_key = 1
	and m._marker_status_key in (1,3)
	and m._marker_key = a._object_key
	and a._mgitype_key = 2
	and a._logicaldb_key = 1
	and a.prefixpart = 'mgi:'
	and a.preferred = 1
	''', None)
db.sql('create index marker_idx on #markers(_marker_key)', None)

#
# coordinates
#
results = db.sql('''    
    select m._marker_key,
           c.strand, 
           convert(int, c.startcoordinate) as startc,
           convert(int, c.endcoordinate) as endc
    from #markers m, MRK_Location_Cache c
    where m._marker_key = c._marker_key
        ''', 'auto')
coords = {}
for r in results:
    key = r['_marker_key']
    value = r
    if not coords.has_key(key):
        coords[key] = []
    coords[key].append(value)

#
# refs
#
results = db.sql('''
        select distinct m._marker_key, r.mgiid
        from #markers m, MRK_Reference r 
        where m._marker_key = r._marker_key
        ''', 'auto')
refs = {}
for r in results:
    key = r['_marker_key']
    value = r['mgiid']
    if not refs.has_key(key):
	refs[key] = []
    refs[key].append(value)

#
# alleles
#
results = db.sql('''
        select distinct m._marker_key, a.accid
        from #markers m, ALL_Marker_Assoc aa, ACC_Accession a
        where m._marker_key = aa._marker_key
	and aa._allele_key = a._object_key
	and a._mgitype_key = 11
	and a._logicaldb_key = 1
	and a.prefixpart = 'mgi:'
	and a.preferred = 1
        ''', 'auto')
alleles = {}
for r in results:
    key = r['_marker_key']
    value = r['accid']
    if not alleles.has_key(key):
	alleles[key] = []
    alleles[key].append(value)

#
# GO annotations
#
results = db.sql('''
        select distinct m._marker_key, a.accid, d.dag
        from #markers m, VOC_Annot aa, ACC_Accession a, VOC_Term t, DAG_Node_View d
        where m._marker_key = aa._object_key
	and aa._annottype_key = 1000
	and aa._term_key = a._object_key
	and a._mgitype_key = 13
	and a.preferred = 1
	and aa._Term_key = t._Term_key
	and t._Term_key = d._Object_key 
	and t._Vocab_key = d._Vocab_key 
        ''', 'auto')
goCannots = {}
goFannots = {}
goPannots = {}
for r in results:
    key = r['_marker_key']
    value = r['accid']

    if r['dag'] == 'Cellular Component':
	goannots = goCannots
    elif r['dag'] == 'Molecular Function':
	goannots = goFannots
    elif r['dag'] == 'Biological Process':
	goannots = goPannots

    if not goannots.has_key(key):
	goannots[key] = []
    goannots[key].append(value)

#
# Phenotype Annotations
#
results = db.sql('''
        select distinct m._marker_key, a.accid
        from #markers m, VOC_Annot aa, ACC_Accession a
        where m._marker_key = aa._object_key
	and aa._annottype_key = 1002
	and aa._term_key = a._object_key
	and a._mgitype_key = 13
	and a.preferred = 1
        ''', 'auto')
phenoannots = {}
for r in results:
    key = r['_marker_key']
    value = r['accid']

    if not phenoannots.has_key(key):
	phenoannots[key] = []
    phenoannots[key].append(value)

#
# OMIM ortholog (human)
#
results = db.sql('''
        select distinct m._marker_key, a.termid
        from #markers m, MRK_OMIM_Cache a
        where m._marker_key = a._marker_key
        ''', 'auto')
omimannots = {}
for r in results:
    key = r['_marker_key']
    value = r['termid']

    if not omimannots.has_key(key):
	omimannots[key] = []
    omimannots[key].append(value)

#
# Nomenclature History
#
results = db.sql('''
        select distinct m._marker_key, h.sequenceNum, e.event
        from #markers m, MRK_History h, MRK_Event e
        where m._marker_key = h._marker_key
	and h._marker_event_key = e._marker_event_key
        ''', 'auto')
nomen = {}
for r in results:
    key = r['_marker_key']
    value = r
    if not nomen.has_key(key):
	nomen[key] = []
    nomen[key].append(value)

#
# report
#
results = db.sql('select * from #markers', 'auto')

for r in results:

    key = r['_marker_key']

#	1: MGI Marker ID
#	2: Symbol
#	3: Name

    fp.write(r['accid'] + TAB)
    fp.write(r['symbol'] + TAB)
    fp.write(r['name'] + TAB)

#	4: Start Coordinate
#	5: End Coordinate
#	6: Strand

    if coords.has_key(key):
        fp.write(mgi_utils.prvalue(coords[key][0]['startc']) + TAB)
        fp.write(mgi_utils.prvalue(coords[key][0]['endc']) + TAB)
        fp.write(mgi_utils.prvalue(coords[key][0]['strand']) + TAB)
    else:
        fp.write(TAB + TAB + TAB)

#	7: # of References
#	8: MGI Ref ID (MGI:xxx|MGI:xxx|...)
	
    if refs.has_key(key):
	i = refs[key]
	fp.write(str(len(refs[key])) + TAB)
	fp.write(string.join(refs[key], '|') + TAB)
    else:
        fp.write('0' + TAB + TAB)

#	9: # of Alleles
#	10: MGI Allele ID (MGI:xxx|MGI:xxx|...)

    if alleles.has_key(key):
	fp.write(str(len(alleles[key])) + TAB)
	fp.write(string.join(alleles[key], '|') + TAB)
    else:
        fp.write('0' + TAB + TAB)

#	11: # of GO annotations (C group)
#	12: GO ID (C group): (GO:xxxx|GO:xxxx|...)

    if goCannots.has_key(key):
	fp.write(str(len(goCannots[key])) + TAB)
	fp.write(string.join(goCannots[key], '|') + TAB)
    else:
        fp.write('0' + TAB + TAB)

#	13: # of GO annotations (F group)
#	14: GO ID (F group): (GO:xxxx|GO:xxxx|...)

    if goFannots.has_key(key):
	fp.write(str(len(goFannots[key])) + TAB)
	fp.write(string.join(goFannots[key], '|') + TAB)
    else:
        fp.write('0' + TAB + TAB)

#	15: # of GO annotations (P group)
#	16: GO ID (P group): (GO:xxxx|GO:xxxx|...)

    if goPannots.has_key(key):
	fp.write(str(len(goPannots[key])) + TAB)
	fp.write(string.join(goPannots[key], '|') + TAB)
    else:
        fp.write('0' + TAB + TAB)

#	17: # of MP annotations
#	18: MP ID: (MP:xxxx|MP:xxxx|...)

    if phenoannots.has_key(key):
	fp.write(str(len(phenoannots[key])) + TAB)
	fp.write(string.join(phenoannots[key], '|') + TAB)
    else:
        fp.write('0' + TAB + TAB)

#	19: # of OMIM annotations (via orthology)
#	20: OMIM ID: (xxxx|xxxx|...)

    if omimannots.has_key(key):
	fp.write(str(len(omimannots[key])) + TAB)
	fp.write(string.join(omimannots[key], '|') + TAB)
    else:
        fp.write('0' + TAB + TAB)

#	21: # of Nomenclature Events (via marker history)
#	22: Nomenclature sequence number:event:  

    if nomen.has_key(key):
	for n in nomen[key]:
	    fp.write(str(n['sequenceNum']) + '||' + n['event'] + '|')
        fp.write(CRT)
    else:
        fp.write('0' + TAB + CRT)

reportlib.finish_nonps(fp)

# re-create a symbolic link between the new file and the current file
os.chdir(reportDir)
os.symlink(reportWithDate + '.rpt', currentReport)

