#!/usr/local/bin/python

'''
#
#
# TR9308: OMIM report for Jeremy Miller
#
# Associations of Mouse Genes with OMIM Phenotypes
# 
# Report:
#
# 1) Human Disease (Term, Symbol, and ID)
# 2) Synonyms
# 3) Associated Genes 
#    (Mouse gene; Human gene; Characteristics of this human
#     disease are associated with mutations in... : ...the mouse gene, ...the
#     human gene, and ...both mouse and human orthologous genes.).
#
# ==> For this header, I don't need any of the links and I don't care how the
# three "Characteristics of this human disease..." columns are presented. You
# can give me one column with the words "Mouse only", "Human only", or "Both;"
# or you can give me three columns with the word "Yes" in the proper column;
# or you can give me two columns corresponding to mouse and human with a "yes"
# or a "no" in each column.
#
# 4) If there is a column that lists the type of human disease that the entry
# corresponds to (such as "neurological disorder," "muscle disease," etc.)
# then I would like a column with that information, but if not, it is not
# necessary.
#
# I have gotten all of these column names from this website:
# http://www.informatics.jax.org/userdocs/omim_detail_help.shtml. Ideally, I
# would like this information for all of the entries--I can filter the data
# myself if necessary, so don't worry about that. Please feel free to contact
# me if any of the above information is not clear. Thanks again!
# Sincerely,
# Jeremy Miller
# jeremymiller@ucla.edu
# 
# History:
#
# 07/05/2011	lec
#	- TR10770/fixed bugs
#
# 05/12/2011	lec
#	- TR9308/add MGI ids/add to public reports
#
# 2008-11-06	Susan McClatchy
#	- created
'''

import sys
import os
import string
import reportlib

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
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

db.useOneConnection(1)

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('# OMIM ID' + CRT)
fp.write('# OMIM Term' + CRT)
fp.write('# Synonyms (|-delimited)' + CRT)
fp.write('# Marker Symbol' + CRT)
fp.write('# MGI Marker Accession ID' + CRT)
fp.write('# Human Marker Symbol' + CRT)
fp.write('# "Mouse only" or "Human only" or "Both"' + CRT*2)

#
# all organisms
#
db.sql('''
       select distinct o._Organism_key, o._Marker_key, o._Genotype_key, o._Term_key, 
       o._OrthologOrganism_key, o._OrthologMarker_key, o.qualifier, o.markerSymbol, 
       o.termID, o.term, o.orthologSymbol
       into #omimall
       from MRK_OMIM_Cache o
       ''', None)

#
# select mouse orthologs with human ortholog
# type = 1
#

db.sql('''
        select distinct m._Organism_key, m._Marker_key, m._Term_key, 
        m._OrthologOrganism_key, m._OrthologMarker_key, m.markerSymbol, 
        m.termID, m.term, m.orthologSymbol, type = 1
        into #omim
        from #omimall o, MRK_OMIM_Cache m 
        where o._Organism_key = 2 
        and o._OrthologMarker_key is not null 
        and o._OrthologOrganism_key = 1 
        and o._Term_key = m._Term_key 
        and o._Marker_key = m._OrthologMarker_key 
        and m._OrthologOrganism_key = 2 
        and o._OrthologMarker_key = m._Marker_key 
        and m._Organism_key = 1 
        and m._Genotype_key is not null 
        and m.qualifier is null

        ''', None)

#
# select mouse orthologs only
# type = 2
#

db.sql('''
	insert into #omim
	select distinct o._Organism_key, o._Marker_key, o._Term_key, 
        o._OrthologOrganism_key, o._OrthologMarker_key, o.markerSymbol, 
        o.termID, o.term, o.orthologSymbol, type = 2
	from #omimall o 
	where o._Organism_key = 1 
	and o._Genotype_key is not null 
	and o.qualifier is null 
	and not exists 
	  (select _Term_key from MRK_OMIM_Cache m 
	  where o._Term_key = m._Term_key 
	  and o._OrthologMarker_key = m._Marker_key 
	  and m._Organism_key = 2)
	''', None)

#
# select human orthologs only
# type = 3
#

db.sql('''
	insert into #omim
	select distinct o._Organism_key, o._Marker_key, o._Term_key, 
        o._OrthologOrganism_key, o._OrthologMarker_key, o.markerSymbol, 
        o.termID, o.term, o.orthologSymbol, type = 3
	from #omimall o
	where o._Organism_key = 2 
	and not exists 
	  (select m._Marker_key from MRK_OMIM_Cache m 
	  where o._Term_key = m._Term_key 
	  and o._Marker_key = m._OrthologMarker_key 
	  and m._Organism_key = 1 
	  and m._Genotype_key is not null 
	  and m.qualifier is null)
	''', None)

#
# index
#
db.sql('create index idx1 on #omim(_Marker_key)', None)
db.sql('create index idx2 on #omim(_Term_key)', None)

#
# MGI ids
#

results = db.sql('''
	select distinct o._Marker_key, a.accID 
        from #omim o, ACC_Accession a 
        where o._Organism_key = 1
	and o._Marker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a._LogicalDB_key = 1 
        and a.prefixPart = 'MGI:' 
        and a.preferred = 1
	union
	select distinct o._Marker_key, a.accID 
        from #omim o, ACC_Accession a 
        where o._OrthologOrganism_key = 1
	and o._OrthologMarker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a._LogicalDB_key = 1 
        and a.prefixPart = 'MGI:' 
        and a.preferred = 1
	''', 'auto')
mgiID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    mgiID[key] = value

#
# synonyms of terms
#

syns = {}
results = db.sql('''
	select distinct o._Term_key, s.synonym 
	from #omim o, MGI_Synonym s
	where o._Term_key = s._Object_key
	''', 'auto')
for r in results:
	key = r['_Term_key']
	value = r['synonym']
	if not syns.has_key(key):
		syns[key] = []
	syns[key].append(value)

#
# to print
#

results = db.sql('select * from #omim order by type, markerSymbol', 'auto')

for r in results:

	termKey = r['_Term_key']
	markerKey = r['_Marker_key']
	orthologSymbol = r['orthologSymbol']
	thisType = r['type']

	if orthologSymbol == None:
	    orthologSymbol = ''

	if mgiID.has_key(r['_Marker_key']):
	    mid = mgiID[r['_Marker_key']]
        else:
	    mid = ''

	fp.write(r['termID'] + TAB)
	fp.write(r['term'] + TAB)

	if syns.has_key(termKey):
	    fp.write(string.join(syns[termKey], '|'))
	fp.write(TAB)

	if thisType == 1:
	    fp.write(r['markerSymbol'] + TAB)
	    fp.write(mid + TAB)
	    fp.write(orthologSymbol + TAB)
	    fp.write('Both' + CRT)

	elif thisType == 2:
	    fp.write(r['markerSymbol'] + TAB)
	    fp.write(mid + TAB)
	    fp.write(orthologSymbol + TAB)
	    fp.write('Mouse only' + CRT)

	else:
	    fp.write(orthologSymbol + TAB)
	    fp.write(mid + TAB)
	    fp.write(r['markerSymbol'] + TAB)
	    fp.write('Human only' + CRT)

db.useOneConnection(0)
reportlib.finish_nonps(fp)

