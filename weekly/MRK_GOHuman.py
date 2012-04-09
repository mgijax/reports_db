#!/usr/local/bin/python

'''
#
# MRK_GOHuman.py 06/07/2002
#
# Report:
#       TR 3762
#
#	Select markers of type 'gene'
#               where the 'GO' association is not IEA
#               where the reference exludes J:60000 (61933), J:72447 (73199)
#		where the marker has a human ortholog
#		where the GO ID is not in (GO:0008150,GO:0003674,GO:0005575)
#		(these are the "unknown" terms within each DAG)
#
#       Report in a tab delimited file with the following columns:
#
#    		MGI:ID
#    		symbol
#    		Name
#    		GO ID
#		GO Term
#		PubMed ID for GO assertion (PMID:)
#		SwissProt ID (SP:)
#		GenBank ID of Human Marker (if it exists)
#
# Usage:
#       MRK_GOHuman.py
#
# History:
#
# lec	01/08/2002
#	- created
#
'''
 
import sys 
import os
import string
import db
import mgi_utils
import reportlib

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

fp = reportlib.init("MRK_GOHuman.rpt", printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])

# all mouse genes with a human ortholog

db.sql('select distinct h1._Marker_key, humanMarker = h2._Marker_key ' + \
	'into #ortholog ' + \
	'from MRK_Homology_Cache h1, MRK_Homology_Cache h2, MRK_Marker m ' + \
	'where h1._Organism_key = 1 ' + \
	'and h1._Marker_key = m._Marker_key ' + \
	'and m._Marker_Type_key = 1 ' + \
	'and m._Marker_Status_key in (1,3) ' + \
	'and h1._Class_key = h2._Class_key ' + \
	'and h2._Organism_key = 2 ', None)
db.sql('create index idx1 on #ortholog(_Marker_key)', None)
print 'query 1 end...%s' % (mgi_utils.date())

#
# select all orthologous mouse genes with annotation where evidence code != IEA (115)
# and the reference exludes J:60000 (61933), J:72447 (73199)
#

db.sql('select o._Marker_key, o.humanMarker, a._Annot_key, e._Refs_key, ' + \
	'e._EvidenceTerm_key, evidenceCode = t.abbreviation ' + \
	'into #temp1 ' + \
	'from #ortholog o, VOC_Annot a, VOC_Evidence e, VOC_Term t ' + \
	'where o._Marker_key = a._Object_key ' + \
	'and a._AnnotType_key = 1000 ' + \
	'and a._Annot_key = e._Annot_key ' + \
	'and e._EvidenceTerm_key != 115 ' + \
	'and e._Refs_key not in (61933, 73199) ' + \
	'and e._EvidenceTerm_key = t._Term_key', None)
db.sql('create index idx1 on #temp1(_Marker_key)', None)
db.sql('create index idx2 on #temp1(_Annot_key)', None)
print 'query 2 end...%s' % (mgi_utils.date())

# and the GO ID is not in (GO:0008150,GO:0003674,GO:0005575)

db.sql('select t.*, a.term, goID = a.accID ' + \
	'into #temp2 ' + \
	'from #temp1 t, VOC_Annot_View a ' + \
	'where t._Annot_key = a._Annot_key ' + \
	'and a.accID not in ("GO:0008150", "GO:0003674", "GO:0005575")', None)
db.sql('create index idx1 on #temp2(_Marker_key)', None)
print 'query 3 end...%s' % (mgi_utils.date())

# marker attributes

db.sql('select t.*, m.symbol, m.name ' + \
	'into #temp3 ' + \
	'from #temp2 t, MRK_Marker m ' + \
	'where t._Marker_key = m._Marker_key', None)
db.sql('create index idx1 on #temp3(_Marker_key)', None)
db.sql('create index idx2 on #temp3(_Refs_key)', None)
db.sql('create index idx3 on #temp3(humanMarker)', None)
print 'query 4 end...%s' % (mgi_utils.date())

# marker accession ids

results = db.sql('select a._Object_key, a.accID ' + \
	'from #temp3 t, ACC_Accession  a ' + \
	'where t._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.preferred = 1 ', 'auto')
markerid = {}
for r in results:
	key = r['_Object_key']
	value = r['accID']
	markerid[key] = value
print 'query 5 end...%s' % (mgi_utils.date())

# pubmedids for references

results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #temp3 t, ACC_Accession a ' + \
	'where t._Refs_key = a._Object_key ' + \
	'and a._MGIType_key = 1 ' + \
	'and a._LogicalDB_key = 29', 'auto')
pmid = {}
for r in results:
	pmid[r['_Object_key']] = 'PMID:' + r['accID']
print 'query 6 end...%s' % (mgi_utils.date())

# sp ids for mouse markers

results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #temp3 t, ACC_Accession a ' + \
	'where t._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 13', 'auto')
spid = {}
for r in results:
	if not spid.has_key(r['_Object_key']):
		spid[r['_Object_key']] = []
	spid[r['_Object_key']].append('SP:' + r['accID'])
print 'query 7 end...%s' % (mgi_utils.date())

# genbank ids for human markers
results = db.sql('select distinct a._Object_key, a.accID ' + \
	'from #temp3 t, ACC_Accession a ' + \
	'where t.humanMarker = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 9', 'auto')
gbid = {}
for r in results:
	if not gbid.has_key(r['_Object_key']):
		gbid[r['_Object_key']] = []
	gbid[r['_Object_key']].append(r['accID'])
print 'query 8 end...%s' % (mgi_utils.date())

# process records

results = db.sql('select distinct * from #temp3 order by symbol', 'auto')
print 'query 9 end...%s' % (mgi_utils.date())

# process each record in the final set
for r in results:
	fp.write(markerid[r['_Marker_key']] + TAB + \
		r['symbol'] + TAB + \
		r['name'] + TAB + \
		r['goID'] + TAB + \
		r['term'] + TAB + \
		r['evidenceCode'] + TAB)

	if pmid.has_key(r['_Refs_key']):
		fp.write(pmid[r['_Refs_key']])
	fp.write(TAB)

	if spid.has_key(r['_Marker_key']):
		fp.write(string.joinfields(spid[r['_Marker_key']], ';'))
	fp.write(TAB)

	if gbid.has_key(r['humanMarker']):
		fp.write(string.joinfields(gbid[r['humanMarker']], ';'))
	fp.write(CRT)

reportlib.finish_nonps(fp)	# non-postscript file
