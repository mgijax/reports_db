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
#		where the GO ID is not in (GO:0000004,GO:0008372,GO:0005554)
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
# Notes:
#	- all reports use mgireport directory for output file
#	- all reports use db default of public login
#	- all reports use server/database default of environment
#	- use lowercase for all SQL commands (i.e. select not SELECT)
#	- all public SQL reports require the header and footer
#	- all private SQL reports require the header
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
import reportlib
import mgi_utils

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

db.useOneConnection(1)

fp = reportlib.init("MRK_GOHuman.rpt", printHeading = 0, outputdir = os.environ['REPORTOUTPUTDIR'])

# all mouse genes with a human ortholog

cmds = []
cmds.append('select distinct m._Marker_key, humanMarker = m2._Marker_key ' + \
	'into #ortholog ' + \
	'from MRK_Marker  m, HMD_Homology h1, HMD_Homology_Marker hm1, ' + \
	'HMD_Homology h2, HMD_Homology_Marker hm2, MRK_Marker m2 ' + \
	'where m._Organism_key = 1 ' + \
	'and m._Marker_Type_key = 1 ' + \
	'and m._Marker_Status_key in (1,3) ' + \
	'and m._Marker_key = hm1._Marker_key ' + \
	'and hm1._Homology_key = h1._Homology_key ' + \
	'and h1._Class_key = h2._Class_key ' + \
	'and h2._Homology_key = hm2._Homology_key ' + \
	'and hm2._Marker_key = m2._Marker_key ' + \
	'and m2._Organism_key = 2 ')
cmds.append('create index idx1 on #ortholog(_Marker_key)')
db.sql(cmds, None)
print 'query 1 end...%s' % (mgi_utils.date())

#
# select all orthologous mouse genes with annotation where evidence code != IEA (115)
# and the reference exludes J:60000 (61933), J:72447 (73199)
#

cmds = []
cmds.append('select o._Marker_key, o.humanMarker, a._Annot_key, e._Refs_key, ' + \
	'e._EvidenceTerm_key, evidenceCode = t.abbreviation ' + \
	'into #temp1 ' + \
	'from #ortholog o, VOC_Annot a, VOC_Evidence e, VOC_Term t ' + \
	'where o._Marker_key = a._Object_key ' + \
	'and a._AnnotType_key = 1000 ' + \
	'and a._Annot_key = e._Annot_key ' + \
	'and e._EvidenceTerm_key != 115 ' + \
	'and e._Refs_key not in (61933, 73199) ' + \
	'and e._EvidenceTerm_key = t._Term_key')
cmds.append('create index idx1 on #temp1(_Marker_key)')
cmds.append('create index idx2 on #temp1(_Annot_key)')
db.sql(cmds, None)
print 'query 2 end...%s' % (mgi_utils.date())

# and the GO ID is not in (GO:0000004,GO:0008372,GO:0005554)

cmds = []
cmds.append('select t.*, a.term, goID = a.accID ' + \
	'into #temp2 ' + \
	'from #temp1 t, VOC_Annot_View a ' + \
	'where t._Annot_key = a._Annot_key ' + \
	'and a.accID not in ("GO:0000004", "GO:0008372", "GO:0005554")')
cmds.append('create index idx1 on #temp2(_Marker_key)')
db.sql(cmds, None)
print 'query 3 end...%s' % (mgi_utils.date())

# marker attributes

cmds = []
cmds.append('select t.*, m.symbol, m.name ' + \
	'into #temp3 ' + \
	'from #temp2 t, MRK_Marker m ' + \
	'where t._Marker_key = m._Marker_key')
cmds.append('create index idx1 on #temp3(_Marker_key)')
cmds.append('create index idx2 on #temp3(_Refs_key)')
cmds.append('create index idx3 on #temp3(humanMarker)')
db.sql(cmds, None)
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

db.useOneConnection(0)
