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
#
#       Report in a tab delimited file with the following columns:
#
#    		MGI:ID
#    		symbol
#    		Name
#    		GO ID
#		GO Term
#		PMID for GO assertion
#		SP ID (SP:)
#		GB ID of Human Marker (if it exists)
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

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

fp = reportlib.init("MRK_GOHuman.rpt", printHeading = 0, outputdir = os.environ['QCREPORTOUTPUTDIR'])

cmds = []

cmds.append('select m._Marker_key, m.symbol, m.name, ma.accID, a.term, goID = a.accID, e._Refs_key ' + \
'into #m1 ' + \
'from MRK_Marker m, MRK_Acc_View ma, VOC_Annot_View a, VOC_Evidence e ' + \
'where m._Species_key = 1 ' + \
'and m._Marker_Type_key = 1 ' + \
'and m._Marker_Status_key = 1 ' + \
'and m._Marker_key = ma._Object_key ' + \
'and ma.prefixPart = "MGI:" ' + \
'and ma.preferred = 1 ' + \
'and m._Marker_key = a._Object_key ' + \
'and a._AnnotType_key = 1000 ' + \
'and a._Annot_key = e._Annot_key ' + \
'and a.accID not in ("GO:0000004", "GO:0008372", "GO:0005554") ' + \
'and e._EvidenceTerm_key != 115 ' + \
'and e._Refs_key not in (61933, 73199) ')

cmds.append('create nonclustered index idx_marker_key on #m1(_Marker_key)')

cmds.append('select distinct m.*, humanMarker = m2._Marker_key ' + \
'into #m2 ' + \
'from #m1 m, HMD_Homology h1, HMD_Homology_Marker hm1, HMD_Homology h2, HMD_Homology_Marker hm2, MRK_Marker m2 ' + \
'where m._Marker_key = hm1._Marker_key ' + \
'and hm1._Homology_key = h1._Homology_key ' + \
'and h1._Class_key = h2._Class_key ' + \
'and h2._Homology_key = hm2._Homology_key ' + \
'and hm2._Marker_key = m2._Marker_key ' + \
'and m2._Species_key = 2 ')

# pubmedids for references
cmds.append('select distinct a._Object_key, a.accID ' + \
'from #m2 m, BIB_Acc_View a ' + \
'where m._Refs_key = a._Object_key ' + \
'and a._LogicalDB_key = 29')

# sp ids for mouse markers
cmds.append('select distinct a._Object_key, a.accID ' + \
'from #m2 m, MRK_Acc_View a ' + \
'where m._Marker_key = a._Object_key ' + \
'and a._LogicalDB_key = 13')

# genbank ids for human markers
cmds.append('select distinct a._Object_key, a.accID ' + \
'from #m2 m, MRK_Acc_View a ' + \
'where m.humanMarker = a._Object_key ' + \
'and a._LogicalDB_key = 9')

cmds.append('select distinct * from #m2 order by symbol')

results = db.sql(cmds, 'auto')

pmid = {}
for r in results[3]:
	pmid[r['_Object_key']] = 'PMID:' + r['accID']

spid = {}
for r in results[4]:
	if not spid.has_key(r['_Object_key']):
		spid[r['_Object_key']] = []

	spid[r['_Object_key']].append('SP:' + r['accID'])

gbid = {}
for r in results[5]:
	if not gbid.has_key(r['_Object_key']):
		gbid[r['_Object_key']] = []

	gbid[r['_Object_key']].append(r['accID'])

for r in results[6]:
	fp.write(r['accID'] + TAB + \
		r['symbol'] + TAB + \
		r['name'] + TAB + \
		r['goID'] + TAB + \
		r['term'] + TAB)

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

