#!/usr/local/bin/python

'''
#
# TR 3655
#
# Report:
#       Tab-delimited file
#
#	1. mouse gene symbol
#	2. mgi acc id
#	3. gene symbol for human ortholog
#	4. locuslink acc id for human gene
#	5. refseq acc id for mouse gene
#	6. refseq acc id for human gene
#	7. sp acc id for mouse gene
#	8. sp acc id for human gene
#	9. if no ref seq or sp id, then GB acc ids for mouse gene
#	10. evidence used to support mouse0human homology
#
# Usage:
#       HMD_HumanSequence.py
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
# lec	05/07/2002
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

fp = reportlib.init(sys.argv[0], printHeading = 0, outputdir = os.environ['REPORTOUTPUTDIR'])

cmds = []

cmds.append('select distinct mouseKey = h1._Marker_key, mouseSym = m1.symbol, ' +
	'humanKey = h2._Marker_key, humanSym = m2.symbol, a.abbrev ' + \
	'into #homology ' +
        'from HMD_Homology r1, HMD_Homology_Marker h1, HMD_Homology_Assay ha, HMD_Assay a, ' + \
        'HMD_Homology r2, HMD_Homology_Marker h2, ' + \
        'MRK_Marker m1, MRK_Marker m2 ' + \
        'where m1._Species_key = 1 ' + \
        'and m1._Marker_key = h1._Marker_key ' + \
        'and h1._Homology_key = r1._Homology_key ' + \
        'and r1._Class_key = r2._Class_key ' + \
        'and r2._Homology_key = h2._Homology_key ' + \
        'and h2._Marker_key = m2._Marker_key ' + \
        'and m2._Species_key = 2 ' + \
	'and h1._Homology_key = ha._Homology_key ' + \
	'and ha._Assay_key = a._Assay_key')

cmds.append('create nonclustered index index_mouseKey on #homology(mouseKey)')
cmds.append('create nonclustered index index_humanKey on #homology(humanKey)')
cmds.append('create nonclustered index index_humanSym on #homology(humanSym)')

# MGI for Mouse
cmds.append('select a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.mouseKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1')

# LocusLink for Human
cmds.append('select distinct a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.humanKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 24 ')

# RefSeq for Mouse
cmds.append('select distinct a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.mouseKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 27 ')

# RefSeq for Human
cmds.append('select distinct _Object_key = h.humanKey, accID = r.refseqID ' + \
	'from #homology h, ACC_Accession a, tempdb..LLRef r ' + \
	'where h.humanKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 24 ' + \
	'and a.accID = r.locusID ')

# SWISSPROT for Mouse
cmds.append('select distinct a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.mouseKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 13 ')

# SWISSPROT for Human
cmds.append('select distinct _Object_key = h.humanKey, accID = r.proteinID ' + \
	'from #homology h, tempdb..LL l, tempdb..LLRef r ' + \
	'where h.humanSym = l.osymbol ' + \
	'and l.taxID = 9606 ' + \
	'and l.locusID = r.locusID')

# GenBank for Mouse
cmds.append('select distinct a._Object_key, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.mouseKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 9 ')

cmds.append('select distinct mouseKey, abbrev from #homology')

cmds.append('select distinct mouseKey, mouseSym, humanKey, humanSym from #homology')

results = db.sql(cmds, 'auto')

mgiID = {}
for r in results[4]:
	mgiID[r['_Object_key']] = r['accID']

llID = {}
for r in results[5]:
	llID[r['_Object_key']] = r['accID']

mrefseqID = {}
for r in results[6]:
	mrefseqID[r['_Object_key']] = r['accID']

hrefseqID = {}
for r in results[7]:
	hrefseqID[r['_Object_key']] = r['accID']

mspID = {}
for r in results[8]:
	if not mspID.has_key(r['_Object_key']):
		mspID[r['_Object_key']] = []
	mspID[r['_Object_key']].append(r['accID'])

hspID = {}
for r in results[9]:
	if r['accID'] != None:
		hspID[r['_Object_key']] = r['accID']

gbID = {}
for r in results[10]:
	if not gbID.has_key(r['_Object_key']):
		gbID[r['_Object_key']] = []
	gbID[r['_Object_key']].append(r['accID'])

habbrev = {}
for r in results[11]:
	if not habbrev.has_key(r['mouseKey']):
		habbrev[r['mouseKey']] = []
	habbrev[r['mouseKey']].append(r['abbrev'])

for r in results[12]:

	fp.write(r['mouseSym'] + TAB)
	fp.write(mgiID[r['mouseKey']] + TAB)
	fp.write(r['humanSym'] + TAB)

	if llID.has_key(r['humanKey']):
		fp.write(llID[r['humanKey']])
	fp.write(TAB)

	if mrefseqID.has_key(r['mouseKey']):
		fp.write(mrefseqID[r['mouseKey']])
	fp.write(TAB)

	if hrefseqID.has_key(r['humanKey']):
		fp.write(hrefseqID[r['humanKey']])
	fp.write(TAB)

	if mspID.has_key(r['mouseKey']):
		fp.write(string.join(mspID[r['mouseKey']], ','))
	fp.write(TAB)

	if hspID.has_key(r['humanKey']):
		fp.write(hspID[r['humanKey']])
	fp.write(TAB)

	if not mrefseqID.has_key(r['mouseKey']) and not mspID.has_key(r['mouseKey']) and gbID.has_key(r['mouseKey']):
		fp.write(string.join(gbID[r['mouseKey']], ','))
	fp.write(TAB)

	if habbrev.has_key(r['mouseKey']):
		fp.write(string.join(habbrev[r['mouseKey']], ','))
	fp.write(CRT)

reportlib.finish_nonps(fp)	# non-postscript file

