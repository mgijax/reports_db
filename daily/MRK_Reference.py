#!/usr/local/bin/python

'''
#
# MRK_Reference.py 01/23/2003
#
# Report:
#       TR 4454
#	Tab-delimited file of:
#		Mouse Marker MGI ID
#		Symbol
#		Name
#		Synonyms
#		PubMed IDs for References 
#	Sorted by Symbol
#
# Usage:
#       MRK_Reference.py
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
# lec	01/23/2003
#	- created
#
'''
 
import sys 
import os
import db
import reportlib
import string

TAB = reportlib.TAB
CRT = reportlib.CRT

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []

cmds.append('select a.accID, m._Marker_key, m.symbol, m.name ' + \
	'into #markers ' + \
	'from MRK_Marker m, MRK_Acc_View a ' + \
	'where m._Species_key = 1 ' + \
	'and m._Marker_Status_key = 1 ' + \
	'and m._Marker_key = a._Object_key ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1')

cmds.append('create nonclustered index idx_marker on #markers(_Marker_key)')

cmds.append('select distinct m._Marker_key, a.accID ' + \
	'from #markers m, MRK_Reference r, BIB_Acc_View a ' + \
	'where m._Marker_key = r._Marker_key ' + \
	'and r._Refs_key = a._Object_key ' + \
	'and a._LogicalDB_key = 29')

cmds.append('select m._Marker_key, o.name ' + \
	'from #markers m, MRK_Other o ' + \
	'where m._Marker_key = o._Marker_key')

cmds.append('select * from #markers order by symbol')

results = db.sql(cmds, 'auto')

pubmed = {}
for r in results[-3]:
	key = r['_Marker_key']
	if not pubmed.has_key(key):
		pubmed[key] = []

	pubmed[key].append(r['accID'])

syn = {}
for r in results[-2]:
	key = r['_Marker_key']
	if not syn.has_key(key):
		syn[key] = []

	syn[key].append(r['name'])

for r in results[-1]:

	# The list should include only publications with PubMed identifiers. (per TR)

	if pubmed.has_key(r['_Marker_key']):

		fp.write(r['accID'] + TAB + \
	         	r['symbol'] + TAB + \
		 	r['name'] + TAB)

		if syn.has_key(r['_Marker_key']):
			fp.write(string.joinfields(syn[r['_Marker_key']], '|'))

		fp.write(TAB)
		fp.write(string.joinfields(pubmed[r['_Marker_key']], '|') + CRT)

reportlib.finish_nonps(fp)	# non-postscript file

