#!/usr/local/bin/python

'''
#
# GO_gene_association.py 01/28/2002
#
# Report:
#       Tab-delimited file of all MGI GO/Marker associations for Stanford
#
# Usage:
#       GO_gene_association.py
#
# Used by:
#	Stanford - central GO site
#	also made publically available on MGI FTP site
#
# Output format:
#
# The GO format has the following columns:
#
#   1.  Database designation (MGI)
#   2.  MGI Marker ID
#   3.  Symbol
#   4.  NOT
#   5.  GO id
#   6.  MGI ID of Reference (in format MGI:MGI:####) (double MGI: necessary)
#   7.  Evidence abbreviation
#   8.  GO Vocab Abbreviation (F, P, C)
#   9.  Gene name
#   10. Gene synonym(s) - list of |-delimited synonyms
#   11. Marker Type (gene, complex, etc.)
#   12. Species (10090)
#
# History:
#
# lec	01/28/2002
#	- new - revision 2.0
#
'''

import sys
import os
import string
import regsub
import db
import reportlib
import mgi_utils

DBABBREV = 'MGI'
SPECIES = 'taxon:10090'

# mapping between MGI Marker Type and what gets printed in the gene association file
markerTypes = {1: 'gene', 10: 'gene'}

TAB = reportlib.TAB
CRT = reportlib.CRT

fp = reportlib.init('gene_association.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

#
# Header information
# This file gets submitted to CVS so don't muck with this header
# Update the Revision number every time you modify this code 
# (like adding a new column)
#

fp.write('!software version: $Revision$\n')
fp.write('!date: %s $\n' % (mgi_utils.date("%m/%d/%Y")))
fp.write('!\n')
fp.write('! from Mouse Genome Database (MGD) & Gene Expression Database (GXD)\n')
fp.write('!\n')

cmds = []

# retreive all dag abbrevations for each term
cmds.append('select distinct _Object_key, dagAbbrev from DAG_Node_View where _Vocab_key = 4')

cmds.append('select a._Term_key, a.term, termID = a.accID, a.isNot, ' + \
	'm.symbol, m.name, m._Marker_key, markerID = ma.accID, m._Marker_Type_key, ' + \
	'e.inferredFrom, eCode = t.abbreviation, ' + \
	'mDate = convert(varchar(10), e.modification_date, 101), refID = b.accID ' + \
	'into #gomarker ' + \
	'from VOC_Annot_View a, MRK_Marker m, MRK_Acc_View ma, ' + \
	'VOC_Evidence e, VOC_Term t, BIB_Acc_View b ' + \
	'where a._AnnotType_key = 1000 ' + \
	'and a._Object_key = m._Marker_key ' + \
	'and m._Marker_key = ma._Object_key ' + \
	'and ma.prefixPart = "MGI:" ' + \
	'and ma.preferred = 1 ' + \
	'and a._Annot_key = e._Annot_key ' + \
	'and e._EvidenceTerm_key = t._Term_key ' + \
	'and e._Refs_key = b._Object_key ' + \
	'and b.prefixPart = "MGI:"')

cmds.append('select distinct m._Marker_key, o.name ' + \
	'from #gomarker m, MRK_Other o ' + \
	'where m._Marker_key = o._Marker_key ' + \
	'order by m._Marker_key')

cmds.append('select * from #gomarker order by symbol')

results = db.sql(cmds, 'auto')

# Get DAG abbreviations

dag = {}
for r in results[0]:

	dag[r['_Object_key']] = r['dagAbbrev']

# Get Marker Synonyms
syns = {}
for r in results[2]:
	
	if syns.has_key(r['_Marker_key']):
		syns[r['_Marker_key']].append(r['name'])
	else:
		syns[r['_Marker_key']] = []
		syns[r['_Marker_key']].append(r['name'])

for r in results[3]:

	fp.write(DBABBREV + TAB)
	fp.write(r['markerID'] + TAB)
	fp.write(r['symbol'] + TAB)

	if r['isNot'] == 1:
		fp.write('Not')

	fp.write(TAB)
	fp.write(r['termID'] + TAB)
	fp.write(DBABBREV + ':' + r['refID'] + TAB)
	fp.write(r['eCode'] + TAB)

	# substitute | for ", " in inferredFrom

	if r['inferredFrom'] != None:
		inferredFrom = regsub.gsub(', ', '|', r['inferredFrom'])
	else:
		inferredFrom = r['inferredFrom']

	fp.write(mgi_utils.prvalue(inferredFrom) + TAB)

	fp.write(dag[r['_Term_key']] + TAB)
	fp.write(r['name'] + TAB)

	if syns.has_key(r['_Marker_key']):
		fp.write(string.join(syns[r['_Marker_key']], '|'))

	fp.write(TAB)

	if markerTypes.has_key(r['_Marker_Type_key']):
		fp.write(markerTypes[r['_Marker_Type_key']] + TAB)
	else:
		fp.write('UNKNOWN' + TAB)

	fp.write(SPECIES + TAB)
	fp.write(r['mDate'] + CRT)
	
reportlib.finish_nonps(fp)

