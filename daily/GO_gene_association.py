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
#   8.  Inferred From
#   9.  GO DAG Abbreviation (F, P, C)
#   10. Gene name
#   11. Gene synonym(s) - list of |-delimited synonyms
#   12. Marker Type (gene)
#   13. Species (taxon:10090)
#   14. Modification Date (YYYYMMDD)
#   15. Assigned By
#
# History:
#
# lec	01/25/2007
#	- TR 8122; don't convert inferredFrom delimiters; use as is in the database
#
# lec	09/14/2006
#	- TR 7904; append GOA annotations
#
# lec	10/19/2005
#	- added PMID, TR 7173
#
# lec	10/04/2005
#	- TR 5188; GO Qualifier
#
# lec	10/03/2005
#	- replace SWALL with UniProt
#
# lec	03/10/2004
#       - only include Markers of type Gene
#
# lec	02/11/2003
#	- TR 4511; add new column "assigned by"
#
# lec	04/02/2002
#	- TR 3527; some terms may not have DAGs; this shouldn't happen
#	but we don't want the report to bomb if it does.
#
# lec	01/28/2002
#	- new - revision 2.0
#
'''

import sys
import os
import string
import re
import db
import reportlib
import mgi_utils

DBABBREV = 'MGI'
SPECIES = 'taxon:10090'

TAB = reportlib.TAB
CRT = reportlib.CRT

db.useOneConnection(1)

fp = reportlib.init('gene_association', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

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

#
# retrieve all dag abbrevations for each term
#
results = db.sql('select distinct _Object_key, dagAbbrev = rtrim(dagAbbrev) from DAG_Node_View where _Vocab_key = 4', 'auto')
dag = {}
for r in results:
	dag[r['_Object_key']] = r['dagAbbrev']

#
# retrieve data set to process
#
db.sql('select a._Term_key, t.term, termID = ta.accID, qualifier = q.synonym, a._Object_key, ' + \
	'e.inferredFrom, e.modification_date, e._EvidenceTerm_key, e._Refs_key, e._ModifiedBy_key, ' + \
	'm.symbol, m.name, markerType = mt.name ' + \
	'into #gomarker ' + \
	'from VOC_Annot a, ACC_Accession ta, VOC_Term t, VOC_Evidence e, MRK_Marker m, MRK_Types mt, MGI_Synonym q ' + \
	'where a._AnnotType_key = 1000 ' + \
	'and a._Annot_key = e._Annot_key ' + \
	'and a._Object_key = m._Marker_key ' + \
	'and m._Marker_Type_key = 1 ' + \
	'and a._Term_key = t._Term_key ' + \
	'and a._Term_key = ta._Object_key ' + \
	'and ta._MGIType_key = 13 ' + \
	'and ta.preferred = 1 ' + \
	'and m._Marker_Type_key = mt._Marker_Type_key ' + \
	'and a._Qualifier_key = q._Object_key ' + \
	'and q._SynonymType_key = 1023', None)
db.sql('create index idx1 on #gomarker(_Object_key)', None)
db.sql('create index idx2 on #gomarker(_EvidenceTerm_key)', None)
db.sql('create index idx3 on #gomarker(_Refs_key)', None)
db.sql('create index idx4 on #gomarker(_ModifiedBy_key)', None)

#
# retrieve synonyms for markers in data set
#
results = db.sql('select distinct g._Object_key, s.synonym ' + \
	'from #gomarker g, MGI_Synonym s, MGI_SynonymType st ' + \
	'where g._Object_key = s._Object_key ' + \
	'and s._MGIType_key = 2 ' + \
	'and s._SynonymType_key = st._SynonymType_key ' + \
	'and st.synonymType = "exact" ' + \
	'order by g._Object_key', 'auto')
syns = {}
for r in results:
	key = r['_Object_key']
	value = r['synonym']
	if not syns.has_key(key):
		syns[key] = []
	syns[key].append(value)

#
# resolve foreign keys
#
db.sql('select g._Refs_key, g._Term_key, g.termID, g.qualifier, g.inferredFrom, ' + \
	'g._Object_key, g.symbol, g.name, g.markerType, ' + \
	'mDate = convert(varchar(10), g.modification_date, 112), ' + \
	'markerID = ma.accID, ' + \
	'refID = b.accID, ' + \
	'eCode = rtrim(t.abbreviation), ' + \
	'modifiedBy = u.login ' + \
	'into #results ' + \
	'from #gomarker g, ACC_Accession ma, ACC_Accession b, VOC_Term t, MGI_User u ' + \
	'where g._Object_key = ma._Object_key ' + \
	'and ma._MGIType_key = 2 ' + \
	'and ma.prefixPart = "MGI:" ' + \
	'and ma._LogicalDB_key = 1 ' + \
	'and ma.preferred = 1 ' + \
	'and g._Refs_key = b._Object_key ' + \
	'and b._MGIType_key = 1 ' + \
	'and b.prefixPart = "MGI:" ' + \
	'and b._LogicalDB_key = 1 ' + \
	'and g._EvidenceTerm_key = t._Term_key ' + \
	'and g._ModifiedBy_key = u._User_key', None)
db.sql('create index idx1 on #results(symbol)', None)
db.sql('create index idx2 on #results(_Refs_key)', None)

#
# resolve PubMed IDs for References
#
pubMed = {}
results = db.sql('select r._Refs_key, a.accID from #results r, ACC_Accession a ' + \
	'where r._Refs_key = a._Object_key ' + \
	'and a._MGIType_key = 1 ' + \
	'and a._LogicalDB_key = 29 ', 'auto')
for r in results:
    key = r['_Refs_key']
    value = r['accID']
    pubMed[key] = value

#
# resolve PubMed IDs for References
#
pubMed = {}
results = db.sql('select r._Refs_key, a.accID from #results r, ACC_Accession a ' + \
        'where r._Refs_key = a._Object_key ' + \
        'and a._MGIType_key = 1 ' + \
        'and a._LogicalDB_key = 29 ', 'auto')
for r in results:
    key = r['_Refs_key']
    value = r['accID']
    pubMed[key] = value

#
# process results
#
results = db.sql('select * from #results order by symbol', 'auto')
for r in results:
	# if we can't find the DAG for the Term, skip it

	if dag.has_key(r['_Term_key']):
		fp.write(DBABBREV + TAB)
		fp.write(r['markerID'] + TAB)
		fp.write(r['symbol'] + TAB)
		fp.write(string.strip(r['qualifier']) + TAB)
		fp.write(r['termID'] + TAB)

                # reference
                referenceID = DBABBREV + ':' + r['refID']
                if pubMed.has_key(r['_Refs_key']):
                    referenceID = referenceID + '|PMID:' + pubMed[r['_Refs_key']]
                fp.write(referenceID + TAB)
		fp.write(r['eCode'] + TAB)
		inferredFrom = re.sub('MGI:','MGI:MGI:',mgi_utils.prvalue(r['inferredFrom']))
		fp.write(inferredFrom + TAB)

		fp.write(dag[r['_Term_key']] + TAB)
		fp.write(r['name'] + TAB)

		if syns.has_key(r['_Object_key']):
			fp.write(string.join(syns[r['_Object_key']], '|'))

		fp.write(TAB)

		fp.write(r['markerType'] + TAB)
		fp.write(SPECIES + TAB)
		fp.write(r['mDate'] + TAB)

		if r['modifiedBy'] == 'swissload':
			fp.write('UniProt')
		elif string.find(r['modifiedBy'], 'GOA_') >= 0:
			modifiedBy = re.sub('GOA_', '', r['modifiedBy'])
			fp.write(modifiedBy)
		else:
			fp.write(DBABBREV)

		fp.write(CRT)
	
#
# append GOA annotations, if they exist
#

try:

    goaFile = open(os.environ['GOAMGI'], 'r')

    for line in goaFile.readlines():
        fp.write(line)
    goaFile.close()

except:

    pass

reportlib.finish_nonps(fp)
db.useOneConnection(0)
