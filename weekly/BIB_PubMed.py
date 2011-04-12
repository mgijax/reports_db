#!/usr/local/bin/python

'''
#
# BIB_PubMed.py 04/25/2002
#
# Report:
#       Tab-delimited file of all references and their pubmed ids
#
# Usage:
#       BIB_PubMed.py
#
# Used by:
#	Steve Grubb for Phenome project
#
# Output format:
#
#   1.  MGI ID for Reference
#   2.  PubMed ID for Reference
#
# doesn't include MGI IDs for References which don't have PubMed IDs
#
# History:
#
# lec	04/25/2002
#	- TR 3626
#
'''

import sys
import os
import db
import reportlib
import mgi_utils

TAB = reportlib.TAB
CRT = reportlib.CRT

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

db.sql('select a.accID, a._Object_key ' + \
	'into #refs ' + \
	'from ACC_Accession a ' + \
	'where a._MGIType_key = 1 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.preferred = 1', None)

db.sql('create unique index index_object_key on #refs(_Object_key)', None)

results = db.sql('select distinct r.accID, pubMedID = b.accID, jnum = a.accID ' + \
	'from #refs r, ACC_Accession b, ACC_Accession a ' + \
	'where r._Object_key = b._Object_key ' + \
	'and b._MGIType_key = 1 ' + \
	'and b._LogicalDB_key = 29 ' + \
	'and r._Object_key = a._Object_key ' + \
	'and a._MGIType_key = 1 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "J:" ' + \
	'and a.preferred = 1', 'auto')

for r in results:

	fp.write(r['accID'] + TAB)
	fp.write(r['pubMedID'] + TAB)
	fp.write(r['jnum'] + CRT)

reportlib.finish_nonps(fp)

