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

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []

cmds.append('select a.accID, a._Object_key ' + \
	'into #refs ' + \
	'from BIB_Acc_View a ' + \
	'where a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a._LogicalDB_key = 1')

cmds.append('create unique index index_object_key on #refs(_Object_key)')

cmds.append('select r.accID, pubMedID = b.accID ' + \
	'from #refs r, BIB_Acc_View b ' + \
	'where r._Object_key = b._Object_key ' + \
	'and b._LogicalDB_key = 29 ')

results = db.sql(cmds, 'auto')

for r in results[-1]:

	fp.write(r['accID'] + TAB)
	fp.write(mgi_utils.prvalue(r['pubMedID']) + CRT)

reportlib.finish_nonps(fp)

