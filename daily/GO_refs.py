#!/usr/local/bin/python

'''
#
# GO_refs.py 03/21/2002
#
# Report:
#       Tab-delimited file of all GO references for NCBI (LocusLink)
#
# Usage:
#       GO_refs.py
#
# Used by:
#	LocusLink for providing GO/Reference links
#	publically available on MGI FTP site
#
# Output format:
#
#   1.  MGI ID for GO Reference
#   2.  PubMed ID for Reference (blank if we don't have the PubMed ID)
#
# History:
#
# lec	03/21/2002
#	- new
#
'''

import sys
import os
import db
import reportlib
import mgi_utils

TAB = reportlib.TAB
CRT = reportlib.CRT

fp = reportlib.init('go_refs', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []

cmds.append('select distinct e._Refs_key, b.accID ' + \
	'into #gorefs ' + \
	'from VOC_Annot a, VOC_Evidence e, BIB_Acc_View b ' + \
	'where a._AnnotType_key = 1000 ' + \
	'and a._Annot_key = e._Annot_key ' + \
	'and e._Refs_key = b._Object_key ' + \
	'and b.prefixPart = "MGI:"')

cmds.append('select g.*, pubMedID = b.accID ' + \
	'from #gorefs g, BIB_Acc_View b ' + \
	'where g._Refs_key = b._Object_key ' + \
	'and b._LogicalDB_key = 29 ' + \
        'union ' + \
        'select g.*, null ' + \
	'from #gorefs g ' + \
	'where not exists (select 1 from BIB_Acc_View b ' + \
	'where g._Refs_key = b._Object_key ' + \
	'and b._LogicalDB_key = 29)')

results = db.sql(cmds, 'auto')

for r in results[-1]:

	fp.write(r['accID'] + TAB)
	fp.write(mgi_utils.prvalue(r['pubMedID']) + CRT)

reportlib.finish_nonps(fp)

