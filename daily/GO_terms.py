#!/usr/local/bin/python

'''
#
# GO_terms.py 01/30/2002
#
# Report:
#       Tab-delimited file of all GO Terms which have been annotated to MGI Markers
#
# Usage:
#       GO_terms.py
#
# Used by:
#	Those who are interested in the GO ID + term (which are not available
#	in the gene_association.mgi file)
#
# Output format:
#
#   1.  GO Vocab Name
#   2.  GO ID
#   3.  GO Term
#
# History:
#
# lec	01/30/2002
#	- new
#
'''

import sys
import os
import db
import reportlib

TAB = reportlib.TAB
CRT = reportlib.CRT

fp = reportlib.init('go_terms', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmd = 'select t.term, t.accID, d.dag ' + \
	'from VOC_Term_View t, DAG_Node_View d ' + \
	'where t._Vocab_key = 4 ' + \
	'and t._Vocab_key = d._Vocab_key ' + \
	'and t._Term_key = d._Object_key ' + \
	'and exists (select 1 from VOC_Annot a ' + \
	'where t._Term_key = a._Term_key) ' + \
	'order by dag, t.accID'

results = db.sql(cmd, 'auto')

for r in results:

	fp.write(r['dag'] + TAB)
	fp.write(r['accID'] + TAB)
	fp.write(r['term'] + CRT)

reportlib.finish_nonps(fp)

