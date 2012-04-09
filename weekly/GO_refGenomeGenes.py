#!/usr/local/bin/python

'''
#
# GO_refGenomeGenes.py 01/30/2009
#
# Report:
#       Tab-delimited file of all GO Genome Reference Genes
#
# Usage:
#       GO_refGenomeGenes.py
#
# Used by:
#	Those who are interested in a list of all the Genome Reference Genes.
#
# Output format:
#
#   1.  accID
#   2.  symbol
#
# History:
#
# mhall	01/30/2009
#	- new
#
'''

import sys
import os
import db
import reportlib

TAB = reportlib.TAB
CRT = reportlib.CRT

fp = reportlib.init('go_ref_genome', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

cmd = '''
      select m.symbol, a.accID 
      from GO_Tracking t, MRK_Marker m, ACC_Accession a 
      where t.isReferenceGene = 1 
      and t._Marker_key = m._Marker_key 
      and m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 1 
      and a.prefixPart = "MGI:" 
      and a.preferred = 1 
      order by m.symbol
      '''

results = db.sql(cmd, 'auto')

for r in results:
	fp.write(r['accID'] + TAB)
	fp.write(r['symbol'] + CRT)

reportlib.finish_nonps(fp)

