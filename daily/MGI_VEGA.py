#!/usr/local/bin/python

'''
#
# TR 8055
#
# Report:
#       Produce a tab-delimited report with the following output fields:
#
#       MGD ID
#       Gene Symbol
#       Gene Name
#       Chromosome
#       cM Position
#       VEGA Gene ID
#
# Usage:
#       MRK_VEGA.py
#
# History:
#
# lec   12/13/2006
#       - created
#
'''
 
import sys 
import os
import db
import reportlib

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

cmd = 'select mgiID = a1.accID, m.symbol, m.name, m.chromosome, o.offset, vegaID = a2.accID ' + \
      'from ACC_Accession a1, ACC_Accession a2, MRK_Marker m, MRK_Offset o ' + \
      'where a1._Object_key = a2._Object_key and ' + \
            'a1._Object_key = m._Marker_key and ' + \
            'm._Marker_key = o._Marker_key and ' + \
            'a1._LogicalDB_key = 1 and ' + \
            'a1._MGIType_key = 2 and ' + \
            'a1.prefixPart = "MGI:" and ' + \
            'a1.preferred = 1 and ' + \
            'a2._LogicalDB_key = 85 and ' + \
            'a2._MGIType_key = 2 and ' + \
            'm._Marker_Type_key = 1 and ' + \
            'o.source = 0 ' + \
      'order by m.chromosome, m.symbol'

results = db.sql(cmd, 'auto')

for r in results:
    fp.write(r['mgiID'] + TAB + 
	     r['symbol'] + TAB + 
	     r['name'] + TAB +
	     str(r['offset']) + TAB +
             r['chromosome'] + TAB + 
             r['vegaID'] + CRT)

reportlib.finish_nonps(fp)
