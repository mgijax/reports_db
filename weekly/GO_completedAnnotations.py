#!/usr/local/bin/python

'''
#
# GO_completedAnnotations.py 01/30/2009
#
# Report:
#       Tab-delimited file of all markers with a completed go reference annotataion.
#
# Usage:
#       GO_completedAnnotations.py
#
# Used by:
#	Those who are interested in all the markers with complete go annotations.
#
# Output format:
#
#   1.  accID
#   2.  symbol
#   3.  date_complete
#
# History:
#
# mhall	01/30/2009
#	- new
#
'''

import sys
import os
import reportlib

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db


TAB = reportlib.TAB
CRT = reportlib.CRT

fp = reportlib.init('go_complete', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

cmd = '''
      select m.symbol, a.accID, 
      convert(char(10), t.completion_date, 112) as date_complete
      from GO_Tracking t, MRK_Marker m, ACC_Accession a 
      where t.completion_date is not null 
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
	fp.write(r['symbol'] + TAB)
	fp.write(str(r['date_complete']) + CRT)

reportlib.finish_nonps(fp)

