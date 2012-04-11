#!/usr/local/bin/python

'''
#
# PRB_CloneLibrary.py
#
# Report:
#       Tab-delimited file of Clone Libraires
#
# Usage:
#       PRB_CloneLibrary.py
#
# History:
#
# lec	04/11/2012
#	- TR 11035/move from sql format to tab-delimited format
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

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

results = db.sql('''
	select substring(name, 1, 100) as name
	from PRB_Source
	where name is not null
	order by name
	''', 'auto')

for r in results:
	fp.write(r['name'] + CRT)

reportlib.finish_nonps(fp)

