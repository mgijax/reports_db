#!/usr/local/bin/python

'''
#
# ES_CellLine.py 
#
# Report:
#       Tab-delimited file
#       ES cell lines and their associated strains
#         - excluding 'Not Applicable','Not Specified','Other (see notes)'
#
# Usage:
#       ES_CellLine.py
#
# Notes:
#
# History:
#
# dbm    04/20/2007 - created (TR 8214)
#
'''

import sys
import os
import mgi_utils
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


#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# Get the cell lines and strains
#
results = db.sql('''
		select a.cellLine, p.strain 
                from ALL_CellLine a, PRB_Strain p 
                where a.isMutant = 0 
                      and a.cellLine not in ('Not Applicable','Not Specified','Other (see notes)') 
                      and a._Strain_key = p._Strain_key
                order by a.cellLine
		''', 'auto')

# Create the report
#
for r in results:
	fp.write(r['cellLine'] + reportlib.TAB + \
	       	 r['strain'] + reportlib.CRT)

reportlib.finish_nonps(fp)
