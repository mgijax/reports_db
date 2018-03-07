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
import db

db.setTrace()

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# Get the cell lines and strains
#
results = db.sql('''
		select ac.cellLine, p.strain, a.accID
                from ALL_CellLine ac, PRB_Strain p, ACC_Accession a
                where ac.isMutant = 0 
                      and ac.cellLine not in ('Not Applicable','Not Specified','Other (see notes)') 
                      and ac._Strain_key = p._Strain_key
                      and p._Strain_key = a._Object_key
                      and a._MGIType_key = 10
                      and a._LogicalDB_key = 1
                      and a.preferred = 1
                order by ac.cellLine
		''', 'auto')

# Create the report
#
for r in results:
	fp.write(r['cellLine'] + reportlib.TAB + \
	       	 r['strain'] + reportlib.TAB + \
	       	 r['accID'] + reportlib.CRT)

reportlib.finish_nonps(fp)
