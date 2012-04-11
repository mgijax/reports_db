#!/usr/local/bin/python

'''
#
# PRB_TissueMouse.py
#
# Report:
#       Tab-delimited file of Probe/Tissues
#
# Usage:
#       PRB_TissueMouse.py
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

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None, \
	title = 'Molecular Segments - All Mouse Tissues')

db.sql('''
	select distinct _Tissue_key
	into #tissues
	from PRB_Probe p, PRB_Source s, VOC_Term t
	where p._SegmentType_key = t._Term_key
	and t.term = 'cDNA'
	and p._Source_key = s._Source_key
	and s._Organism_key = 1
	''', None)

results = db.sql('''
	select p.tissue 
	from #tissues t, PRB_Tissue p
	where t._Tissue_key = p._Tissue_key
	order by tissue
	''', 'auto')

for r in results:
	fp.write(r['tissue'] + CRT)

reportlib.finish_nonps(fp)

