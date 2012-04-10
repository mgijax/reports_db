#!/usr/local/bin/python

'''
#
# MRK_Dump2.py 11/16/98
#
# Report:
#       Tab-delimited file of MGI Mouse Markers
#       excluding Withdrawns Symbols
#
# Usage:
#       MRK_Dump2.py
#
# History:
#
# 07/20/2011	lec
#	- ansi standard; change str(offst) to just offset
#
# 05/30/2002	lec
#	- TR 3736; add Marker Type
#
# 01/13/98	lec
#	- added comments
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


fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

db.sql('''
	select m._Marker_key, m.symbol, m.name, m.chromosome, t.name as markerType, offset 
	into #markers 
	from MRK_Marker m, MRK_Types t, MRK_Offset o 
	where m._Organism_key = 1 
	and m._Marker_Status_key in (1,3) 
	and m._Marker_Type_key = t._Marker_Type_key 
	and m._Marker_key = o._Marker_key 
	and o.source = 0
	''', None)
db.sql('create index idx1 on #markers(_Marker_key)', None)
db.sql('create index idx2 on #markers(symbol)', None)

results = db.sql('''
	select m.symbol, m.symbol, m.name, m.chromosome, m.markerType, m.offset, a.accID 
	from #markers m, ACC_Accession a 
	where m._Marker_key = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 1 
	and a.prefixPart = "MGI:" 
	and a.preferred = 1 
	order by m.symbol
	''', 'auto')

for r in results:
	fp.write(r['accID'] + reportlib.TAB + \
	         r['symbol'] + reportlib.TAB + \
		 r['name'] + reportlib.TAB + \
		 str(r['offset']) + reportlib.TAB + \
		 r['chromosome'] + reportlib.TAB + \
		 r['markerType'] + reportlib.CRT)

reportlib.finish_nonps(fp)

