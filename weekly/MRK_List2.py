#!/usr/local/bin/python

'''
#
# MRK_List2.py 07/01/2008
#
# Report:
#       TR 9120
#	Tab-delimited version of existing report, MRK_List2.sql
#	"Genetic Marker List (sorted alphabetically/excludes withdrawns)"
#
# Usage:
#       MRK_List2.py
#
# History:
#
# jer	07/01/2008
#	- created
#
'''
 
import sys 
import os
import string
import db
import reportlib

TAB = reportlib.TAB
CRT = reportlib.CRT

#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

headers = [
    "MGI Accession ID",
    "Chr", 
    "cM Position", 
    "Symbol", 
    "Status", 
    "Name", 
    "Type"]

fp.write(TAB.join(headers) + CRT)

cmd = '''
    select _Marker_key, _Marker_Status_key, _Marker_Type_key, 
	   symbol, name = substring(name,1,150), chromosome
    into #markers
    from MRK_Marker
    where _Organism_key = 1
    and _Marker_Status_key in (1, 3)
    '''
db.sql(cmd, 'auto')

cmd = '''
    create index idx1 on #markers(_Marker_key)
    '''
db.sql(cmd, 'auto')

cmd = '''
    select m.*, 
    markerStatus = upper(substring(s.status, 1, 1)),
    markerType = substring(t.name,1,25),
    cmPosition =
        case
        when o.offset >= 0 then str(o.offset, 10, 2)
        when o.offset = -999.0 then "       N/A"
        when o.offset = -1.0 then "  syntenic"
        end
    into #markersAll
    from #markers m, MRK_Status s, MRK_Types t, MRK_Offset o
    where m._Marker_key = o._Marker_key
    and o.source = 0
    and m._Marker_Type_key = t._Marker_Type_key
    and m._Marker_Status_key = s._Marker_Status_key
    '''
db.sql(cmd, 'auto')

cmd = '''
    create index idx1 on #markersAll(_Marker_key)
    '''
db.sql(cmd, 'auto')

cmd = '''
    create index idx2 on #markersAll(symbol)
    '''
db.sql(cmd, 'auto')


cmd = '''
    select a.accID "MGI Accession ID", m.chromosome "Chr", 
	m.cmPosition "cM Position", 
	m.symbol "Symbol", m.markerStatus "Status", m.name "Name", 
	m.markerType "Type"
    from #markersAll m, ACC_Accession a
    where m._Marker_key = a._Object_key
    and a._MGIType_key = 2
    and a.prefixPart = "MGI:"
    and a._LogicalDB_key = 1
    and a.preferred = 1
    order by m.symbol
    '''
results = db.sql(cmd, 'auto')

for r in results:
    line = []
    for h in headers:
	val = r[h]
	if val is None:
	    line.append("NULL")
	else:
	    line.append(str(r[h]))
    fp.write(TAB.join(line)+CRT)

reportlib.finish_nonps(fp)	# non-postscript file
db.useOneConnection(0)

