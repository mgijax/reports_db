#!/usr/local/bin/python

'''
#
# PRB_ESTMapped.py 06/29/2000
# TR 1734
#
# Report:
#       Tab-delimited file
#       Mapped ESTs
#
# Usage:
#       PRB_ESTMapped.py
#
# Used by:
#       Incyte
#
# Notes:
#
# History:
#
# lec	06/29/2000
#	- created
#
'''
 
import sys
import os
import db
import reportlib

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []

# Select all Markers for ESTs which have a map position

cmds.append('select p._Probe_key, m._Marker_key, m.symbol, m.name, m.chromosome, o.offset ' + \
	'into #marker ' + \
	'from PRB_Probe p, PRB_Marker pm, MRK_Marker m, MRK_Offset o ' + \
	'where p.name like "IMAGE clone%"' + \
	'and p._Probe_key = pm._Probe_key ' + \
	'and pm._Marker_key = m._Marker_key ' + \
	'and m.chromosome != "UN" ' + \
	'and m._Marker_key = o._Marker_key ' + \
	'and o.source = 0')

# Select MGI Acc ID for ESTs

cmds.append('select m.*, mgiID = a.accID ' + \
      'into #acc ' + \
      'from #marker m, ACC_Accession a ' + \
      'where m._Probe_key = a._Object_key ' + \
      'and a._MGIType_key = 3 ' + \
      'and a.prefixPart = "MGI:"' + \
      'and a._LogicalDB_key = 1 ' + \
      'and a.preferred = 1')

# Select GenBank ID for ESTs

cmds.append('select m.*, genbankID = a.accID ' + \
      'from #acc m, ACC_Accession a ' + \
      'where m._Probe_key = a._Object_key ' + \
      'and a._MGIType_key = 3 ' + \
      'and a._LogicalDB_key = 9 ' + \
      'order by symbol')

results = db.sql(cmds, 'auto')

for r in results[2]:

	if r['offset'] == -1.0:
		offset = 'syntenic'
	elif r['offset'] == -999.0:
		offset = 'N/A'
	else:
		offset = str(r['offset'])

	fp.write(r['mgiID'] + reportlib.TAB)
	fp.write(r['genbankID'] + reportlib.TAB)
	fp.write(r['chromosome'] + reportlib.TAB)
	fp.write(offset + reportlib.TAB)
	fp.write(r['symbol'] + reportlib.TAB)
	fp.write(r['name'] + reportlib.CRT)

reportlib.finish_nonps(fp)

