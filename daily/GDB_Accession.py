#!/usr/local/bin/python

'''
#
# GDB_Accession.py 11/16/98
#
# Report:
#       Tab-delimited file of MGI and GDB Markers and Accession numbers
#	for existing Mouse/Human homologies.
#
# Usage:
#       GDB_Accession.py
#
# Used by:
#	The folks at GDB to provide HTML links from GDB Web detail pages
#	to MGI Web detail pages by MGI or GDB Accession numbers.
#	Report any changes in format/content to John Campbell (johnc@gdb.org).
#
# Notes:
#
# 1.  Read marker information into temp table
# 2.  Load MGI and GDB Acc# into dictionaries
# 3.  Sort temp file and process
#
# History:
#
# lec	07/02/2002
#	- TR 3855; add human LocusLink ID column
#
# lec	01/13/98
#	- added comments
#
'''
 
import sys
import os
import string
import db
import reportlib

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []

cmds.append('select distinct gdbSymbol = m1.symbol, gdbKey = m1._Marker_key, ' + \
            'mgiSymbol = m2.symbol, mgiName = m2.name, mgiKey = m2._Marker_key ' + \
            'into #homology ' + \
            'from HMD_Homology h1, HMD_Homology h2, ' + \
            'HMD_Homology_Marker hm1, HMD_Homology_Marker hm2, ' + \
            'MRK_Marker m1, MRK_Marker m2 ' + \
            'where m1._Organism_key = 2 ' + \
            'and m1._Marker_key = hm1._Marker_key ' + \
            'and hm1._Homology_key = h1._Homology_key ' + \
            'and h1._Class_key = h2._Class_key ' + \
            'and h2._Homology_key = hm2._Homology_key ' + \
            'and hm2._Marker_key = m2._Marker_key ' + \
            'and m2._Organism_key = 1')

cmds.append('select a.accID, a._Object_key from MRK_Acc_View a, #homology h ' + 
	    'where h.gdbKey = a._Object_key and a.prefixPart = "GDB:" and a._LogicalDB_key = 2 and a.preferred = 1')

cmds.append('select a.accID, a._Object_key from MRK_Acc_View a, #homology h ' + 
	    'where h.mgiKey = a._Object_key and a.prefixPart = "MGI:" and a._LogicalDB_key = 1 and a.preferred = 1')

cmds.append('select a.accID, a._Object_key from MRK_Acc_View a, #homology h ' + 
	    'where h.gdbKey = a._Object_key and a._LogicalDB_key = 24')

cmds.append('select * from #homology order by mgiSymbol')

results = db.sql(cmds, 'auto')

gdb = {}
for r in results[1]:
	gdb[r['_Object_key']] = r['accID']

mgi = {}
for r in results[2]:
	mgi[r['_Object_key']] = r['accID']

ll = {}
for r in results[3]:
	ll[r['_Object_key']] = r['accID']


for r in results[4]:

	if mgi.has_key(r['mgiKey']):

		fp.write(mgi[r['mgiKey']] + reportlib.TAB + \
	         	r['mgiSymbol'] + reportlib.TAB + \
	         	r['mgiName'][0:50] + reportlib.TAB)

		if gdb.has_key(r['gdbKey']):
			fp.write(gdb[r['gdbKey']])

		fp.write(reportlib.TAB + \
	         	r['gdbSymbol'] + reportlib.TAB)

		if ll.has_key(r['gdbKey']):
			fp.write(ll[r['gdbKey']])

		fp.write (reportlib.CRT)

reportlib.finish_nonps(fp)

