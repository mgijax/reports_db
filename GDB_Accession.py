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
# lec	01/13/98
#	- added comments
#
'''
 
import sys
import os
import string
import mgdlib
import reportlib

def parseGDB(_tuple):
	global gdb

	gdb[_tuple['_Object_key']] = _tuple['accID']

def parseMGI(_tuple):
	global mgi

	mgi[_tuple['_Object_key']] = _tuple['accID']

def parseHomology(_tuple):

	try:
		gdbID = gdb[_tuple['gdbKey']]
	except:
		gdbID = 'None'

	# Might hit a Marker which doesn't have an Accession number

	try:
		fp.write(mgi[_tuple['mgiKey']] + reportlib.TAB + \
	         	_tuple['mgiSymbol'] + reportlib.TAB + \
	         	_tuple['mgiName'][0:50] + reportlib.TAB + \
	         	gdbID + reportlib.TAB + \
	         	_tuple['gdbSymbol'] + reportlib.CRT)
	except:
		pass

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

gdb = {}
mgi = {}

cmds = []
parsers = []

cmds.append('select distinct gdbSymbol = m1.symbol, gdbKey = m1._Marker_key, ' + \
            'mgiSymbol = m2.symbol, mgiName = m2.name, mgiKey = m2._Marker_key ' + \
            'into #homology ' + \
            'from HMD_Homology h1, HMD_Homology h2, ' + \
            'HMD_Homology_Marker hm1, HMD_Homology_Marker hm2, ' + \
            'MRK_Marker m1, MRK_Marker m2 ' + \
            'where m1._Species_key = 2 ' + \
            'and m1._Marker_key = hm1._Marker_key ' + \
            'and hm1._Homology_key = h1._Homology_key ' + \
            'and h1._Class_key = h2._Class_key ' + \
            'and h2._Homology_key = hm2._Homology_key ' + \
            'and hm2._Marker_key = m2._Marker_key ' + \
            'and m2._Species_key = 1')

cmds.append('select a.accID, a._Object_key from MRK_Acc_View a, #homology h ' + 
	    'where h.gdbKey = a._Object_key and a.prefixPart = "GDB:" and a.preferred = 1')

cmds.append('select a.accID, a._Object_key from MRK_Acc_View a, #homology h ' + 
	    'where h.mgiKey = a._Object_key and a.prefixPart = "MGI:" and a.preferred = 1')

cmds.append('select * from #homology order by mgiSymbol')

parsers.append(None)
parsers.append(parseGDB)
parsers.append(parseMGI)
parsers.append(parseHomology)
mgdlib.sql(cmds, parsers)
reportlib.finish_nonps(fp)

