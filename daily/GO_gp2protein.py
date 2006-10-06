#!/usr/local/bin/python

'''
#
# GO_gp2protein.py (TR 4877, originally TR 3659)
#
# Report:
#       Tab-delimited file
#
# Usage:
#       GO_gp2protein.py
#
# Output
#
#	A tab-delimited file in this format:
#	field 1: MGI Marker ID
#	field 2: SWP:#####;SWP:#####;...
# Used by:
#
# Notes:
#
# History:
#	Created for TR3659.  
#	Report the MGI ID, SwissProt ID(s), and GO ID(s) for all mouse
#	markers that have associated SwissProt sequences.

'''
 
import sys
import string
import os
import db
import reportlib

#
# Main
#

TAB = reportlib.TAB
CRT = reportlib.CRT

fp = reportlib.init('gp2protein', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# Retrieve Markers with GO Annotations that have SP IDs

cmds = []

cmds.append('select distinct m._Marker_key, a.accID, spID = a2.accID ' + \
	'into #markers ' + \
	'from MRK_Marker m, ACC_Accession a, ACC_Accession a2 ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1 ' + \
	'and m._Marker_key = a2._Object_key ' + \
	'and a2._MGIType_key = 2 ' + \
	'and a2._LogicalDB_key = 13 ' + \
	'and exists (select 1 from VOC_Annot v ' + \
	'where v._AnnotType_key = 1000 ' + \
	'and m._Marker_key = v._Object_key) ' + \
	'order by a.accID')

cmds.append('select _Marker_key, spID from #markers')
cmds.append('select distinct _Marker_key, accID from #markers')

results = db.sql(cmds, 'auto')

spIDs = {}
for r in results[-2]:
    key = r['_Marker_key']
    value = 'SWP:' + r['spID']
    if not spIDs.has_key(key):
	spIDs[key] = []
    spIDs[key].append(value)

for r in results[-1]:

    fp.write(r['accID'] + reportlib.TAB + \
	string.join(spIDs[r['_Marker_key']], ';') + CRT)

reportlib.finish_nonps(fp)

