#!/usr/local/bin/python

'''
#
# MRK_SwissProt_Sequence.py 03/25/99
#
# Report:
#       Tab-delimited file
#       All Mouse Genes and their Nucleotide Sequence Accession numbers
#	(when available).  Includes Withdrawn Markers.
#
# Usage:
#       MRK_SwissProt_Sequence.py
#
# Used by:
#       Swiss-Prot - Amos Bairoch - TR 444
#
# Notes:
#
# History:
#
# lec	02/15/2000
#	- TR 1278; include new Complex/Cluster/Region Marker type
#
# lec	03/25/1999
#	- created - see TR 444
#
'''
 
import sys
import os
import db
import mgi_utils
import reportlib

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []

# Retrieve all Gene Markers which have MGI Acc IDs
# Union
# Retrieve all withdrawn Gene Markers (which don't have MGI Acc IDs, except for splits)

cmds.append('select m._Marker_key, m.mgiID, m.symbol, m.name, m.chromosome ' +
      'into #markers ' +
      'from MRK_Mouse_View m ' +
      'where m._Marker_Type_key = 1 ' +
      'union ' +
      'select m._Marker_key, mgiID = null, m.symbol, m.name, m.chromosome ' +
      'from MRK_Marker m ' +
      'where m._Marker_Type_key = 1 ' +
      'and m._Organism_key = 1 ' +
      'and m._Marker_Status_key = 2 ' +
      'and not exists (select a.* from MRK_Acc_View a where a._Object_key = m._Marker_key) ' +
      'order by m.symbol, m.mgiID')

# Retrieve any Nucleotide Seq IDs which exist for Markers

cmds.append('select m.*, a.accID ' +
	'from #markers m, ACC_Accession a ' +
        'where m._Marker_key = a._Object_key ' +
        'and a._MGIType_key = 2 ' + \
        'and a._LogicalDB_key = 9 ' + \
	'union ' +
	'select m.*, accID = null ' +
	'from #markers m ' +
	'where not exists (select a.* from ACC_Accession a ' +
        'where m._Marker_key = a._Object_key ' +
        'and a._MGIType_key = 2 ' + \
        'and a._LogicalDB_key = 9) ' + \
        'order by m.symbol, m.mgiID, a.accID')

results = db.sql(cmds, 'auto')

prevMarker = ''
sequence = ''

for r in results[1]:

	if r['mgiID'] is None or r['accID'] is None:
		if len(sequence) > 0:
			fp.write(mgi_utils.prvalue(sequence) + reportlib.CRT)
		sequence = ''

		fp.write(mgi_utils.prvalue(r['symbol']) + reportlib.TAB + \
	                 mgi_utils.prvalue(r['name']) + reportlib.TAB + \
	        	 mgi_utils.prvalue(r['mgiID']) + reportlib.TAB + \
	                 mgi_utils.prvalue(r['chromosome']) + reportlib.TAB + reportlib.CRT)

		prevMarker = r['mgiID']

	elif prevMarker != r['mgiID']:
		if len(sequence) > 0:
			fp.write(mgi_utils.prvalue(sequence) + reportlib.CRT)
		sequence = ''

		fp.write(mgi_utils.prvalue(r['symbol']) + reportlib.TAB + \
	                 mgi_utils.prvalue(r['name']) + reportlib.TAB + \
	        	 mgi_utils.prvalue(r['mgiID']) + reportlib.TAB + \
	                 mgi_utils.prvalue(r['chromosome']) + reportlib.TAB)

		prevMarker = r['mgiID']

	if len(sequence) > 0:
		sequence = sequence + ' '

	if r['accID'] is not None:
        	sequence = sequence + r['accID']


fp.write(sequence + reportlib.CRT)	# Don't forget to write out the last one
reportlib.finish_nonps(fp)

