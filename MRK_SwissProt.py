#!/usr/local/bin/python

'''
#
# MRK_SwissProt.py 03/02/99
#
# Report:
#       Tab-delimited file
#       Mouse Markers and their SWISS-PROT Accession numbers.
#
# Usage:
#       MRK_SwissProt.py
#
# Generated from:
#       Editing Interface Nightly Reports
#
# Used by:
#       Those establishing relationships between MGI Markers
#	and Protein Sequences in SWISS-PROT.
#
# Notes:
#
# History:
#
# lec	03/02/1999
#	- modified to use one query only; consistent w/ MRK_Sequence.py
#
# lec	07/08/98
#	- created; requested by User Support
#
'''
 
import sys
import os
import mgdlib
import accessionlib
import reportlib

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['FTPREPORTDIR'])

MGIType = accessionlib.get_MGIType_key('Marker')
Sequence = accessionlib.get_LogicalDB_key('SWISS-PROT')

# Retrieve MGI Accession number, Marker symbol, name, etc.

cmd = 'select m.mgiID, m.symbol, m.name, m.chromosome, m.offset, a.accID ' + \
      'from MRK_Mouse_View m, ACC_Accession a ' + \
      'where m._Marker_key = a._Object_key ' + \
      ' and a._MGIType_key = ' + `MGIType` + \
      ' and a._LogicalDB_key = ' + `Sequence` + \
      'order by m.symbol, m.mgiID'

results = mgdlib.sql(cmd, 'auto')

prevMarker = ''
sequence = ''

for r in results:

	if prevMarker != r['mgiID']:
		if len(sequence) > 0:
			fp.write(mgdlib.prvalue(sequence) + reportlib.CRT)
		sequence = ''

		if r['offset'] == -1.0:
			offset = 'syntenic'
		elif r['offset'] == -999.0:
			offset = 'N/A'
		else:
			offset = str(r['offset'])

		fp.write(mgdlib.prvalue(r['mgiID']) + reportlib.TAB + \
	        	 mgdlib.prvalue(r['symbol']) + reportlib.TAB + \
	                 mgdlib.prvalue(r['name']) + reportlib.TAB + \
	                 mgdlib.prvalue(offset) + reportlib.TAB + \
	                 mgdlib.prvalue(r['chromosome']) + reportlib.TAB)

		prevMarker = r['mgiID']

	if len(sequence) > 0:
		sequence = sequence + ' '

        sequence = sequence + r['accID']

fp.write(sequence + reportlib.CRT)	# Don't forget to write out the last one
reportlib.finish_nonps(fp)

