#!/usr/local/bin/python

'''
#
# MRK_Sequence.py 03/02/99
#
# Report:
#       Tab-delimited file
#       Mouse Markers and their Nucleotide Sequence Accession numbers.
#
# Usage:
#       MRK_Sequence.py
#
# Generated from:
#       Editing Interface Nightly Reports
#
# Used by:
#       Those establishing relationships between MGI Markers
#	and Nucleotide Sequences.
#
# Notes:
#
# History:
#
# lec	03/02/1999
#	- TR 130; now use direct Marker-Seq ID relationships
#
# lec	01/18/1999
#	- missing Segment records where relationship is null
#
# lec	01/13/98
#	- added comments section
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

fp = reportlib.init(sys.argv[0], None, os.environ['FTPREPORTDIR'])

Marker = accessionlib.get_MGIType_key('Marker')
Sequence = accessionlib.get_LogicalDB_key('Sequence DB')

# Retrieve MGI Accession number, Marker symbol, name, offset and 
# all associated Nucleotide Seq IDs for Markers

cmd = 'select distinct m.mgiID, m.symbol, m.name, m.chromosome, m.offset, a.accID ' + \
      'from MRK_Mouse_View m, ACC_Accession a ' + \
      'where m._Marker_key = a._Object_key ' + \
      ' and a._MGIType_key = ' + `Marker` + \
      ' and a._LogicalDB_key = ' + `Sequence` + \
      'order by m.symbol, m.mgiID, a.accID'

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

