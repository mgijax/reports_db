#!/usr/local/bin/python

'''
#
# MRK_Sequence.py 03/02/99
#
# Report:
#       Tab-delimited file
#       Mouse Markers and their Nucleotide Sequence Accession numbers.
#	and their UniGene Accession numbers (TR 1631).
#
# Usage:
#       MRK_Sequence.py
#
# Used by:
#       Those establishing relationships between MGI Markers
#	and Nucleotide Sequences.
#
# Notes:
#
# History:
#
# lec	05/30/1999
#	- TR 1631; add UniGene Accession numbers
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
import db
import mgi_utils
import reportlib

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

# Retrieve MGI Accession number, Marker symbol, name, offset and 
# all associated Nucleotide Seq IDs for Markers

cmd = 'select distinct m.mgiID, m.symbol, m.name, m.chromosome, m.offset, a.accID, type = "N" ' + \
      'from MRK_Mouse_View m, ACC_Accession a ' + \
      'where m._Marker_key = a._Object_key ' + \
      'and a._MGIType_key = 2 ' + \
      'and a._LogicalDB_key = 9 ' + \
      'union ' + \
      'select distinct m.mgiID, m.symbol, m.name, m.chromosome, m.offset, a.accID, type = "U" ' + \
      'from MRK_Mouse_View m, ACC_Accession a ' + \
      'where m._Marker_key = a._Object_key ' + \
      'and a._MGIType_key = 2 ' + \
      'and a._LogicalDB_key = 23 ' + \
      'order by m.symbol, m.mgiID, type, a.accID'

results = db.sql(cmd, 'auto')

prevMarker = ''
sequence = ''
unigene = ''

for r in results:

	if prevMarker != r['mgiID']:
		if len(sequence) > 0:
			fp.write(mgi_utils.prvalue(sequence) + reportlib.TAB)
		sequence = ''

		if len(unigene) > 0:
			fp.write(mgi_utils.prvalue(unigene))
		unigene = ''

		if prevMarker != '':
			fp.write(reportlib.CRT)

		if r['offset'] == -1.0:
			offset = 'syntenic'
		elif r['offset'] == -999.0:
			offset = 'N/A'
		else:
			offset = str(r['offset'])

		fp.write(mgi_utils.prvalue(r['mgiID']) + reportlib.TAB + \
	        	 mgi_utils.prvalue(r['symbol']) + reportlib.TAB + \
	                 mgi_utils.prvalue(r['name']) + reportlib.TAB + \
	                 mgi_utils.prvalue(offset) + reportlib.TAB + \
	                 mgi_utils.prvalue(r['chromosome']) + reportlib.TAB)

		prevMarker = r['mgiID']

	if r['type'] == 'N':
		if len(sequence) > 0:
			sequence = sequence + ' '

        	sequence = sequence + r['accID']

	if r['type'] == 'U':
		if len(unigene) > 0:
			unigene = unigene + ' '

        	unigene = unigene + r['accID']

fp.write(sequence + reportlib.TAB)	# Don't forget to write out the last one
fp.write(unigene + reportlib.CRT)	# Don't forget to write out the last one
reportlib.finish_nonps(fp)

