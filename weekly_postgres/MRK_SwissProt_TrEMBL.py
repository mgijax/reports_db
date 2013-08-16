#!/usr/local/bin/python

'''
#
# MRK_SwissProt_TrEMBL.py 03/26/2003
#
# Report:
#       Tab-delimited file
#       Mouse Markers and their SWISS-PROT and TrEMBL Accession numbers.
#
# Usage:
#       MRK_SwissProt_TrEMBL.py
#
# Used by:
#       Those establishing relationships between MGI Markers
#	and Protein Sequences in SWISS-PROT/TrEMBL.
#
# History:
#
# lec	03/26/2003
#	- TR 4270; copied from MRK_SwissProt.py
#
'''
 
import sys
import os
import string
import mgi_utils
import reportlib
import pg_db
db = pg_db
db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE()

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# Retrieve MGI Accession number, Marker symbol, name, etc.

results = db.sql('''
       select m.mgiID, m.symbol, m.status, m.name, m.chromosome, m.offset, a.accID 
       from MRK_Mouse_View m, ACC_Accession a 
       where m._Marker_key = a._Object_key 
       and a._MGIType_key = 2 
       and a._LogicalDB_key in (13, 41) 
       order by m.symbol, m.mgiID
       ''', 'auto')

prevMarker = ''
sequence = ''

for r in results:

	if prevMarker != r['mgiID']:
		if len(sequence) > 0:
			fp.write(mgi_utils.prvalue(sequence) + reportlib.CRT)
		sequence = ''

		if r['offset'] == -1.0:
			offset = 'syntenic'
		elif r['offset'] == -999.0:
			offset = 'N/A'
		else:
			offset = str(r['offset'])

		fp.write(mgi_utils.prvalue(r['mgiID']) + reportlib.TAB + \
	        	 mgi_utils.prvalue(r['symbol']) + reportlib.TAB + \
			 mgi_utils.prvalue(string.upper(r['status'][0])) + reportlib.TAB + \
	                 mgi_utils.prvalue(r['name']) + reportlib.TAB + \
	                 mgi_utils.prvalue(offset) + reportlib.TAB + \
	                 mgi_utils.prvalue(r['chromosome']) + reportlib.TAB)

		prevMarker = r['mgiID']

	if len(sequence) > 0:
		sequence = sequence + ' '

        sequence = sequence + r['accID']

fp.write(sequence + reportlib.CRT)	# Don't forget to write out the last one
reportlib.finish_nonps(fp)
