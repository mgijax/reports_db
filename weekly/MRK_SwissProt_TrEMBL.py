
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
import db

db.setTrace()

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# Retrieve MGI Accession number, Marker symbol, name, etc.

results = db.sql('''
       select m.mgiID, m.symbol, m.status, m.name, m.chromosome, m.cmoffset, a.accID 
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

                if r['cmoffset'] == -1.0:
                        cmoffset = 'syntenic'
                elif r['cmoffset'] == -999.0:
                        cmoffset = 'N/A'
                else:
                        cmoffset = str(r['cmoffset'])

                fp.write(mgi_utils.prvalue(r['mgiID']) + reportlib.TAB + \
                         mgi_utils.prvalue(r['symbol']) + reportlib.TAB + \
                         mgi_utils.prvalue(str.upper(r['status'][0])) + reportlib.TAB + \
                         mgi_utils.prvalue(r['name']) + reportlib.TAB + \
                         mgi_utils.prvalue(cmoffset) + reportlib.TAB + \
                         mgi_utils.prvalue(r['chromosome']) + reportlib.TAB)

                prevMarker = r['mgiID']

        if len(sequence) > 0:
                sequence = sequence + ' '

        sequence = sequence + r['accID']

fp.write(sequence + reportlib.CRT)	# Don't forget to write out the last one
reportlib.finish_nonps(fp)
