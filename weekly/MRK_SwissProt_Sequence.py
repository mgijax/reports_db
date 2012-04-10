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
import string
import mgi_utils
import reportlib

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db


#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# Retrieve all Gene Markers which have MGI Acc IDs
# Union
# Retrieve all withdrawn Gene Markers (which don't have MGI Acc IDs, except for splits)

db.sql('''
      (
      select m._Marker_key, a.accID as mgiID, m.symbol, m.name, m.chromosome 
      into #markers 
      from MRK_Marker m, ACC_Accession a 
      where m._Marker_Type_key = 1 
      and m._Organism_key = 1 
      and m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 1 
      and a.prefixPart = "MGI:" 
      and a.preferred = 1 
      union 
      select m._Marker_key, null as mgiID, m.symbol, m.name, m.chromosome 
      from MRK_Marker m 
      where m._Marker_Type_key = 1 
      and m._Organism_key = 1 
      and m._Marker_Status_key = 2 
      and not exists (select 1 from ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2) 
      )
      order by symbol, mgiID
      ''', None)

# Retrieve any Nucleotide Seq IDs which exist for Markers

results = db.sql('''
	select m._Marker_key, a.accID 
	from #markers m, ACC_Accession a 
        where m._Marker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a._LogicalDB_key = 9 
	order by m._Marker_key, a.accID
	''', 'auto')
seqIDs = {}
for r in results:
    if not seqIDs.has_key(r['_Marker_key']):
	seqIDs[r['_Marker_key']] = []
    seqIDs[r['_Marker_key']].append(r['accID'])

results = db.sql('select * from #markers order by symbol, mgiID', 'auto')
for r in results:

	fp.write(mgi_utils.prvalue(r['symbol']) + reportlib.TAB + \
	         mgi_utils.prvalue(r['name']) + reportlib.TAB + \
	         mgi_utils.prvalue(r['mgiID']) + reportlib.TAB + \
	         mgi_utils.prvalue(r['chromosome']) + reportlib.TAB)

	if seqIDs.has_key(r['_Marker_key']):
		fp.write(string.join(seqIDs[r['_Marker_key']], ' '))

	fp.write(reportlib.CRT)

reportlib.finish_nonps(fp)

