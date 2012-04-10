#!/usr/local/bin/python

'''
#
# MRK_Synonym.py 07/01/2008
#
# Report:
#       TR 9120
#	Tab-delimited version of existing report, MRK_Synonym.sql
#	"Marker Symbols and their Synonyms (ordered by Marker Symbol)"
#
# Usage:
#       MRK_Synonym.py
#
# History:
#
# jer	07/01/2008
#	- created
#
'''
 
import sys 
import os
import string
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


TAB = reportlib.TAB
CRT = reportlib.CRT

#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

headers = [
    "MGI Accession ID", 
    "Chr",
    "cM Position",
    "Marker Symbol",
    "Synonym"
    ]

fp.write(TAB.join(headers) + CRT)

cmd = '''
    select a.accID, m.chromosome,
        case
        when o.offset >= 0 then str(o.offset,10,2)
        when o.offset = -999.0 then "       N/A"
        when o.offset = -1.0 then "  syntenic"
        end as cmPosition,
    m.symbol,
    substring(s.synonym,1,90) as synonym
    from MRK_Marker m, ACC_Accession a, MGI_Synonym s, MGI_SynonymType st, MRK_Offset o 
    where m._Organism_key = 1 
    and m._Marker_key = a._Object_key 
    and a._MGIType_key = 2
    and a.prefixPart = "MGI:" 
    and a._LogicalDB_key = 1
    and a.preferred = 1 
    and m._Marker_key = o._Marker_key 
    and o.source = 0 
    and m._Marker_key = s._Object_key
    and s._MGIType_key = 2
    and s._SynonymType_key = st._SynonymType_key
    and st.synonymType = "exact"
    order by m.symbol
    '''
results = db.sql(cmd, 'auto')

for r in results:
    fp.write(r['accID'] + TAB)
    fp.write(r['chromosome'] + TAB)
    fp.write(r['cmPosition'] + TAB)
    fp.write(r['symbol'] + TAB)
    fp.write(r['synonym'] + CRT)

reportlib.finish_nonps(fp)	# non-postscript file
db.useOneConnection(0)

