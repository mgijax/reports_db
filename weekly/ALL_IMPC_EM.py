
'''
# ALL_IMPC_EM.py
#
# Report:
#       Tab-delimited file of all IMPC EM Alleles
#
# Usage:
#       ALL_IMPC.py
#
# Used by:
#       IMPC - Alba Gomez
#
# Output Format:
#
#   1. Allele ID
#   2. Allele symbol
#   3. Allelel name
#   4. Marker ID
#   5. Marker symbol
#   6. colony ID - if there is one
#
# 10/8/2018	sc
#       - TR12115 Load IMPC EM CRISPR alleles
#
'''

import sys
import os
import reportlib
import db

db.setTrace()

CRT = reportlib.CRT
TAB = reportlib.TAB

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('Allele ID' + TAB)
fp.write('Allele Symbol' + TAB)
fp.write('Allele Name' + TAB)
fp.write('Marker ID' + TAB)
fp.write('Marker Symbol' + TAB)
fp.write('Colony ID' + CRT)

db.sql('''select a._Allele_key, a.symbol as aSymbol, a.name, 
        a1.accid as alleleID, m.symbol as mSymbol, a2.accid as markerID
        into temporary table impc
        from ALL_Allele a, MRK_MArker m, ACC_Accession a1, ACC_Accession a2
        where a._collection_key  = 24755824 --IMPC
        and a._Allele_Type_key = 11927650 --Endonuclease-mediated
        and a._Allele_Status_key in (847114, 3983021)
        and a._Marker_key = m._Marker_key
        and a._Allele_key = a1._Object_key
        and a1._MGIType_key = 11
        and a1.preferred = 1
        and a1._LogicalDB_key = 1
        and a1.prefixPart = 'MGI:'
        and a._Marker_key = a2._Object_key
        and a2._MGIType_key = 2
        and a2.preferred = 1
        and a2._LogicalDB_key = 1
        and a2.prefixPart = 'MGI:' ''', None)

db.sql('''create index idx1 on impc(_Allele_key) ''', None)

db.sql('''select _Object_key as _Allele_key, note
        into temporary table cidNote
        from MGI_Note n
        where n._Notetype_key = 1041
        ''', None)

db.sql('''create index idx2 on cidNote(_Allele_key)''', None)

results = db.sql('''select i.*, n.note as cidNote
        from impc i
        left outer join cidNote n on (i._Allele_key = n._Allele_key) 
        order by i.aSymbol''', 'auto')

for r in results:	
    alleleID = r['alleleID']
    aSymbol = r['aSymbol']
    aName = r['name']
    markerID = r['markerID'] 
    mSymbol = r['mSymbol']
    cidNote = r['cidNote']
    if cidNote == None:
        cidNote = ''

    fp.write('%s%s%s%s%s%s%s%s%s%s%s%s' % (alleleID, TAB, aSymbol, TAB, aName, TAB, markerID, TAB, mSymbol, TAB, cidNote, CRT))

reportlib.finish_nonps(fp)
