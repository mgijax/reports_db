
'''
#
# MP_EMAPA.py
#
# Report:
#       TR12867 Report of MP-EMAOA mappings for Peter Robinson
#
# Usage:
#       MP_EMAPA.py
#
# History:
#
# lec	08/29/2018
#	- TR12943/move report into reports_db product
#
# sc	05/09/2018
#	- created
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
results = db.sql('''
        select a1.accid as mpID, t1.term as mpTerm, a2.accid as emapaID, t2.term as emapaTerm
        from MGI_Relationship r, ACC_Accession a1, VOC_Term t1, ACC_Accession a2, VOC_Term t2
        where r._Category_key = 1007
        and r._Object_key_1 = a1._Object_key
        and a1._MGIType_key = 13
        and a1._LogicalDB_key = 34
        and a1.preferred = 1
        and r._Object_key_1 = t1._Term_key
        and r._Object_key_2 = a2._Object_key
        and a2._MGIType_key = 13
        and a2._LogicalDB_key = 169
        and a2.preferred = 1
        and r._Object_key_2 = t2._Term_key
        order by mpID, emapaID
        ''', 'auto')
    
for r in results:
    fp.write(r['mpID'] + TAB + r['mpTerm'] + TAB + r['emapaID'] + TAB + r['emapaTerm'] + CRT)

reportlib.finish_nonps(fp)	# non-postscript file
