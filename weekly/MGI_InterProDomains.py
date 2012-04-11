#!/usr/local/bin/python

'''
#
# MGI_InterProDomains.py
#
# Report:
#       Tab-delimited file of InterPro Domains
#
# Usage:
#       MGI_InterProDomains.py
#
# History:
#
# lec	04/11/2012
#	- TR 11035/move from sql format to tab-delimited format
#
'''
 
import sys
import os
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

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

results = db.sql('''
	select a.accID, t.term
	from VOC_Vocab v, VOC_Term t, ACC_Accession a 
	where v.name = 'InterPro Domains' 
	and v._Vocab_key = t._Vocab_key 
	and t._Term_key = a._Object_key 
	and a._MGIType_key = 13
	order by t.sequenceNum
	''', 'auto')

for r in results:
	fp.write(r['accID'] + TAB)
	fp.write(r['term'] + CRT)

reportlib.finish_nonps(fp)

