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
# lec	06/14/2013
#	- interpro domain may contain > 1 marker
#
# lec	04/30/2012	
#	- TR11035/merge with MRK_InterPro.py
#
# lec	04/11/2012
#	- TR 11035/move from sql format to tab-delimited format
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
	db.setAutoTranslate(False)
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

#
# mouse markers associated with interpro domain
#
results = db.sql('''
        select distinct a1.accID, m.symbol, a._Term_key
        from MRK_Marker m, ACC_Accession a1, VOC_Annot a
        where m._Organism_key = 1 
        and m._Marker_key = a1._Object_key 
        and a1._MGIType_key = 2 
        and a1._LogicalDB_key = 1 
        and a1.prefixPart = 'MGI:' 
        and a1.preferred = 1 
        and m._Marker_key = a._Object_key 
        and a._AnnotType_key = 1003 
        ''', 'auto')
markers = {}
for r in results:
    key = r['_Term_key']
    value = r
    if not markers.has_key(key):
        markers[key] = []
    markers[key].append(r)

#
# interpro domain
#
results = db.sql('''
	select a.accID, t.term, t._Term_key
	from VOC_Vocab v, VOC_Term t, ACC_Accession a 
	where v.name = 'InterPro Domains' 
	and v._Vocab_key = t._Vocab_key 
	and t._Term_key = a._Object_key 
	and a._MGIType_key = 13
	order by t.sequenceNum
	''', 'auto')

for r in results:

	key = r['_Term_key']

	if markers.has_key(key):
            for k in markers[key]:
	        fp.write(r['accID'] + TAB)
	        fp.write(r['term'] + TAB)
	        fp.write(k['accID'] + TAB)
	        fp.write(k['symbol'] + CRT)
	else:
	    fp.write(r['accID'] + TAB)
	    fp.write(r['term'] + TAB)
	    fp.write(TAB + CRT)

reportlib.finish_nonps(fp)

