#!/usr/local/bin/python

'''
#
# VOC_MammalianPhenotype.py
#
# Report:
#       TR 6670
#
# Usage:
#       VOC_MammalianPhenotype.py
#
# History:
#
# lec	03/18/2005
#	- created
#
'''
 
import sys 
import os
import string
import reportlib
import db

db.setTrace()

CRT = reportlib.CRT
TAB = reportlib.TAB

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

results = db.sql('''
	select m._Term_key, m.note
	from VOC_Term m 
	where m._Vocab_key = 5 
	order by m._Term_key
	''', 'auto')
mpnotes = {}
for r in results:
    key = r['_Term_key']
    value = r['note']
    mpnotes[key] = value

results = db.sql('''
	select a.accID, t._Term_key, t.term from VOC_Term t, ACC_Accession a 
	where t._Vocab_key = 5 
	and t._Term_key = a._Object_key 
	and a._MGIType_key = 13 
	and a.preferred = 1 
	order by a.accID
	''', 'auto')

for r in results:

    key = r['_Term_key']

    fp.write(r['accID'] + TAB + r['term'] + TAB)

    if mpnotes.has_key(key):
	fp.write(mpnotes[key])

    fp.write(CRT)

reportlib.finish_nonps(fp)	# non-postscript file

