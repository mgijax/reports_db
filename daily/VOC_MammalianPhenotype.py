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
# Notes:
#	- all reports use mgireport directory for output file
#	- all reports use db default of public login
#	- all reports use server/database default of environment
#	- use lowercase for all SQL commands (i.e. select not SELECT)
#	- all public SQL reports require the header and footer
#	- all private SQL reports require the header
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
import db
import reportlib

CRT = reportlib.CRT
TAB = reportlib.TAB

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

results = db.sql('select t._Term_key, t.note, t.sequenceNum ' + \
	'from VOC_Text t, VOC_Term m ' + \
	'where m._Vocab_key = 5 ' + \
	'and m._Term_key = t._Term_key ' + \
	'order by t._Term_key, t.sequenceNum', 'auto')
mpnotes = {}
for r in results:
    key = r['_Term_key']
    value = r['note']

    if not mpnotes.has_key(key):
	mpnotes[key] = []

    mpnotes[key].append(value)

results = db.sql('select a.accID, t._Term_key, t.term from VOC_Term t, ACC_Accession a ' + \
	'where t._Vocab_key = 5 ' + \
	'and t._Term_key = a._Object_key ' + \
	'and a._MGIType_key = 13 ' + \
	'and a.preferred = 1 ' + \
	'order by a.accID', 'auto')

for r in results:

    key = r['_Term_key']

    fp.write(r['accID'] + TAB + r['term'] + TAB)

    if mpnotes.has_key(key):
	fp.write(string.join(mpnotes[key],''))

    fp.write(CRT)

reportlib.finish_nonps(fp)	# non-postscript file

