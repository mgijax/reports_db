#!/usr/local/bin/python

'''
#
# Report:
#       Tab-delimited file
#
#	1. MGI Strain
#	2. Literature Strain
#
# Usage:
#       MGI_StrainTranslations.py
#
# History:
#
# lec	08/02/2005
#	- Donna Maglott
#
'''
 
import sys 
import os
import string
import db
import reportlib

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], printHeading = 0, outputdir = os.environ['REPORTOUTPUTDIR'])

results = db.sql('select s.strain, t.badName ' + \
	'from MGI_Translation t, PRB_Strain s ' + \
	'where t._TranslationType_key = 1007 ' + \
	'and t._Object_key = s._Strain_key', 'auto')

for r in results:

	fp.write(r['strain'] + TAB + r['badName'] + CRT)

reportlib.finish_nonps(fp)	# non-postscript file
db.useOneConnection(0)
