#!/usr/local/bin/python

'''
#
# VOC_PhenoSlim.py 04/24/2002
#
# Report:
#       TR 3623 - PhenoSlim Terms
#
# Usage:
#       VOC_PhenoSlim.py
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
# lec	04/24/2002
#	- created
#
'''
 
import sys 
import os
import string
import db
import reportlib
import SimpleVocab

CRT = reportlib.CRT
TAB = reportlib.TAB

fpSQL = reportlib.init(sys.argv[0], fileExt = '.sql.rpt', outputdir = os.environ['REPORTOUTPUTDIR'])
fpTAB = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

phenoslim = SimpleVocab.PhenoSlimVocab()
for t in SimpleVocab.getPhenoslimTabDelim(phenoslim, headings = []):
	fpTAB.write(t + CRT)

for t in SimpleVocab.getPhenoslimText(phenoslim):
	fpSQL.write(t + CRT)

fpSQL.write(CRT + '(' + `len(phenoslim)` + ' rows affected)' + CRT)

reportlib.trailer(fpSQL)
reportlib.finish_nonps(fpSQL)	# non-postscript file

