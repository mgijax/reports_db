#!/usr/local/bin/python

'''
#
# HMD_OMIM.py
#
# Report:
#	TR 6039
#
# Usage:
#       HMD_OMIM.py
#
# Notes:
#       - all reports use db default of public login
#       - all reports use server/database default of environment
#       - use lowercase for all SQL commands (i.e. select not SELECT)
#       - all public SQL reports require the header and footer
#       - all private SQL reports require the header
#
# History:
#
# lec   07/20/2004
#       - created
#
'''

import sys
import os
import string
import regex
import db
import mgi_utils
import reportlib

CRT = reportlib.CRT
SPACE = reportlib.SPACE

reportTitle = 'Human and Mouse Orthology with Human OMIM IDs'
	
fp = reportlib.init(sys.argv[0], reportTitle, os.environ['REPORTOUTPUTDIR'])
	
fp.write(string.ljust('Mouse MGI Acc ID', 30))
fp.write(SPACE)
fp.write(string.ljust('Mouse Symbol', 25))
fp.write(SPACE)
fp.write(string.ljust('Human Symbol', 25))
fp.write(SPACE)
fp.write(string.ljust('OMIM ID', 10))
fp.write(SPACE)
fp.write(CRT)

fp.write(string.ljust('------------------', 30))
fp.write(SPACE)
fp.write(string.ljust('------------', 25))
fp.write(SPACE)
fp.write(string.ljust('------------', 25))
fp.write(SPACE)
fp.write(string.ljust('-------', 10))
fp.write(SPACE)
fp.write(CRT)

cmd = 'select distinct humanSymbol = m1.symbol, mouseMGI = a.accID, mouseSymbol = m2.symbol, l.mim ' + \
	'from HMD_Homology r1, HMD_Homology_Marker h1, HMD_Homology r2, HMD_Homology_Marker h2, ' + \
	'MRK_Marker m1, MRK_Marker m2, ACC_Accession a, ACC_Accession ha, radar_2..DP_LL l ' + \
	'where m1._Organism_key = 2 ' + \
	'and m1._Marker_key = h1._Marker_key ' + \
	'and h1._Homology_key = r1._Homology_key ' + \
	'and r1._Class_key = r2._Class_key ' + \
	'and r2._Homology_key = h2._Homology_key ' + \
	'and h2._Marker_key = m2._Marker_key ' + \
	'and m2._Organism_key = 1 ' + \
	'and m2._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.preferred = 1 ' + \
	'and m1._Marker_key = ha._Object_key ' + \
	'and ha._MGIType_key = 2 ' + \
	'and ha._LogicalDB_key = 24 ' + \
	'and ha.accID = l.locusID ' + \
	'and l.mim is not null ' + \
	'order by mouseSymbol'

results = db.sql(cmd, 'auto')
count = 0
for r in results:
    fp.write(string.ljust(r['mouseMGI'], 30))
    fp.write(SPACE)
    fp.write(string.ljust(r['mouseSymbol'], 25))
    fp.write(SPACE)
    fp.write(string.ljust(r['humanSymbol'], 25))
    fp.write(SPACE)
    fp.write(string.ljust(str(r['mim']), 10))
    fp.write(SPACE)
    fp.write(CRT)
    count = count + 1

fp.write(CRT + '(%d rows affected)' % (count) + CRT)
reportlib.trailer(fp)
reportlib.finish_nonps(fp)

