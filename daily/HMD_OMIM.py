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
	
db.useOneConnection(1)
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

cmds = []
cmds.append('select distinct mouseKey = h1._Marker_key, mouseSymbol = m1.symbol, ' +
	'humanKey = h2._Marker_key, humanSymbol = m2.symbol ' + \
	'into #homology ' +
        'from HMD_Homology r1, HMD_Homology_Marker h1, ' + \
        'HMD_Homology r2, HMD_Homology_Marker h2, ' + \
        'MRK_Marker m1, MRK_Marker m2 ' + \
        'where m1._Organism_key = 1 ' + \
        'and m1._Marker_key = h1._Marker_key ' + \
        'and h1._Homology_key = r1._Homology_key ' + \
        'and r1._Class_key = r2._Class_key ' + \
        'and r2._Homology_key = h2._Homology_key ' + \
        'and h2._Marker_key = m2._Marker_key ' + \
        'and m2._Organism_key = 2 ')

cmds.append('create nonclustered index idx1 on #homology(mouseKey)')
cmds.append('create nonclustered index idx2 on #homology(humanKey)')
cmds.append('create nonclustered index idx3 on #homology(mouseSymbol)')

db.sql(cmds, None)

##

results = db.sql('select h.mouseKey, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.preferred = 1 ', 'auto')
mgiIDs = {}
for r in results:
    key = r['mouseKey']
    value = r['accID']
    mgiIDs[key] = value

results = db.sql('select h.humanKey, l.mim ' + \
	'from #homology h, ACC_Accession a, radar_2..DP_LL l ' + \
	'where h._Marker_key = h._Object_key ' + \
	'and h._MGIType_key = 2 ' + \
	'and h._LogicalDB_key = 24 ' + \
	'and h.accID = l.locusID ' + \
	'and l.mim is not null ', 'auto')
mimIDs = {}
for r in results:
    key = r['humanKey']
    value = r['accID']
    mimIDs[key] = value

results = db.sql('select h.* from #homology order by mouseSymbol', 'auto')
count = 0
for r in results:
    if not mimIDs.has_key(h.humanKey):
	continue

    fp.write(string.ljust(mgiIDs[r['mouseKey']], 30))
    fp.write(SPACE)
    fp.write(string.ljust(r['mouseSymbol'], 25))
    fp.write(SPACE)
    fp.write(string.ljust(r['humanSymbol'], 25))
    fp.write(SPACE)
    fp.write(string.ljust(mimIDs[r['humanKey']]), 10))
    fp.write(SPACE)
    fp.write(CRT)
    count = count + 1

fp.write(CRT + '(%d rows affected)' % (count) + CRT)
reportlib.trailer(fp)
reportlib.finish_nonps(fp)
db.useOneConnection(0)

