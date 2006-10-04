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
# lec	10/25/2005
#	- MGI 3.5; now loading Human RefSeqs directly into MGI
#
# lec	01/04/2004
#	- TR 5939; LocusLink->EntrezGene
#
# lec   07/20/2004
#       - created
#
'''

import sys
import os
import string
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

db.sql('select distinct mouseKey = h1._Marker_key, mouseSymbol = m1.symbol, ' +
	'humanKey = h2._Marker_key, humanSymbol = m2.symbol ' + \
	'into #homology ' +
        'from MRK_Homology_Cache h1, MRK_Homology_Cache h2, ' + \
        'MRK_Marker m1, MRK_Marker m2 ' + \
        'where h1._Organism_key = 1 ' + \
        'and h1._Class_key = h2._Class_key ' + \
        'and h1._Marker_key = m1._Marker_key ' + \
        'and h2._Marker_key = m2._Marker_key ', None)

db.sql('create nonclustered index idx1 on #homology(mouseKey)', None)
db.sql('create nonclustered index idx2 on #homology(humanKey)', None)
db.sql('create nonclustered index idx3 on #homology(mouseSymbol)', None)

##

results = db.sql('select h.mouseKey, a.accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.mouseKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.preferred = 1 ', 'auto')
mgiIDs = {}
for r in results:
    key = r['mouseKey']
    value = r['accID']
    mgiIDs[key] = value

results = db.sql('select h.humanKey, accID ' + \
	'from #homology h, ACC_Accession a ' + \
	'where h.humanKey = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 15 ', 'auto')
mimIDs = {}
for r in results:
    key = r['humanKey']
    value = r['accID']
    if not mimIDs.has_key(key):
	mimIDs[key] = []
    mimIDs[key].append(value)

results = db.sql('select h.* from #homology h order by h.mouseSymbol', 'auto')
count = 0
for r in results:
    if not mimIDs.has_key(r['humanKey']):
	continue
    fp.write(string.ljust(mgiIDs[r['mouseKey']], 30))
    fp.write(SPACE)
    fp.write(string.ljust(r['mouseSymbol'], 25))
    fp.write(SPACE)
    fp.write(string.ljust(r['humanSymbol'], 25))
    fp.write(SPACE)
    fp.write(string.ljust(string.join(mimIDs[r['humanKey']], ','), 10))
    fp.write(SPACE)
    fp.write(CRT)
    count = count + 1

fp.write(CRT + '(%d rows affected)' % (count) + CRT)
reportlib.trailer(fp)
reportlib.finish_nonps(fp)
db.useOneConnection(0)

