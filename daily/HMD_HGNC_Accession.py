#!/usr/local/bin/python

'''
#
# HMD_HGNC_Accession.py
#
# Report:
#       Tab-delimited file of MGI and HGNC Markers and Accession numbers
#	for existing Mouse/Human orthologies.
#
# Usage:
#       HMD_HGNC_Accession.py
#
# Used by:
#
# Notes:
#
# History:
#
# lec	01/04/2005
#	- TR 6456; HGNC
#	- TR 5939; LocusLink->EntrezGene
#
# lec	07/02/2002
#	- TR 3855; add human LocusLink ID column
#
# lec	01/13/98
#	- added comments
#
'''
 
import sys
import os
import string
import db
import reportlib

#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []
cmds.append('select distinct hgncSymbol = m1.symbol, hgncKey = m1._Marker_key, ' + \
            'mgiSymbol = m2.symbol, mgiName = m2.name, mgiKey = m2._Marker_key ' + \
            'into #homology ' + \
            'from HMD_Homology h1, HMD_Homology h2, ' + \
            'HMD_Homology_Marker hm1, HMD_Homology_Marker hm2, ' + \
            'MRK_Marker m1, MRK_Marker m2 ' + \
            'where m1._Organism_key = 2 ' + \
            'and m1._Marker_key = hm1._Marker_key ' + \
            'and hm1._Homology_key = h1._Homology_key ' + \
            'and h1._Class_key = h2._Class_key ' + \
            'and h2._Homology_key = hm2._Homology_key ' + \
            'and hm2._Marker_key = m2._Marker_key ' + \
            'and m2._Organism_key = 1')

cmds.append('create index idx1 on #homology(hgncKey)')
cmds.append('create index idx2 on #homology(mgiKey)')
db.sql(cmds, None)

results = db.sql('select a.accID, a._Object_key from ACC_Accession a, #homology h ' + 
	    'where h.hgncKey = a._Object_key ' + \
	    'and a._MGIType_key = 2 ' + \
	    'and a._LogicalDB_key = 64 ' + \
	    'and a.preferred = 1', 'auto')
hgnc = {}
for r in results:
	hgnc[r['_Object_key']] = r['accID']


results = db.sql('select a.accID, a._Object_key from ACC_Accession a, #homology h ' + 
	    'where h.mgiKey = a._Object_key ' + \
	    'and a._MGIType_key = 2 ' + \
	    'and a.prefixPart = "MGI:" and a._LogicalDB_key = 1 and a.preferred = 1', 'auto')
mgi = {}
for r in results:
	mgi[r['_Object_key']] = r['accID']


results = db.sql('select a.accID, a._Object_key from ACC_Accession a, #homology h ' + 
	    'where h.hgncKey = a._Object_key ' + \
	    'and a._MGIType_key = 2 ' + \
	    'and a._LogicalDB_key = 55', 'auto')
eg = {}
for r in results:
	eg[r['_Object_key']] = r['accID']

results = db.sql('select * from #homology order by mgiSymbol', 'auto')
for r in results:

	if mgi.has_key(r['mgiKey']):

		fp.write(mgi[r['mgiKey']] + reportlib.TAB + \
	         	r['mgiSymbol'] + reportlib.TAB + \
	         	r['mgiName'][0:50] + reportlib.TAB)

		if hgnc.has_key(r['hgncKey']):
			fp.write(hgnc[r['hgncKey']])

		fp.write(reportlib.TAB + \
	         	r['hgncSymbol'] + reportlib.TAB)

		if eg.has_key(r['hgncKey']):
			fp.write(eg[r['hgncKey']])

		fp.write (reportlib.CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
