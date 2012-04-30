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
# History:
#
# lec	07/15/2005
#	- added Mouse EntrezGene ID per Michael at HUGO
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


#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

db.sql('''
	select distinct m1.symbol as hgncSymbol, 
			m1._Marker_key as hgncKey, 
                        m2.symbol as mgiSymbol, 
			m2.name as mgiName, 
			m2._Marker_key as mgiKey
        into #homology 
        from MRK_Homology_Cache h1, MRK_Homology_Cache h2, 
        MRK_Marker m1, MRK_Marker m2 
        where h1._Organism_key = 2 
        and h1._Class_key = h2._Class_key 
        and h2._Organism_key = 1 
	and h1._Marker_key = m1._Marker_key 
	and h2._Marker_key = m2._Marker_key
	''', None)

db.sql('create index idx1 on #homology(hgncKey)', None)
db.sql('create index idx2 on #homology(mgiKey)', None)

# human HGNC
results = db.sql('''
	select a.accID, a._Object_key 
	from ACC_Accession a, #homology h 
	where h.hgncKey = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 64 
	and a.preferred = 1
	''', 'auto')
hgnc = {}
for r in results:
	hgnc[r['_Object_key']] = r['accID']

# mouse MGI
results = db.sql('''
	select a.accID, a._Object_key 
	from ACC_Accession a, #homology h 
	where h.mgiKey = a._Object_key 
	and a._MGIType_key = 2 
	and a.prefixPart = 'MGI:' 
	and a._LogicalDB_key = 1 
	and a.preferred = 1
	''', 'auto')
mgi = {}
for r in results:
	mgi[r['_Object_key']] = r['accID']

# mouse EntrezGene
results = db.sql('''
	select a.accID, a._Object_key 
	from ACC_Accession a, #homology h 
	where h.mgiKey = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 55
	''', 'auto')
meg = {}
for r in results:
	meg[r['_Object_key']] = r['accID']

# human EntrezGene
results = db.sql('''
	select a.accID, a._Object_key 
	from ACC_Accession a, #homology h 
	where h.hgncKey = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 55
	''', 'auto')
heg = {}
for r in results:
	heg[r['_Object_key']] = r['accID']

results = db.sql('select * from #homology order by mgiSymbol', 'auto')
for r in results:

	if mgi.has_key(r['mgiKey']):

		fp.write(mgi[r['mgiKey']] + reportlib.TAB + \
	         	r['mgiSymbol'] + reportlib.TAB + \
	         	r['mgiName'] + reportlib.TAB)

		if meg.has_key(r['mgiKey']):
			fp.write(meg[r['mgiKey']])
		fp.write(reportlib.TAB)

		if hgnc.has_key(r['hgncKey']):
			fp.write(hgnc[r['hgncKey']])
		fp.write(reportlib.TAB)

		fp.write(r['hgncSymbol'] + reportlib.TAB)

		if heg.has_key(r['hgncKey']):
			fp.write(heg[r['hgncKey']])

		fp.write (reportlib.CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
