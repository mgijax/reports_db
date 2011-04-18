#!/usr/local/bin/python

'''
#
# PRB_CloneLibrarycDNA.py
#
# Report:
#       TR 7131
#
# Usage:
#       PRB_CloneLibrarycDNA.py
#
# History:
#
# lec	09/20/2005
#	- created
#
'''
 
import sys 
import os
import db
import reportlib
import string

CRT = reportlib.CRT
TAB = reportlib.TAB

#
# Main
#

fp = reportlib.init(sys.argv[0], fileExt = '.sql.rpt', title = 'MGI cDNA Clone Libraries', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = "MGI")
fp.write(string.ljust('library name', 95))
fp.write(string.ljust('clone collection', 50))
fp.write('clone ids' + CRT)
fp.write(string.ljust('------------', 95))
fp.write(string.ljust('----------------', 50))
fp.write('----------------' + CRT)

db.sql('select _Source_key, name = substring(ps.name, 1, 90), cloneCollection = s.name ' + \
    'into #libraries ' + \
    'from MGI_Set s, MGI_SetMember sm, PRB_Source ps ' + \
    'where s._MGIType_key = 5 ' + \
    'and s.name in ("IMAGE", "NIA", "RIKEN", "RIKEN (FANTOM)") ' + \
    'and s._Set_key = sm._Set_key ' + \
    'and sm._Object_key = ps._Source_key ', None)

accIDs = {}
results = db.sql('select distinct l._Source_key, a.accID, ll.name from #libraries l, ACC_Accession a, ACC_LogicalDB ll ' + \
	'where l._Source_key = a._Object_key ' + \
	'and a._MGIType_key = 5 ' + \
	'and a._LogicalDB_key = ll._LogicalDB_key', 'auto')
for r in results:
    key = r['_Source_key']
    tokens = string.split(r['name'], ' ')
    value = tokens[0] + ':' + string.strip(r['accID'])

    if not accIDs.has_key(key):
	accIDs[key] = []

    accIDs[key].append(value)

results = db.sql('select distinct _Source_key, name, cloneCollection from #libraries order by name', 'auto') 

cloneLib = {}
for r in results:
    library = r['name']
    collection = r['cloneCollection']

    if not cloneLib.has_key(library):
	cloneLib[library] = []
    cloneLib[library].append(collection)

results = db.sql('select distinct _Source_key, name from #libraries order by name', 'auto')
for r in results:
    fp.write(string.ljust(r['name'], 95))
    fp.write(string.ljust(string.join(cloneLib[r['name']], ','), 50))

    if accIDs.has_key(r['_Source_key']):
        fp.write(string.ljust(string.join(accIDs[r['_Source_key']], ','), 50))

    fp.write(CRT)

fp.write(CRT + '(%d rows affected)' % (len(results)) + CRT)
reportlib.finish_nonps(fp)
