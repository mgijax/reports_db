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
# lec	09/20/2005
#	- created
#
'''
 
import sys 
import os
import db
import reportlib
import string
import regsub

CRT = reportlib.CRT
TAB = reportlib.TAB

#
# Main
#

fp = reportlib.init(sys.argv[0], fileExt = '.sql.rpt', title = 'MGI cDNA Clone Libraries', outputdir = os.environ['REPORTOUTPUTDIR'])
fp.write(string.ljust('library name', 95))
fp.write('clone collection' + CRT)
fp.write(string.ljust('------------', 95))
fp.write('----------------' + 2*CRT)

db.sql('select name = substring(ps.name, 1, 90), cloneCollection = s.name ' + \
    'into #libraries ' + \
    'from MGI_Set s, MGI_SetMember sm, PRB_Source ps ' + \
    'where s._MGIType_key = 5 ' + \
    'and s.name in ("IMAGE", "NIA", "RIKEN", "RIKEN (FANTOM)") ' + \
    'and s._Set_key = sm._Set_key ' + \
    'and sm._Object_key = ps._Source_key ', None)

results = db.sql('select name, cloneCollection from #libraries order by name', 'auto') 

cloneLib = {}
for r in results:
    library = r['name']
    collection = r['cloneCollection']

    if not cloneLib.has_key(library):
	cloneLib[library] = []
    cloneLib[library].append(collection)

results = db.sql('select distinct name from #libraries order by name', 'auto')
for r in results:
    fp.write(string.ljust(r['name'], 95))
    fp.write(string.join(cloneLib[r['name']], ',') + CRT)

reportlib.trailer(fp)
reportlib.finish_nonps(fp)
