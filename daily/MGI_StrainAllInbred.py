#!/usr/local/bin/python

'''
#
# MGI_StrainAll.py 04/01/2004
#
# Report:
#       Tab-delimited file
#	Official Strain Nomenclature for All Inbred Strains
#
# Usage:
#       MGI_StrainAllInbred.py
#
# Used by:
#
# Notes:
#
# History:
#
# lec	04/01/2004
#	- TR 5638
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

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []

# Retrieve all standard, public Strains

cmds.append('select _Strain_key, strain into #strains from PRB_Strain where standard = 1 and private = 0')

# Retrieve MGI Accession number

cmds.append('select distinct a._Object_key, a.accID ' + \
	'from #strains s, PRB_Strain_Acc_View a ' + \
	'where s._Strain_key = a._Object_key ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1')

cmds.append('select * from #strains order by strain')

results = db.sql(cmds, 'auto')

mgiIDs = {}

for r in results[-2]:
	mgiIDs[r['_Object_key']] = r['accID']
	
for r in results[-1]:
	fp.write(mgiIDs[r['_Strain_key']] + reportlib.TAB + \
	         r['strain'] + reportlib.CRT)

reportlib.finish_nonps(fp)

