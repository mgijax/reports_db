#!/usr/local/bin/python

'''
#
# PRB_Strain.py 01/03/2003
#
# Report:
#       Tab-delimited file
#       Public Strains
#
# Usage:
#       PRB_Strain.py
#
# Used by:
#       
# Notes:
#
# History:
#
# lec	01/03/2003
#	- TR 4378
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

# Retrieve all Super Strains

cmds.append('select s._Strain_key, s.strain ' + \
	'into #strains ' + \
	'from PRB_Strain s, VOC_Annot a ' + \
	'where a._AnnotType_key = 1003 ' + \
	'and a._Object_key = s._Strain_key ')

# Retrieve MGI Accession number

cmds.append('select distinct a._Object_key, a.accID ' + \
	'from #strains s, PRB_Strain_Acc_View a ' + \
	'where s._Strain_key = a._Object_key ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.preferred = 1')

cmds.append('select * from #strains order by strain')

results = db.sql(cmds, 'auto')

mgiIDs = {}

for r in results[-2]:
	mgiIDs[r['_Object_key']] = r['accID']
	
for r in results[-1]:
	fp.write(mgiIDs[r['_Strain_key']] + reportlib.TAB)
	fp.write(r['strain'] + reportlib.CRT)

reportlib.finish_nonps(fp)

