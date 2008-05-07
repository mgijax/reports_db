#!/usr/local/bin/python

'''
#
# MGI_Strain_Strain_Standard.py
#
# Report:
#       Tab-delimited file
#       All Strains  w/ Standard = true, Private = false
#	Dislay fields: MGI ID, Strain Name, Strain Types
#
#	Report A: sorted by alpha by strain name
#	Report B: sorted by strain type, then by strain name
#
# Usage:
#       MGI_Strain_Standard.py
#
# Used by:
#       Internal Report
#
# Notes:
#
# History:
#
# lec	05/06/2008
#	- TR 8511
#
'''
 
import sys
import os
import string
import db
import mgi_utils
import reportlib

#
# Main
#

fp1 = reportlib.init('MGI_Strain_Standard_A', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = "MGI")
fp2 = reportlib.init('MGI_Strain_Standard_B', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = "MGI")

# Retrieve all Strains w/ Standard = true, Private = false

db.sql('select s.strain, s.strainType, a.accID ' + \
	'into #strain ' + \
	'from PRB_Strain_View s, ACC_Accession a ' + \
	'where s.standard = 1 ' + \
	'and s.private = 0 ' + \
	'and s._Strain_key = a._Object_key ' + \
	'and a._MGIType_key = 10 ' + \
	'and a.prefixPart = "MGI:" ' +
	'and a.preferred = 1', None)

results = db.sql('select * from #strain order by strain', 'auto')
for r in results:

	fp1.write(r['accID'] + reportlib.TAB)
	fp1.write(r['strain'] + reportlib.TAB)
	fp1.write(r['strainType'] + reportlib.CRT)

results = db.sql('select * from #strain order by strainType, strain', 'auto')
for r in results:

	fp2.write(r['accID'] + reportlib.TAB)
	fp2.write(r['strain'] + reportlib.TAB)
	fp2.write(r['strainType'] + reportlib.CRT)

reportlib.finish_nonps(fp1)
reportlib.finish_nonps(fp2)

