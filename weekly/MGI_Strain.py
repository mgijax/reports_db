#!/usr/local/bin/python

'''
#
# MGI_Strain.py
#
# Report:
#       Tab-delimited file
#       All Strains  w/ Standard = true, Private = false
#	Dislay fields: MGI ID, Strain Name, Strain Types
#
#	Report A: sorted by alpha by strain name
#
# Usage:
#       MGI_Strain.py
#
# Used by:
#       Internal Report
#
# Notes:
#
# History:
#
# lec	04/24/2012
#	- TR11035/postgres cleanup/remove report b (sorted by strain type)
#
# lec	05/06/2008
#	- TR 8511
#
'''
 
import sys
import os
import reportlib
import db

db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE()

#
# Main
#

fp1 = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# Retrieve all Strains w/ Standard = true, Private = false

db.sql('''
	select s.strain, s.strainType, a.accID 
	into #strain 
	from PRB_Strain_View s, ACC_Accession a 
	where s.standard = 1 
	and s.private = 0 
	and s._Strain_key = a._Object_key 
	and a._MGIType_key = 10 
	and a.prefixPart = 'MGI:' 
	and a.preferred = 1
	''', None)

results = db.sql('select * from #strain order by strain', 'auto')
for r in results:

	fp1.write(r['accID'] + reportlib.TAB)
	fp1.write(r['strain'] + reportlib.TAB)
	fp1.write(r['strainType'] + reportlib.CRT)

reportlib.finish_nonps(fp1)

