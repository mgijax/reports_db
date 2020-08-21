
'''
#
# MGI_Nonstandard_Strain.py
#
# Report:
#       Tab-delimited file
#       All Strains  w/ Standard = true, Private = false
#	Exclude 'genetic backgrounds'
#	Exclude strainPrefix = yes
#	Exclude  these strains
#	1. 'Not Applicable' (MGI:5649511), 
#	2. 'Not Specified' (MGI:4867032), 
#	3. 'Not Curated' (MGI:5651265), 
#	4. 'Not Loaded' (MGI:5652408), 
#	5. 'Not Resolved' (MGI:5652861)
@
#	Columns:
#	1. MGI ID 
#	2. Strain Name
#	3. Strain Type
#
#	Sort by strain name
#
# Usage:
#       GI_Nonstandard_Strain.py
#
# Notes:
#
# History:
#
# sc	07/11/2016
#	- TR12339
#
'''
 
import sys
import os
import reportlib
import db

db.setTrace()

#
# Main
#

fp1 = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# Retrieve all Strains w/ Standard = true, Private = false

db.sql('''
        select s.strain, s.strainType, a.accID 
        into temporary table strain 
        from PRB_Strain_view s, ACC_Accession a 
        where s.standard = 0
        and s.private = 0 
        and not (s.strain  like 'involves%' 
        or s.strain like 'either:%' 
        or s.strain like 'chimera involves:%' 
        or s.strain like 'mixed%' 
        or s.strain like 'various%')
        and s.geneticbackground = 0
        and s._Strain_key = a._Object_key 
        and a._MGIType_key = 10 
        and a.prefixPart = 'MGI:' 
        and a.preferred = 1
        and a.accID not in  ('MGI:5649511', 'MGI:4867032', 'MGI:5651265', 'MGI:5652408', 'MGI:5652861')

        ''', None)

results = db.sql('select * from strain order by strain', 'auto')
for r in results:

        fp1.write(r['accID'] + reportlib.TAB)
        fp1.write(r['strain'] + reportlib.TAB)
        fp1.write(r['strainType'] + reportlib.CRT)

reportlib.finish_nonps(fp1)
