'''
#
# Report:
#       Tab-delimited file
#
#       1. MGI Strain
#       2. Literature Strain
#
# Usage:
#       MGI_StrainTranslations.py
#
# History:
#
# sc    03/21/20 python 3 upgrade
#
# lec   08/02/2005
#       - Donna Maglott
#
'''
 
import sys 
import os
import reportlib
import db

db.setTrace()

CRT = reportlib.CRT
TAB = reportlib.TAB

#
# Main
#

fp = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])

results = db.sql('''
select s.strain, t.badName 
from MGI_Translation t, PRB_Strain s 
where t._TranslationType_key = 1007 
and t._Object_key = s._Strain_key
''', 'auto')

for r in results:
        fp.write(r['strain'] + TAB + r['badName'] + CRT)

reportlib.finish_nonps(fp)      # non-postscript file
