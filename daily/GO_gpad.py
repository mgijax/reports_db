#!/usr/local/bin/python

'''
#
# GO_gpad.py
#
# Report:
#       see: wiki/mediawiki/index.php/sw:GOload 
#	contains link to GO/GPAD format
#
# Output format:
#
#   1. DB			MGI
#   2. DB Object ID		12345
#   3. Qualifier		
#   4. GO ID			GO:xxxx
#   5. DB:Reference(s)		PMID:xxxx
#   6. Evidence code		ECO:xxxx
#   7. With (or)From		(optional)
#   8. Interacting taxon ID	(optional)
#   9. Date			YYYYMMDD
#   10. Assigned by
#   11. Annotation Extension	(optional)
#   12. Annotation Properties	(optional)
#
# History:
#
# 06/27/2016	lec
#	- TR12349/12345/GPAD/GPI
#
'''

import sys
import os
import mgi_utils
import reportlib
import db

db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE(False)

TAB = reportlib.TAB
CRT = reportlib.CRT

db.useOneConnection(1)

fp = reportlib.init('mgi_association', fileExt = '.gpad', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('!gpa-version: 1.1\n') 
fp.write('!date: %s $\n' % (mgi_utils.date("%m/%d/%Y")))
fp.write('!\n')
fp.write('! from Mouse Genome Database (MGD) & Gene Expression Database (GXD)\n')
fp.write('!\n')

db.sql('''select distinct a._Object_key, a._Term_key, 
	    ma.prefixPart, ma.numericPart, ma.accID as markerID,
	    t1.term, ta.accID as termID, 
	    t2.term as qualifier,
            e.inferredFrom,
            e._Refs_key, e._ModifiedBy_key,
	    to_char(e.modification_date, 'yyyyMMdd') as moddate,
	    u.login as assignedBy
        into temporary table gomarker1 
        from VOC_Annot a,  
             ACC_Accession ta, 
             VOC_Term t1,  
             VOC_Evidence e,  
             VOC_Term t2,
	     ACC_Accession ma,
	     MGI_User u
        where a._AnnotType_key = 1000 
        and a._Annot_key = e._Annot_key 
        and a._Term_key = t1._Term_key 
        and a._Term_key = ta._Object_key 
        and ta._MGIType_key = 13  
        and ta.preferred = 1 
        and a._Qualifier_key = t2._Term_key 
	and a._Object_key = ma._Object_key
	and ma._MGIType_key = 2
	and ma._LogicalDB_key = 1
	and ma.prefixPart = 'MGI:'
	and ma.preferred = 1
	and e._ModifiedBy_key = u._User_key
and a._Object_key = 13433
        ''', None)

db.sql('create index gomarker1_idx1 on gomarker1(_Object_key)', None)

results = db.sql('select * from gomarker1', 'auto')

for r in results:
	#   1. DB			MGI
	fp.write(r['prefixPart'] + TAB)

	#   2. DB Object ID		12345
	fp.write(r['markerID'] + TAB)

	#   3. Qualifier		
        fp.write('(3)' + TAB)

	#   4. GO ID			GO:xxxx
	fp.write(r['termID'] + TAB)

	#   5. DB:Reference(s)		PMID:xxxx
        fp.write('(5)' + TAB)

	#   6. Evidence code		ECO:xxxx
        fp.write('(6)' + TAB)

	#   7. With (or)From		(optional)
        fp.write(mgi_utils.prvalue(r['inferredFrom']) + TAB)

	#   8. Interacting taxon ID	(optional)
        fp.write('(8)' + TAB)

	#   9. Date			YYYYMMDD
        fp.write(r['moddate'] + TAB)

	#   10. Assigned by
        fp.write(r['assignedBy'] + TAB)

	#   11. Annotation Extension	(optional)
        fp.write('(11)' + TAB)

	#   12. Annotation Properties	(optional)
        fp.write('(12)' + CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)

