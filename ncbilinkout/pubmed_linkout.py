#!/usr/local/bin/python

'''
#
# Report:
#       TR8776 - NCBI LinkOut of PubMed Ids
#
# Usage:
#       PubMed_LinkOut.py
#
# History:
#
# lec   09/30/2008
#       - generation of NCBI LinkOut file (pubmed-mgi.xml)
#
'''
 
import sys
import os
import string
import db
import mgi_utils
import reportlib

CRT = reportlib.CRT
TAB = reportlib.TAB

db.useOneConnection(1)

fpLinkOut = reportlib.init('pubmed-mgd', fileExt = '.xml', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fpLinkOut.write('<!DOCTYPE LinkSet PUBLIC "-//NLM//DTD LinkOut //EN" "LinkOut.dtd"\n[' + CRT)
fpLinkOut.write('<!ENTITY icon "' + os.environ['NCBILINKOUT_ICON'] + '">' + CRT)
fpLinkOut.write('<!ENTITY base "' + os.environ['NCBILINKOUT_BASE'] + '">' + CRT)
fpLinkOut.write(']>' + CRT)

fpLinkOut.write('''<LinkSet>
   <Link>
      <LinkId></LinkId>
         <ProviderId>2002</ProviderId>
         <IconUrl>&icon;</IconUrl>
         <ObjectSelector>
            <Database>PubMed</Database>
            <ObjectList>
''')

# retrieve all PubMed Ids

results = db.sql('select a.accID ' + \
      'from ACC_Accession a ' + \
      'where a._MGIType_key = 1 ' + \
      'and a._LogicalDB_key = 29 ' + \
      'order by a.accID', 'auto')

#results = db.sql('select a.accID, pubMedID = b.accID ' + \
#       'from ACC_Accession a, ACC_Accession b ' + \
#       'where a._MGIType_key = 1 ' + \
#       'and a._LogicalDB_key = 1 ' + \
#       'and a.prefixPart = "MGI:" ' + \
#       'and a._LogicalDB_key = 1 ' + \
#       'and a.preferred = 1 ' + \
#       'and a._Object_key = b._Object_key ' + \
#       'and b._MGIType_key = 1 ' + \
#       'and b._LogicalDB_key = 29 ' + \
#       'order by a.accID', 'auto')

count = 1
for r in results:
    fpLinkOut.write(2 * TAB + '<ObjId>' + r['accID'] + '</ObjId>' + CRT)
    count = count + 1

fpLinkOut.write('''         </ObjectList>
         </ObjectSelector>
         <ObjectUrl>
            <Base>&base;</Base>
            <Rule>id=&lo.id;</Rule>
         </ObjectUrl>
      </Link>
</LinkSet>
''')

reportlib.finish_nonps(fpLinkOut)
db.useOneConnection(0)

