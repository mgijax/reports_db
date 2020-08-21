'''
#
# Report:
#       TR8776 - NCBI LinkOut of Protein (UniProt) Ids
#
# Usage:
#       Protein_LinkOut.py
#
# History:
#
# sc    03/21/20 python 3 upgrade
#
# lec   09/30/2008
#       - generation of NCBI LinkOut file (protein-mgi.xml)
#       - copied from MRK_SwissProt.py
#
'''
 
import sys
import os
import string
import mgi_utils
import reportlib
import db

db.setTrace()

CRT = reportlib.CRT
TAB = reportlib.TAB

db.useOneConnection(1)

fpLinkOut = reportlib.init('protein-mgd', fileExt = '.xml', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fpLinkOut.write('<!DOCTYPE LinkSet PUBLIC "-//NLM//DTD LinkOut //EN" "LinkOut.dtd"\n[' + CRT)
fpLinkOut.write('<!ENTITY icon "' + os.environ['NCBILINKOUT_ICON'] + '">' + CRT)
fpLinkOut.write('<!ENTITY base "' + os.environ['NCBILINKOUT_BASE_MARKER'] + '">' + CRT)
fpLinkOut.write(']>' + CRT)
fpLinkOut.write('''<LinkSet>
   <Link>
      <LinkId></LinkId>
         <ProviderId>2002</ProviderId>
         <IconUrl>&icon;</IconUrl>
         <ObjectSelector>
            <Database>Protein</Database>
            <ObjectList>
''')

# retrieve all UniProt ID associated with mouse markers

sequence = {}

results = db.sql('select a.accID ' + \
      'from MRK_Marker m, ACC_Accession a ' + \
      'where m._Organism_key = 1 ' + \
      'and m._Marker_key = a._Object_key ' + \
      'and a._MGIType_key = 2 ' + \
      'and a._LogicalDB_key = 13 ' + \
      'order by a.accID', 'auto')

count = 1
for r in results:
    fpLinkOut.write(2 * TAB + '<Query>' + r['accID'] + '[pacc]</Query>' + CRT)
    count = count + 1

fpLinkOut.write('''         </ObjectList>
         </ObjectSelector>
         <ObjectUrl>
            <Base>&base;</Base>
            <Rule>&lo.pacc;</Rule>
         </ObjectUrl>
      </Link>
</LinkSet>
''')

reportlib.finish_nonps(fpLinkOut)
db.useOneConnection(0)

