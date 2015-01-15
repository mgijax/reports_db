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
# sc	12/12/2014
#	- TR11846 - updated
#
# lec   09/30/2008
#       - generation of NCBI LinkOut file (pubmed-mgi.xml)
#
'''
 
import sys
import os
import string
import mgi_utils
import reportlib
import pg_db
db = pg_db
db.setTrace()
db.setAutoTranslateBE()

CRT = reportlib.CRT
TAB = reportlib.TAB

fileName = 'pubmed-mgd-1'

db.useOneConnection(1)

# remove old file names
os.system('rm -rf ' + os.environ['REPORTOUTPUTDIR'] + "/" + fileName + "*")

# UID file descriptor
fpUid = reportlib.init('pubmed-mgd', fileExt = '.uid', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# xml file descriptor
fpLinkOut = reportlib.init(fileName, fileExt = '.xml', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

count = 1

fpLinkOut.write('<!DOCTYPE LinkSet PUBLIC "-//NLM//DTD LinkOut //EN" "LinkOut.dtd"\n[' + CRT)
fpLinkOut.write('<!ENTITY icon "' + os.environ['NCBILINKOUT_ICON'] + '">' + CRT)
fpLinkOut.write('<!ENTITY base "' + os.environ['NCBILINKOUT_BASE'] + 'reference/' + '">' + CRT)
fpLinkOut.write(']>' + CRT + '<LinkSet>' + CRT)
fpLinkOut.write('<Link>' + CRT)
fpLinkOut.write(TAB + '<LinkId>' + str(count) + '</LinkId>' + CRT)
fpLinkOut.write(TAB + '<ProviderId>2002</ProviderId>' + CRT)
fpLinkOut.write(TAB + '<IconUrl>&icon;</IconUrl>' + CRT)
fpLinkOut.write(TAB + '<ObjectSelector>' + CRT)
fpLinkOut.write(2* TAB + '<Database>PubMed</Database>' + CRT)
fpLinkOut.write(2* TAB + '<ObjectList>' + CRT)
fpLinkOut.write(3* TAB + '<FileName fieldname="uid">pubmed-mgd.uid</FileName>' + CRT)
fpLinkOut.write(2* TAB + '</ObjectList>' + CRT)
fpLinkOut.write(TAB + '</ObjectSelector>' + CRT)
fpLinkOut.write(TAB + '<ObjectUrl>' + CRT)
fpLinkOut.write(2* TAB + '<Base>&base;</Base>' + CRT)
fpLinkOut.write(2* TAB + '<Rule>&lo.id;</Rule>' + CRT)
fpLinkOut.write(TAB + '</ObjectUrl>' + CRT)
fpLinkOut.write('</Link>' + CRT)
fpLinkOut.write('</LinkSet>' + CRT)

# retrieve all PubMed Ids
results = db.sql('select a.accID, b.accID as pubMedID ' + \
       'from ACC_Accession a, ACC_Accession b ' + \
       'where a._MGIType_key = 1 ' + \
       'and a._LogicalDB_key = 1 ' + \
       'and a.prefixPart = "MGI:" ' + \
       'and a._LogicalDB_key = 1 ' + \
       'and a.preferred = 1 ' + \
       'and a._Object_key = b._Object_key ' + \
       'and b._MGIType_key = 1 ' + \
       'and b._LogicalDB_key = 29 ' + \
       'order by a.accID', 'auto')

for r in results:
    fpUid.write(r['pubMedID'] + CRT)

reportlib.finish_nonps(fpLinkOut)
reportlib.finish_nonps(fpUid)
db.useOneConnection(0)

