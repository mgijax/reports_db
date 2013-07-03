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
import mgi_utils
import reportlib
import pg_db
db = pg_db
db.setTrace()
db.setAutoTranslateBE()

CRT = reportlib.CRT
TAB = reportlib.TAB

maxfileCounter = int(os.environ['NCBILINKOUT_COUNT'])
fileName = 'pubmed-mgd-'

db.useOneConnection(1)

# remove old file names
os.system('rm -rf ' + os.environ['REPORTOUTPUTDIR'] + "/" + fileName + "*")

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

fileCounter = 1
count = 1
for r in results:

    if count == 1:
        newfile = fileName + str(fileCounter)
        fpLinkOut = reportlib.init(newfile, fileExt = '.xml', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
        fpLinkOut.write('<!DOCTYPE LinkSet PUBLIC "-//NLM//DTD LinkOut //EN" "LinkOut.dtd"\n[' + CRT)
        fpLinkOut.write('<!ENTITY icon "' + os.environ['NCBILINKOUT_ICON'] + '">' + CRT)
        fpLinkOut.write('<!ENTITY base "' + os.environ['NCBILINKOUT_BASE'] + '">' + CRT)
        fpLinkOut.write(']>' + CRT + '<LinkSet>' + CRT)

    fpLinkOut.write('<Link>' + CRT)
    fpLinkOut.write(TAB + '<LinkId>' + str(count) + '</LinkId>' + CRT)
    fpLinkOut.write(TAB + '<ProviderId>2002</ProviderId>' + CRT)
    fpLinkOut.write(TAB + '<IconUrl>&icon;</IconUrl>' + CRT)
    fpLinkOut.write(TAB + '<ObjectSelector>' + CRT)
    fpLinkOut.write(2* TAB + '<Database>PubMed</Database>' + CRT)
    fpLinkOut.write(2* TAB + '<ObjectList>' + CRT)
    fpLinkOut.write(3* TAB + '<ObjId>' + r['pubMedID'] + '</ObjId>' + CRT)
    fpLinkOut.write(2* TAB + '</ObjectList>' + CRT)
    fpLinkOut.write(TAB + '</ObjectSelector>' + CRT)
    fpLinkOut.write(TAB + '<ObjectUrl>' + CRT)
    fpLinkOut.write(2* TAB + '<Base>http://www.informatics.jax.org/searches/accession_report.cgi?id=' 
		+ r['accID']
		+ '</Base>' + CRT)
    fpLinkOut.write(TAB + '</ObjectUrl>' + CRT)
    fpLinkOut.write('</Link>' + CRT)

    count = count + 1

    # skip to a new file every 'maxfileCounter' times...

    if count == maxfileCounter:
        fpLinkOut.write('</LinkSet>' + CRT)
        reportlib.finish_nonps(fpLinkOut)
	count = 1
        fileCounter = fileCounter + 1

fpLinkOut.write('</LinkSet>' + CRT)
reportlib.finish_nonps(fpLinkOut)
db.useOneConnection(0)

