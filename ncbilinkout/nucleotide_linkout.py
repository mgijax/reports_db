#!/usr/local/bin/python

'''
#
# Report:
#	TR 8776 - NCBI LinkOut of Nucleotide -> Marker associations
#	This contains only GenBank (_LogicalDB_key = 9) sequences
#
# Usage:
#       Nucleotide_LinkOut.py
#
# lec	09/30/2008
#	- generation of NCBI LinkOut file (nucleotide-mgd.xml)
#	- copied from MRK_Sequence.py
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
fileName = 'nucleotide-mgd-'

db.useOneConnection(1)
db.set_sqlLogFunction(db.sqlLogAll)

# remove old file names
#os.remove(os.environ['REPORTOUTPUTDIR'] + "/" + fileName + "*")

# deleted sequences

db.sql('select s._Sequence_key into #deleted from SEQ_Sequence s where s._SequenceStatus_key = 316343', None)
db.sql('create index idx1 on #deleted(_Sequence_key)', None)

db.sql('select a.accID, a._LogicalDB_key into #deletedIDs from #deleted d, ACC_Accession a ' + \
    'where d._Sequence_key = a._Object_key ' + \
    'and a._MGIType_key = 19', None)
db.sql('create index idx1 on #deletedIDs(accID)', None)
db.sql('create index idx2 on #deletedIDs(_LogicalDB_key)', None)

# all official/interim mouse markers that have at least one Sequence ID

db.sql('select m._Marker_key, m.symbol ' + \
	'into #markers ' + \
	'from MRK_Marker m ' + \
	'where m._Organism_key = 1 ' + \
	'and m._Marker_Status_key in (1,3) ' + \
	'and exists (select 1 from ACC_Accession a where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 and a._LogicalDB_key in (9))', None)
db.sql('create index idx1 on #markers(_Marker_key)', None)
db.sql('create index idx2 on #markers(symbol)', None)

# MGI ids

results = db.sql('select distinct m._Marker_key, a.accID ' + \
      'from #markers m, ACC_Accession a ' + \
      'where m._Marker_key = a._Object_key ' + \
      'and a._MGIType_key = 2 ' + \
      'and a._LogicalDB_key = 1 ' + \
      'and a.prefixPart = "MGI:" ' + \
      'and a.preferred = 1', 'auto')
mgiID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    mgiID[key] = value

# GenBank ids

results = db.sql('select distinct m._Marker_key, a.accID ' + \
      'from #markers m, ACC_Accession a ' + \
      'where m._Marker_key = a._Object_key ' + \
      'and a._MGIType_key = 2 ' + \
      'and a._LogicalDB_key = 9 ' + \
      'and not exists (select 1 from #deletedIDs d where a.accID = d.accID and a._LogicalDB_key = d._LogicalDB_key)', 'auto')
gbID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if not gbID.has_key(key):
	gbID[key] = []
    gbID[key].append(value)

# process

results = db.sql('select * from #markers order by _Marker_key', 'auto')
fileCounter = 1
count = 1

for r in results:
    key = r['_Marker_key']

    if not gbID.has_key(key) and not rsID.has_key(key):
	    continue

    # create a file

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
    fpLinkOut.write(2* TAB + '<Database>Nucleotide</Database>' + CRT)
    fpLinkOut.write(2* TAB + '<ObjectList>' + CRT)

    if gbID.has_key(key):
	for i in gbID[key]:
	    fpLinkOut.write(3 * TAB + '<Query>' + i + '[pacc]</Query>' + CRT)

    fpLinkOut.write(2* TAB + '</ObjectList>' + CRT)
    fpLinkOut.write(TAB + '</ObjectSelector>' + CRT)
    fpLinkOut.write(TAB + '<ObjectUrl>' + CRT)
    fpLinkOut.write(2* TAB + '<Base>&base;</Base>' + CRT)
    fpLinkOut.write(2* TAB + '<Rule>id=' + mgiID[key] + '</Rule>' + CRT)
    fpLinkOut.write(2* TAB + '<UrlName>' + r['symbol'] + '</UrlName>' + CRT)
    fpLinkOut.write(TAB + '</ObjectUrl>' + CRT)
    fpLinkOut.write('</Link>' + CRT)

    count = count + 1

    if count == 10000:
        fpLinkOut.write('</LinkSet>' + CRT)
        reportlib.finish_nonps(fpLinkOut)
	count = 1
        fileCounter = fileCounter + 1

db.useOneConnection(0)

