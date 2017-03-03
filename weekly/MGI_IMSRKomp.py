#!/usr/local/bin/python

'''
#
# TR11615/IMSR/KOMP report
#
# Input:
# 	IMSR data file: 
#	/data/downloads/www.findmice.org/report/kompCounts.txt
#
# Output:
#	1: gene symbol
#	2: mgi id
#	3: imsr strains = field 3 from the kompCounts.txt
#	4: other mutants = total number of alleles for the gene
#
# lec   10/24/2014
#       - TR11750/postres complient
#
'''
 
import sys 
import os
import reportlib
import db

db.setTrace()

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB

#
# Main
#

db.useOneConnection(1)

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('Marker Symbol\tMGI Marker Accession ID\tTotal number of IMSR Strains\tTotal number of Alleles\n')

inFile = open('/data/downloads/www.findmice.org/report/kompCounts.txt', 'rU')

#
# field 1 : MGI id (allele or marker)
# field 2 : type
# field 3 : count
#
# if field 2 == 'MRK:UN', then field 1 == Marker MGI id
#
imsrSet = {}
for line in inFile.readlines():
	tokens = line[:-1].split('\t')
	if tokens[1] == 'MRK:UN':
		imsrSet[tokens[0]] = []
		imsrSet[tokens[0]].append(tokens[2])

results = db.sql('''
	select m._Marker_key, m.symbol, ma.accID, count(al._Allele_key) as mcount
        from MRK_Marker m, ALL_Allele al, ACC_Accession ma
        where al._Allele_Status_key in (847114, 3983021)
        and al.isWildType = 0
        and al._Marker_key = m._Marker_key
        and al._Marker_key = ma._Object_key
        and ma._MGIType_key = 2
        and ma._LogicalDB_key = 1
        and ma.prefixPart = 'MGI:'
        and ma.preferred = 1
        group by m._Marker_key, m.symbol, ma.accID
	''', 'auto')
for r in results:

	acount = 0
	if r['accID'] in imsrSet:
		acount = imsrSet[r['accID']][0]

	fp.write(r['symbol'] + TAB)
	fp.write(r['accID'] + TAB)
	fp.write(str(acount) + TAB)
	fp.write(str(r['mcount']) + CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)

