#!/usr/local/bin/python

'''
#
# MGI_Knockout.py
#
# Reports:
#
#       For NIH
#
#	1.  Tab-delimited report where a row in the report represents:
#
#		a targeted (knockout) allele's genotype that is homozygous for the mutation
#		OR
#		a targeted (knockout) allele with no genotypes
#
#		field 1: MGI Marker ID
#		field 2: Marker symbol
#		field 3: Marker name
#		field 4: MGI Allele ID
#		field 5: Allele symbol
#		field 6: Allele name
#		field 7: ES cell line where KO was made
#		field 8: ES cell line strain
# 		field 9: genetic background (genotype strain)
#		field 10: IMSR strain names (|-delimited)
#
#	2.  Tab-delimited report where a row in the report represents:
#
#		(a targeted (knockout) allele's genotype is homozygous for the mutation
#		OR
#		a targeted (knockout) allele with no genotypes)
#		AND
#		at least one of the gene's knockout alleles is in IMSR
#
#	3.  Tab-delimited report where a row in the report represents:
#
#		(a targeted (knockout) allele's genotype is homozygous for the mutation
#		OR
#		a targeted (knockout) allele with no genotypes)
#		AND
#		none of the gene's knockout alleles is in IMSR
#
# Usage:
#       MGI_Knockout.py
#
# Notes:
#	- all reports use mgireport directory for output file
#	- all reports use db default of public login
#	- all reports use server/database default of environment
#	- use lowercase for all SQL commands (i.e. select not SELECT)
#	- all public SQL reports require the header and footer
#	- all private SQL reports require the header
#
# History:
#
# lec	01/23/2006
#	- created, TR 7418
#
'''
 
import sys 
import os
import db
import regsub
import string
import reportlib

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE
IMSR = os.environ['IMSR_DBNAME']

legend = '''
#
#	field 1: MGI Marker ID
#	field 2: Marker symbol
#	field 3: Marker name
#	field 4: MGI Allele ID
#	field 5: Allele symbol
#	field 6: Allele name
#	field 7: ES cell line where KO was made
#	field 8: ES cell line strain
# 	field 9: genetic background (genotype strain)
#	field 10: IMSR strain name
#
'''

def printAllele():

    s = r['markerID'] + TAB + \
	r['markerSym'] + TAB + \
	r['markerName'] + TAB + \
	r['alleleID'] + TAB + \
	r['alleleSym'] + TAB + \
	r['alleleName'] + TAB + \
	r['cellLine'] + TAB + \
	r['strain'] + TAB

    return s

#
# Main
#

fullreport = 'MGI_Knockout_Full'
publicreport = 'MGI_Knockout_Public'
notpublicreport = 'MGI_Knockout_NotPublic'

fp1 = reportlib.init(fullreport, outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)
fp2 = reportlib.init(publicreport, outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)
fp3 = reportlib.init(notpublicreport, outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

fp1.write(legend)
fp2.write(legend)
fp3.write(legend)

#
# select all targeted (knockout) alleles
#

db.sql('select a._Allele_key, alleleSym = a.symbol, alleleName = a.name, alleleID = aa.accID, ' + \
	'markerSym = m.symbol, markerName = m.name, markerID = ma.accID, ' + \
	'cl.cellLine, s.strain ' + \
	'into #knockouts ' + \
	'from ALL_Allele a, ACC_Accession aa, MRK_Marker m, ACC_Accession ma, ALL_CellLine cl, PRB_Strain s ' + \
	'where a._Allele_Type_key = 847116 ' + \
	'and a._Allele_key = aa._Object_key ' + \
	'and aa._MGIType_key = 11 ' + \
	'and aa._LogicalDB_key = 1 ' + \
	'and aa.prefixPart = "MGI:" ' + \
	'and aa.preferred = 1 ' + \
	'and a._Marker_key = m._Marker_key ' + \
	'and m._Marker_key = ma._Object_key ' + \
	'and ma._MGIType_key = 2 ' + \
	'and ma._LogicalDB_key = 1 ' + \
	'and ma.prefixPart = "MGI:" ' + \
	'and ma.preferred = 1 ' + \
	'and a._ESCellLine_key = cl._CellLine_key ' + \
	'and a._Strain_key = s._Strain_key', None)

db.sql('create index idx1 on #knockouts(_Allele_key)', None)
db.sql('create index idx2 on #knockouts(markerSym)', None)

#
# select those that have no genotypes
#
results = db.sql('select k._Allele_key from #knockouts k ' + \
	'where not exists (select 1 from GXD_AlleleGenotype g where k._Allele_key = g._Allele_key)', 'auto')
noGeno = []
for r in results:
    key = r['_Allele_key']
    noGeno.append(key)
    
#
# select genotypes that are homozygous for the mutation
#
results = db.sql('select distinct k._Allele_key, s.strain ' + \
	'from #knockouts k, GXD_Genotype g, GXD_AllelePair p, PRB_Strain s ' + \
	'where k._Allele_key = p._Allele_key_1 ' + \
	'and k._Allele_key = p._Allele_key_2 ' + \
	'and p._Genotype_key = g._Genotype_key ' + \
	'and g._Strain_key = s._Strain_key', 'auto')
geno = {}
for r in results:
    key = r['_Allele_key']
    value = r['strain']
    if not geno.has_key(key):
        geno[key] = []
    geno[key].append(value)
    
#
# select those alleles that are in IMSR and that occur "alone" in a strain
#
results = db.sql('select distinct ac.accID, ls.label ' + \
	'from %s..Accession ac, %s..Label ls, %s..SGAAssoc sga ' % (IMSR, IMSR, IMSR) + \
	'where sga._Strain_key = ls._Object_key ' + \
	'and ls._IMSRType_key = 1 ' + \
	'and ls.labelType = "N" ' + \
	'and sga._Allele_key = ac._Object_key ' + \
	'and ac._IMSRType_key = 3 ' + \
	'and 1 = (select count(distinct _Allele_key) ' + \
	'from %s..SGAAssoc ' % (IMSR) + \
	'where sga._Strain_key = _Strain_key)', 'auto')
imsr = {}
for r in results:
    key = r['accID']
    value = r['label']
    if not imsr.has_key(key):
	imsr[key] = []
    imsr[key].append(value)

#
# process results
#
results = db.sql('select * from #knockouts order by markerSym', 'auto')

for r in results:

    printRecord = ''
    key = r['_Allele_key']
    imsrKey = r['alleleID']

    if imsr.has_key(imsrKey):
        printIMSR = string.join(imsr[imsrKey], '|')
    else:
	printIMSR = TAB

    # record w/ genotype

    if geno.has_key(key):
	for g in geno[key]:
	    printRecord = printAllele() + g + TAB + printIMSR + CRT

    # record w/out genotype

    if key in noGeno:
	printRecord = printAllele() + TAB + printIMSR + CRT

    # print record

    fp1.write(printRecord)

    if imsr.has_key(imsrKey):
	fp2.write(printRecord)
    else:
	fp3.write(printRecord)

reportlib.finish_nonps(fp1)	# non-postscript file
reportlib.finish_nonps(fp2)	# non-postscript file
reportlib.finish_nonps(fp3)	# non-postscript file

