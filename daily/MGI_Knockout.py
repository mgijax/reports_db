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
#		Allele Types: targeted (knockout), targeted (reporter), targeted (other)
#
#		the allele's genotype that is homozygous for the mutation
#		OR
#		the allele has no genotypes
#
#		field 1: MGI Gene ID
#		field 2: Gene symbol
#		field 3: Gene name
#		field 4: Allele Types
#		field 5: Allele Symbols
#		field 6: Representative Transcript
#		field 7: Representative Gene Model
#		field 8: IMSR strain names (|-delimited)
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
import mgi_html
import table

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE
IMSR = os.environ['IMSR_DBNAME']

def printMarker(r):

    s = '<tr>' + \
	'<td>' + r['accID'] + '</td>\n' + \
        '<td>%s%s%s</td>\n' % (reportlib.create_accession_anchor(r['accID']), r['symbol'], reportlib.close_accession_anchor()) + \
	'<td>' + r['name'] + '</td>\n'

    if repgen.has_key(r['_Marker_key']):
	s = s + '<td>%s</td>' % (repgen[r['_Marker_key']])
    else:
	s = s + '<td>&nbsp;</td>'

    if reptran.has_key(r['_Marker_key']):
	s = s + '<td>%s</td>' % (reptran[r['_Marker_key']])
    else:
	s = s + '<td>&nbsp;</td>'

    s = s + '<td nowrap>%s</td>' % (string.join(alleleTypes[r['_Marker_key']], '<br>'))
    return s

def printAllele(a):

    symbol = regsub.gsub('<', '<sup>', a['symbol'])
    symbol = regsub.gsub('>', '</sup>', symbol)

    s = '%s%s%s<br>' % (reportlib.create_accession_anchor(a['accID']), symbol, reportlib.close_accession_anchor())

    return s

def printHeader(fp, title):

    fp.write('</pre>\n')
    fp.write('<H2>%s</H2>' % (title))
    fp.write('<TABLE BORDER=3 WIDTH=100%>')
    fp.write('<th align=left valign=top>MGI Gene ID</th>')
    fp.write('<th align=left valign=top>Gene Symbol</th>')
    fp.write('<th align=left valign=top>Gene Name</th>')
    fp.write('<th align=left valign=top>Ensembl ID</th>')
    fp.write('<th align=left valign=top>RefSeq/Transcript</th>')
    fp.write('<th align=left valign=top>Allele Types</th>')
    fp.write('<th align=left valign=top>Alleles/Phenotype/Disease</th>')
    fp.write('<th align=left valign=top>IMSR Strain</th>')
    return

#
# Main
#

fullreport = 'MGI_Knockout_Full'
publicreport = 'MGI_Knockout_Public'
notpublicreport = 'MGI_Knockout_NotPublic'

fp1 = reportlib.init(fullreport, printHeading = 0, outputdir = os.environ['REPORTOUTPUTDIR'], isHTML = 1)
fp2 = reportlib.init(publicreport, printHeading = 0, outputdir = os.environ['REPORTOUTPUTDIR'], isHTML = 1)
fp3 = reportlib.init(notpublicreport, printHeading = 0, outputdir = os.environ['REPORTOUTPUTDIR'], isHTML = 1)

printHeader(fp1, 'MGI Full KnockOut Report')
printHeader(fp2, 'MGI Public KnockOut Report')
printHeader(fp3, 'MGI Non-Public KnockOut Report')

#
# select alleles
#

db.sql('select a._Marker_key, a._Allele_key, a.symbol, t.term, aa.accID ' + \
	'into #knockouts ' + \
	'from ALL_Allele a, ACC_Accession aa, VOC_Term t ' + \
	'where a._Allele_Status_key = 847114 ' + \
	'and a._Allele_Type_key in (847116, 847119, 847120) ' + \
	'and a.name not like "%Lexicon%" ' + \
	'and a.name not like "%Deltagen%" ' + \
	'and a._Allele_key = aa._Object_key ' + \
	'and aa._MGIType_key = 11 ' + \
	'and aa._LogicalDB_key = 1 ' + \
	'and aa.prefixPart = "MGI:" ' + \
	'and aa.preferred = 1 ' + \
	'and a._Allele_Type_key = t._Term_key ', None)

db.sql('create index idx1 on #knockouts(_Marker_key)', None)
db.sql('create index idx2 on #knockouts(_Allele_key)', None)

#
# select unique set of markers
#

db.sql('select distinct m._Marker_key, m.symbol, name = substring(m.name,1,75), ma.accID ' + \
	'into #markers ' + \
	'from #knockouts k, MRK_Marker m, ACC_Accession ma ' + \
	'where k._Marker_key = m._Marker_key ' + \
	'and m._Marker_key = ma._Object_key ' + \
	'and ma._MGIType_key = 2 ' + \
	'and ma._LogicalDB_key = 1 ' + \
	'and ma.prefixPart = "MGI:" ' + \
	'and ma.preferred = 1 ', None)

db.sql('create index idx1 on #markers(symbol)', None)

#
# select those alleles that are in IMSR and that occur "alone" in a strain
#
results = db.sql('select distinct m._Marker_key, ac.accID, ls.label ' + \
	'from #markers m, #knockouts k, %s..Accession ac, %s..Label ls, %s..SGAAssoc sga ' % (IMSR, IMSR, IMSR) + \
	'where m._Marker_key = k._Marker_key ' + \
	'and k.accID = ac.accID ' + \
	'and sga._Strain_key = ls._Object_key ' + \
	'and ls._IMSRType_key = 1 ' + \
	'and ls.labelType = "N" ' + \
	'and sga._Allele_key = ac._Object_key ' + \
	'and ac._IMSRType_key = 3 ' + \
	'and 1 = (select count(distinct _Allele_key) ' + \
	'from %s..SGAAssoc ' % (IMSR) + \
	'where sga._Strain_key = _Strain_key)', 'auto')
imsr = {}
for r in results:
    key = r['_Marker_key']
    value = regsub.gsub('<', '<sup>', r['label'])
    value = regsub.gsub('>', '</sup>', value)
    value = '%s%s%s' % (reportlib.create_imsrstrain_anchor(r['label']), value, reportlib.close_accession_anchor())
    if not imsr.has_key(key):
	imsr[key] = []
    if value not in imsr[key]:
        imsr[key].append(value)

#
# cache gene/allele data
#

alleles = {}
results = db.sql('select * from #knockouts', 'auto')
for r in results:
    key = r['_Marker_key']
    value = r
    if not alleles.has_key(key):
	alleles[key] = []
    alleles[key].append(value)

#
# unique set of Allele Types per Gene
#

alleleTypes = {}
results = db.sql('select distinct m._Marker_key, a.term from #markers m, #knockouts a ' + 
	'where m._Marker_key = a._Marker_key', 'auto')
for r in results:
    key = r['_Marker_key']
    value = r['term']
    if not alleleTypes.has_key(key):
	alleleTypes[key] = []
    alleleTypes[key].append(value)

#
# representative genomic for genes
#
repgen = {}
results = db.sql('select m._Marker_key, a.accID ' + \
	'from #markers m, SEQ_Marker_Cache smc, ACC_Accession a ' + \
	'where m._Marker_key = smc._Marker_key ' + \
	'and smc._LogicalDB_key = 60 ' + \
	'and smc._Sequence_key = a._Object_key ' + \
	'and a._MGIType_key = 19 ' + \
	'and a.preferred = 1', 'auto')
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    repgen[key] = value

#
# representative transcript for genes
#
reptran = {}
results = db.sql('select m._Marker_key, a.accID ' + \
	'from #markers m, SEQ_Marker_Cache smc, ACC_Accession a ' + \
	'where m._Marker_key = smc._Marker_key ' + \
	'and smc._Qualifier_key = 615420 ' + \
	'and smc._Sequence_key = a._Object_key ' + \
	'and a._MGIType_key = 19 ' + \
	'and a.preferred = 1', 'auto')
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    reptran[key] = value

#
# process results
#
results = db.sql('select * from #markers order by symbol', 'auto')

for r in results:

    mstr = ''
    astr = ''
    key = r['_Marker_key']
    mstr = printMarker(r)

    # for each allele....

    for a in alleles[key]:
	astr = astr + printAllele(a)

    if imsr.has_key(key):
	istr = '<td>%s</td>' % (string.join(imsr[key], '<br>'))
    else:
	istr = '<td>&nbsp;</td>'

    toPrint = '%s<td>%s</td>%s' % (mstr, astr, istr)

    fp1.write(toPrint)

    if imsr.has_key(key):
	fp2.write(toPrint)
    else:
	fp3.write(toPrint)

fp1.write('</TABLE>')
fp1.write('<pre>')

reportlib.finish_nonps(fp1, isHTML = 1)	# non-postscript file
reportlib.finish_nonps(fp2, isHTML = 1)	# non-postscript file
reportlib.finish_nonps(fp3, isHTML = 1)	# non-postscript file

