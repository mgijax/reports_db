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
import mgi_html
import table

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE
IMSR = os.environ['IMSR_DBNAME']

heading = ['MGI Marker ID', 'Marker Symbol', 'Marker Name', 
	   'MGI Allele ID', 'Allele Symbol', 'Allele Name',
	   'ES Cell Line', 'ES Cell Line Strain', 'Genetic Background', 'IMSR Strain']

def printAllele():

    alleleSymbol = regsub.gsub('<', '&lt;', r['alleleSym'])
    alleleSymbol = regsub.gsub('>', '&gt;', alleleSymbol)

    strain = regsub.gsub('<', '&lt;', r['strain'])
    strain = regsub.gsub('>', '&gt;', strain)

    s = '<tr>' + \
        '<td>%s%s%s</td>\n' % (reportlib.create_accession_anchor(r['markerID']), r['markerID'], reportlib.close_accession_anchor()) + \
	'<td>' + r['markerSym'] + '</td>\n' + \
	'<td>' + r['markerName'] + '</td>\n' + \
        '<td>%s%s%s</td>' % (reportlib.create_accession_anchor(r['alleleID']), r['alleleID'], reportlib.close_accession_anchor()) + \
	'<td>' + alleleSymbol + '</td>\n' + \
	'<td>' + r['alleleName'] + '</td>\n' + \
	'<td>' + r['cellLine'] + '</td>\n' + \
	'<td>' + strain + '</td>\n'

    return s

    s = '%s%-30s%s ' % (reportlib.create_accession_anchor(r['markerID']), r['markerID'], reportlib.close_accession_anchor()) + \
	'%-30s ' % (r['markerSym']) + \
	'%-50s ' % (r['markerName']) + \
        '%s%-30s%s ' % (reportlib.create_accession_anchor(r['alleleID']), r['alleleID'], reportlib.close_accession_anchor()) + \
	'%-40s ' % (alleleSymbol) + \
	'%-50s ' % (r['alleleName']) + \
	'%-30s ' % (r['cellLine']) + \
	'%-50s ' % (strain)

    return s

def printHeader(fp, title):

    fp.write('</pre>\n')
    fp.write('<H2>%s</H2>' % (title))
    fp.write('<TABLE BORDER=3 WIDTH=100%>')
    fp.write('<th align = left valign=top>MGI Marker ID</th>')
    fp.write('<th align = left valign=top>Marker Symbol</th>')
    fp.write('<th align = left valign=top>Marker Name</th>')
    fp.write('<th align = left valign=to valign=topp>MGI Allele ID</th>')
    fp.write('<th align = left valign=top>Allele Symbol</th>')
    fp.write('<th align = left valign=top>Allele Name</th>')
    fp.write('<th align = left valign=top>ES Cell Line</th>')
    fp.write('<th align = left valign=top>ES Cell Line Strain</th>')
    fp.write('<th align = left valign=top>Genetic Background</th>')
    fp.write('<th align = left valign=top>IMSR Strain</th>')
    return

    fp.write('%-30s ' % ('MGI Marker ID'))
    fp.write('%-30s ' % ('Marker Symbol'))
    fp.write('%-50s ' % ('Marker Name'))
    fp.write('%-30s ' % ('MGI Allele ID'))
    fp.write('%-34s ' % ('Allele Symbol'))
    fp.write('%-50s ' % ('Allele Name'))
    fp.write('%-30s ' % ('ES Cell Line'))
    fp.write('%-50s ' % ('ES Cell Line Strain'))
    fp.write('%-50s ' % ('Genetic Background'))
    fp.write('%-30s ' % ('IMSR Strain'))
    fp.write(CRT*2)

#
# Main
#

fullreport = 'MGI_Knockout_Full'
publicreport = 'MGI_Knockout_Public'
notpublicreport = 'MGI_Knockout_NotPublic'

fp1 = reportlib.init(fullreport, printHeading = 0, outputdir = os.environ['REPORTOUTPUTDIR'], isHTML = 1)
fp2 = reportlib.init(publicreport, printHeading = 0, outputdir = os.environ['REPORTOUTPUTDIR'], isHTML = 1)
fp3 = reportlib.init(notpublicreport, printHeading = 0, outputdir = os.environ['REPORTOUTPUTDIR'], isHTML = 1)

printHeader(fp1, 'Full KnockOut Report')
printHeader(fp2, 'Public KnockOut Report')
printHeader(fp3, 'Non-Public KnockOut Report')

#
# select all targeted (knockout) alleles
#

db.sql('select a._Allele_key, alleleSym = a.symbol, alleleName = substring(a.name,1,50), alleleID = aa.accID, ' + \
	'markerSym = m.symbol, markerName = substring(m.name,1,50), markerID = ma.accID, ' + \
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
    value = regsub.gsub('<', '&lt;', r['label'])
    value = regsub.gsub('>', '&gt;', value)
    value = '%s%s%s' % (reportlib.create_imsrstrain_anchor(r['label']), value, reportlib.close_accession_anchor())
    if not imsr.has_key(key):
	imsr[key] = []
    imsr[key].append(value)

#
# process results
#
results = db.sql('select * from #knockouts order by markerSym', 'auto')

tbody = []

for r in results:

    printRecord = ''
    key = r['_Allele_key']
    imsrKey = r['alleleID']

    if imsr.has_key(imsrKey):
	printIMSR = '<td>%s</td>' % (string.join(imsr[imsrKey], ','))
    else:
	printIMSR = '<td>&nbsp;</td>'

    # record w/ genotype

    if geno.has_key(key):
	for g in geno[key]:
	    printRecord = printAllele() + '<td>%s</td>\n' % (g) + printIMSR + '</tr>\n'
            fp1.write(printRecord)

    # record w/out genotype

    if key in noGeno:
	printRecord = printAllele() + '<td>&nbsp;</td>\n' + printIMSR + '</tr>\n'
        fp1.write(printRecord)

    if imsr.has_key(imsrKey):
	fp2.write(printRecord)
    else:
	fp3.write(printRecord)

fp1.write('</TABLE>')
fp1.write('<pre>')

reportlib.finish_nonps(fp1, isHTML = 1)	# non-postscript file
reportlib.finish_nonps(fp2, isHTML = 1)	# non-postscript file
reportlib.finish_nonps(fp3, isHTML = 1)	# non-postscript file

