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
#		AND
#		at least one of the gene's knockout alleles is in IMSR
#
#	3.  Tab-delimited report where a row in the report represents:
#
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

allBLOG = '''
This report provides a list of all genes which have one or more 
published knockout or conditional knockout alleles. Each gene symbol is 
linked to its respective MGI Gene Detail page; each allele is linked to 
its MGI Phenotype and Allele Detail page; and a link is provided to the 
International Mouse Strain Resource (IMSR) strain if a repository holds 
mice carrying one or more of the listed knockout alleles.
<P>
Note that this list is a subset of phenotypic alleles for which there is 
data in MGI. It does not include gene-trapped alleles, targeted 
(Floxed/Frt) or targeted (knock-in) alleles, or those arising 
spontaneously or induced by chemicals or radiation.  To search the full 
spectrum of phenotypes or alleles in MGI, use the 
<a href="http://www.informatics.jax.org/searches/allele_form.shtml">Phenotype and Allele Query Form</a>.
To search repositories for specific strains carrying mutations of all 
types, use the <a href="http://www.informatics.jax.org/imsr/IMSRSearchForm.jsp">IMSR Search Form</a>.
<P>
'''

publicBLOG = '''
This report provides a list of genes which have one or more published 
knockout or conditional knockout alleles and at least one knockout 
allele is available through a public repository. Each gene symbol is 
linked to its respective MGI Gene Detail page; each allele is linked to 
its MGI Phenotype and Allele Detail page; and a link is provided to the 
International Mouse Strain Resource (IMSR) strains for mice carrying 
these knockout alleles.
<P> 
Note that this list is a subset of phenotypic alleles for which there is 
data in MGI and a subset of strains available from IMSR-participating 
repositories. It does not include gene-trapped alleles, targeted 
(Floxed/Frt) or targeted (knock-in) alleles, or those arising 
spontaneously or induced by chemicals or radiation.  To search the full 
spectrum of phenotypes or alleles in MGI, use the
<a href="http://www.informatics.jax.org/searches/allele_form.shtml">Phenotype and Allele Query Form</a>.
To search repositories for specific strains carrying mutations of all 
types, use the <a href="http://www.informatics.jax.org/imsr/IMSRSearchForm.jsp">IMSR Search Form</a>.
<P> 
'''

nonpublicBLOG = '''
Use this list to nominate genes that are known to have been 
knocked-out, but where mice carrying these knockouts are not available 
in public repositories.
<P>
This report provides a list of all genes which have one or more 
published knockout or conditional knockout alleles, but mice carrying 
these knockouts are not widely available. Each gene symbol is linked to 
its respective MGI Gene Detail page; each allele is linked to its MGI 
Phenotype and Allele Detail page.
<P>
Note that this list is a subset of phenotypic alleles for which there is 
data in MGI. It does not include gene-trapped alleles, targeted 
(Floxed/Frt) or targeted (knock-in) alleles, or those arising 
spontaneously or induced by chemicals or radiation.  To search the full 
spectrum of phenotypes or alleles in MGI, use the
<a href="http://www.informatics.jax.org/searches/allele_form.shtml">Phenotype and Allele Query Form</a>.
To search repositories for specific strains carrying mutations of all 
types, use the <a href="http://www.informatics.jax.org/imsr/IMSRSearchForm.jsp">IMSR Search Form</a>.
<P>
'''

def printMarkerHTML(r):

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

def printMarkerTAB(r):

    s = r['accID'] + TAB + \
        r['symbol'] + TAB + \
	r['name'] + TAB

    if repgen.has_key(r['_Marker_key']):
	s = s + repgen[r['_Marker_key']]
    s = s + TAB

    if reptran.has_key(r['_Marker_key']):
	s = s + reptran[r['_Marker_key']]
    s = s + TAB

    s = s + string.join(alleleTypes[r['_Marker_key']], '|') + TAB
    return s

def printAlleleHTML(a):

    symbol = regsub.gsub('<', 'beginss', a['symbol'])
    symbol = regsub.gsub('>', 'endss', symbol)
    symbol = regsub.gsub('beginss', '<sup>', symbol)
    symbol = regsub.gsub('endss', '</sup>', symbol)

    s = '%s%s%s<br>' % (reportlib.create_accession_anchor(a['accID']), symbol, reportlib.close_accession_anchor())

    return s

def printAlleleTAB(a):

    s = a['accID'] + '|'

    return s

def printHeaderHTML(fp, title, blog):

    fp.write('</pre>\n')
    fp.write('<H2>%s</H2>' % (title))
    fp.write(blog)
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

fp1 = reportlib.init(fullreport, printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'], isHTML = 1)
fp2 = reportlib.init(publicreport, printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'], isHTML = 1)
fp3 = reportlib.init(notpublicreport, printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'], isHTML = 1)

fp4 = reportlib.init(fullreport, printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])
fp5 = reportlib.init(publicreport, printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])
fp6 = reportlib.init(notpublicreport, printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])

printHeaderHTML(fp1, 'MGI All Knockouts Report.', allBLOG)
printHeaderHTML(fp2, 'Genes with Knockouts available through Public Repositories.', publicBLOG)
printHeaderHTML(fp3, 'Genes with Knockouts that are not yet available through Public Repositories.', nonpublicBLOG)

#
# select alleles
#

db.sql('select a._Marker_key, a._Allele_key, a.symbol, a.name, t.term, aa.accID ' + \
	'into #knockouts ' + \
	'from ALL_Allele a, ACC_Accession aa, VOC_Term t ' + \
	'where a._Allele_Status_key = 847114 ' + \
	'and a._Allele_Type_key in (847116, 847119, 847120) ' + \
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
imsrHTML = {}
imsrTAB = {}
for r in results:

    key = r['_Marker_key']
    value = regsub.gsub('<', 'beginss', r['label'])
    value = regsub.gsub('>', 'endss', value)
    value = regsub.gsub('beginss', '<sup>', value)
    value = regsub.gsub('endss', '</sup>', value)
    value = '%s%s%s' % (reportlib.create_imsrstrain_anchor(r['label']), value, reportlib.close_accession_anchor())

    if not imsrHTML.has_key(key):
	imsrHTML[key] = []
    if value not in imsrHTML[key]:
        imsrHTML[key].append(value)

    value = r['label']
    if not imsrTAB.has_key(key):
	imsrTAB[key] = []
    if value not in imsrTAB[key]:
        imsrTAB[key].append(value)

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

    mstrHTML = ''
    astrHTML = ''
    mstrTAB = ''
    astrTAB = ''

    isLexicon = 0
    isDeltagen = 0

    key = r['_Marker_key']
    mstrHTML = printMarkerHTML(r)
    mstrTAB = printMarkerTAB(r)

    # for each allele....

    for a in alleles[key]:

	if string.find(a['name'], 'Lexicon') >= 0:
	    isLexicon = 1

	if string.find(a['name'], 'lexicon') >= 0:
	    isLexicon = 1

	if string.find(a['name'], 'Deltagen') >= 0:
	    isDeltagen = 1

	astrHTML = astrHTML + printAlleleHTML(a)
	astrTAB = astrTAB + printAlleleTAB(a)

    if imsrHTML.has_key(key):
	istrHTML = '<td>%s</td>' % (string.join(imsrHTML[key], '<br>'))
	istrTAB = string.join(imsrTAB[key], '|') + CRT
    else:
	istrHTML = '<td>&nbsp;</td>'
	istrTAB = CRT

    fp1.write('%s<td>%s</td>%s' % (mstrHTML, astrHTML, istrHTML))
    fp4.write('%s%s\t%s' % (mstrTAB, astrTAB, istrTAB))

    if imsrHTML.has_key(key):
        fp2.write('%s<td>%s</td>%s' % (mstrHTML, astrHTML, istrHTML))
        fp5.write('%s%s\t%s' % (mstrTAB, astrTAB, istrTAB))
    else:
	if not isLexicon and not isDeltagen:
          fp3.write('%s<td>%s</td>%s' % (mstrHTML, astrHTML, istrHTML))
          fp6.write('%s%s%s' % (mstrTAB, astrTAB, istrTAB))

fp1.write('</TABLE>')
fp1.write('<pre>')
fp2.write('</TABLE>')
fp2.write('<pre>')
fp3.write('</TABLE>')
fp3.write('<pre>')

reportlib.finish_nonps(fp1, isHTML = 1)	# non-postscript file
reportlib.finish_nonps(fp2, isHTML = 1)	# non-postscript file
reportlib.finish_nonps(fp3, isHTML = 1)	# non-postscript file

reportlib.finish_nonps(fp4)	# non-postscript file
reportlib.finish_nonps(fp5)	# non-postscript file
reportlib.finish_nonps(fp6)	# non-postscript file

