#!/usr/local/bin/python

'''
#
# MGI_Recombinase_Full
#
# Report:
#       Cre Alleles
#
#	Tab-delimited report of all ALL_Cre_Cache records
#
#	field 1: Allele ID
#	field 2: Allele Symbol
#	field 3: Allele Name
#	field 4: Allele Type
#	field 5: Driver note
#	field 6: Detected In 
#		pipe-deimited list of anatomical system with expressed assay results (expressed = true)
#	field 7: Absent in
#		pipe-deimited list of anatomical system with expressed assay results (expressed = false)
#		if system already exists in field 6, skip it
#	field 8: IMSR Strain
#		pipe-delimited list of IMSR lines available (cell lines and non-cell lines)
#
# This generates both an HTML version and a tab-delimited version of the report.
# It also superscripts the allele symbol, imsr label
#
# Usage:
#       MGI_Recombinase_Full.py
#
# History:
#
# lec	09/29/2009
#	- created
#
'''
 
import sys 
import os
import db
import re
import string
import reportlib
import regsub

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

WI_URL = os.environ['WI_URL']
IMSR = os.environ['IMSR_DBNAME']

introBLOG = '''
<p>
This report provides a list of all recombinases-containing alleles in the MGI database.
<p>
Each allele symbol is linked to its respective MGI Allele Detail page, containing phenotypic
and disease model data; each Anatomical System for an allele is linked to its MGI Recombinase Detail page.
A link is provided to the International Mouse Strain Resource (IMSR) strain if a repository holds mice
carrying the listed allele.
<p>
To search for floxed, frt, or other recombinase target-containing alleles in MGI, use the
Phenotype and <A HREF="%ssearches/allele_form.shtml">Alleles Query Form</A>.
<p>
To search repositories for specific strains carrying mutations of all types, use the 
<A HREF="http://www.findmice.org//IMSRSearchForm.jsp">IMSR Search Form</A>.
<p>
''' % (WI_URL)

BEGTD = '<td style="white-space:nowrap;">'
ENDTD = '</td>'
BREAK = '<br>'
BLANKFIELD = '%s&nbsp;%s' % (BEGTD, ENDTD)

ALLELE_ANCHOR = '<A HREF="%ssearches/accession_report.cgi?id=%s" target="_blank">'
CRE_ANCHOR = '<A HREF="%sjavawi2/servlet/WIFetch?page=creSpecificity&alleleKey=%s&systemKey=%s" target="_blank">'
CLOSE_ANCHOR = '</A>'

def printHeaderHTML():
    #
    # write header to HTML file
    #

    fpHTML.write('</pre>\n')
    fpHTML.write('<H2>MGI Recombinase Alleles</H2>')
    fpHTML.write(introBLOG)
    fpHTML.write('<TABLE BORDER=3 WIDTH=100%>')
    fpHTML.write('''
    <th align=left valign=top>Allele ID</th>
    <th align=left valign=top>Allele Symbol</th>
    <th align=left valign=top>Allele Name</th>
    <th align=left valign=top>Allele Type</th>
    <th align=left valign=top>Driver</th>
    <th align=left valign=top>Detected in</th>
    <th align=left valign=top>Absent in</th>
    <th align=left valign=top>IMSR Strain</th>
    ''')

def printHeaderTAB():
    #
    # write record to tab-delimited file
    #

    fpTAB.write('# MGI - Recombinase Alleles\n')
    fpTAB.write('# field 1: Allele ID\n')
    fpTAB.write('# field 2: Allele Symbol\n')
    fpTAB.write('# field 3: Allele Name\n')
    fpTAB.write('# field 4: Allele Type\n')
    fpTAB.write('# field 5: Driver\n')
    fpTAB.write('# field 6: Detected in (anatomical systems with expressed assay results)\n')
    fpTAB.write('# field 7: Absent in (anatomical systems with notexpressed assay results)\n')
    fpTAB.write('# field 8: IMSR Strain (list of IMSR lines available (cell lines, cryo materials, and live)\n')
    fpTAB.write(CRT*2)

def writeHTML(r):
    #
    # write record to HTML file
    #

    key = r['_Allele_key']

    # superscript the symbol

    symbol = regsub.gsub('<', 'beginss', r['symbol'])
    symbol = regsub.gsub('>', 'endss', symbol)
    symbol = regsub.gsub('beginss', '<sup>', symbol)
    symbol = regsub.gsub('endss', '</sup>', symbol)

    s = '<tr>' + \
        BEGTD + r['accID'] + ENDTD + '\n' + \
        BEGTD + ALLELE_ANCHOR % (WI_URL, r['accID']) + symbol + CLOSE_ANCHOR + ENDTD + \
        BEGTD + r['name'] + ENDTD + \
        BEGTD + r['alleleType'] + ENDTD + \
        BEGTD + r['driverNote'] + ENDTD

    if expressedHTML.has_key(key):
	s = s + BEGTD + '%s' % (string.join(expressedHTML[key], BREAK)) + ENDTD
    else:
        s = s + BLANKFIELD
      
    if notexpressedHTML.has_key(key):
	s = s + BEGTD + '%s' % (string.join(notexpressedHTML[key], BREAK)) + ENDTD
    else:
        s = s + BLANKFIELD
      
    if imsrHTML.has_key(key):
	s = s + BEGTD + '%s' % (string.join(imsrHTML[key], BREAK)) + ENDTD
    else:
        s = s + BLANKFIELD
      
    s = s + CRT

    fpHTML.write(s)

def writeTAB(r):
    #
    # write record to tab-delimited file
    #

    key = r['_Allele_key']

    fpTAB.write(r['accID'] + TAB + \
                 r['symbol'] + TAB + \
		 r['name'] + TAB + \
		 r['alleleType'] + TAB + \
		 r['driverNote'] + TAB)

    if expressedTAB.has_key(key):
      fpTAB.write(string.join(expressedTAB[key], '|'))
    fpTAB.write(TAB)

    if notexpressedTAB.has_key(key):
      fpTAB.write(string.join(notexpressedTAB[key], '|'))
    fpTAB.write(TAB)

    if imsrTAB.has_key(key):
        fpTAB.write(string.join(imsrTAB[key], '|'))
    fpTAB.write(CRT)

#
# Main
#

fpHTML = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None, isHTML = 1)
fpTAB = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

printHeaderHTML()
printHeaderTAB()

#
# select all Cre Alleles
#

db.sql('''
       select distinct c._Allele_key, c.symbol, c.name, c.alleleType, c.driverNote, a.accID
       into #cre
       from ALL_Cre_Cache c, ACC_Accession a
       where c._Allele_key = a._Object_key
       and a._MGIType_key = 11
       and a._LogicalDB_key = 1 
       and a.prefixPart = "MGI:" 
       and a.preferred = 1
       ''', None)

db.sql('create index idx1 on #cre(_Allele_key)', None)

# anatomical systems of expressed structures

db.sql('''
      select distinct cc._Allele_key, cc._System_key, cc.system
      into #expressed
      from #cre c, ALL_Cre_Cache cc
      where c._Allele_key = cc._Allele_key
      and cc.expressed = 1
      ''', None)
db.sql('create index idx1 on #expressed(_Allele_key)', None)

expressedTAB = {}
expressedHTML = {}
results = db.sql('select * from #expressed', 'auto')
for r in results:
    key = r['_Allele_key']
    value = r['system']

    if not expressedTAB.has_key(key):
        expressedTAB[key] = []
    expressedTAB[key].append(value)

    value = CRE_ANCHOR % (WI_URL, key, r['_System_key']) + r['system'] + CLOSE_ANCHOR
    if not expressedHTML.has_key(key):
        expressedHTML[key] = []
    expressedHTML[key].append(value)

# anatomical systems of not-expressed structures
# where the system does not exist in the expressed group

results = db.sql('''
      select distinct cc._Allele_key, cc._System_key, cc.system
      from #cre c, ALL_Cre_Cache cc
      where c._Allele_key = cc._Allele_key
      and cc.expressed = 0
      and not exists (select 1 from #expressed e
			where cc._Allele_key = e._Allele_key
			and cc.system = e.system)
      ''', 'auto')
notexpressedTAB = {}
notexpressedHTML = {}
for r in results:
    key = r['_Allele_key']
    value = r['system']

    if not notexpressedTAB.has_key(key):
        notexpressedTAB[key] = []
    notexpressedTAB[key].append(value)

    value = CRE_ANCHOR % (WI_URL, key, r['_System_key']) + r['system'] + CLOSE_ANCHOR
    if not notexpressedHTML.has_key(key):
        notexpressedHTML[key] = []
    notexpressedHTML[key].append(value)

#
# select those alleles that are in IMSR
#

results = db.sql('''
        select distinct c._Allele_key, c.accID, ls.label 
        from #cre c, %s..Accession ac, %s..Label ls, %s..SGAAssoc sga 
        where c.accID = ac.accID 
        and ac._Object_key = sga._Allele_key
        and ac._IMSRType_key = 3 
        and sga._Strain_key = ls._Object_key 
        and ls._IMSRType_key = 1 
        and ls.labelType = "N" 
        ''' % (IMSR, IMSR, IMSR), 'auto')
imsrTAB = {}
imsrHTML = {}
for r in results:
    key = r['_Allele_key']
    value = r['label']

    if not imsrTAB.has_key(key):
        imsrTAB[key] = []
    if value not in imsrTAB[key]:
        imsrTAB[key].append(value)

    value = regsub.gsub('<', 'beginss', r['label'])
    value = regsub.gsub('>', 'endss', value)
    value = regsub.gsub('beginss', '<sup>', value)
    value = regsub.gsub('endss', '</sup>', value)
    value = '%s%s%s' % (reportlib.create_imsrstrain_anchor(r['label']), value, reportlib.close_accession_anchor())

    if not imsrHTML.has_key(key):
        imsrHTML[key] = []
    if value not in imsrHTML[key]:
        imsrHTML[key].append(value)

#
# process results
#
results = db.sql('select * from #cre order by symbol', 'auto')

for r in results:
    writeHTML(r)
    writeTAB(r)

#
# clean up
#

fpHTML.write('</TABLE>')
fpHTML.write('<pre>')

reportlib.finish_nonps(fpHTML, isHTML = 1)  # non-postscript file
reportlib.finish_nonps(fpTAB)	# non-postscript file

