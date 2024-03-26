
'''
#
# MGI_Recombinase_Full
#
# Report:
#       Cre Alleles
#
#	Tab-delimited report of all ALL_Cre_Cache records
#
#	field 1: Driver note
#	field 2: Allele Symbol
#	field 3: Name
#	field 4: Detected In 
#		pipe-deimited list of cre system with expressed assay results (expressed = true)
#	field 5: Absent in
#		pipe-deimited list of cre system with expressed assay results (expressed = false)
#		if system already exists in field 6, skip it
#	field 6: IMSR Strain
#		pipe-delimited list of IMSR lines available (cell lines and non-cell lines)
#	field 7: Allele ID
#
# This generates both an HTML version and a tab-delimited version of the report.
# It also superscripts the allele symbol, imsr label
#
# Usage:
#       MGI_Recombinase_Full.py
#
# History:
#
# lec	03/24/2016
#	- TR12223/gxd anatomy II/replace AD with cre system
#
# kstone 12/18/2014
# 	- removed IMSR db query, replaced with Solr csv file
#
# lec   10/24/2014
#       - TR11750/postgres complient
#
# lec   03/28/2012
#       - TR 11027; WI_URL changes for MGI 5.0
#
# lec	01/19/2010
#	- TR 10042; added JAX print header
#
# lec	09/29/2009
#	- created
#
'''
 
import csv
import sys 
import os
import reportlib
import urllib.request, urllib.parse, urllib.error
import db

db.setTrace()

CRT = reportlib.CRT
TAB = reportlib.TAB

WI_URL = os.environ['WI_URL']

IMSR_CSV = os.environ['IMSR_STRAINS_CSV']

# functions
def createImsrAlleleToStrainDict():
        """
        Reads in the IMSR allStrains.csv file and generates
        a map of allele ID to strain name
        """
        strainNameMap = {}
        csvfile = open(IMSR_CSV, 'r')
        reader = csv.reader(csvfile)
        for row in reader:
                allele_ids = row[0]
                strain_name = row[2]
                if allele_ids and strain_name:
                        for id in allele_ids.split(','):
                                strainNameMap.setdefault(id, []).append(strain_name)

        return strainNameMap

introBLOG = '''
<p>
This report provides a list of all recombinase-containing alleles in the MGI database.
<p>
Each allele symbol is linked to its respective MGI Allele Detail page, containing phenotypic
and disease model data; each Cre System for an allele is linked to its MGI Recombinase Detail page.
A link is provided to the International Mouse Strain Resource (IMSR) strain if a repository holds mice
carrying the listed allele.
<p>
To search for floxed, frt, or other recombinase target-containing alleles in MGI, use the
<A HREF="%ssearches/allele_form.shtml">Phenotype and Alleles Query Form</A>.
<p>
To search repositories for specific strains carrying mutations of all types, use the 
<A HREF="http://www.findmice.org//IMSRSearchForm.jsp">IMSR Search Form</A>.
<p>
''' % (WI_URL)

BEGTD = '<td style="white-space:nowrap;"><font size="-1">'
BEGTDWRAP = '<td><font size="-1">'
ENDTD = '</font></td>'
BREAK = '<br>'
BLANKFIELD = '%s&nbsp;%s' % (BEGTD, ENDTD)

CRE_ANCHOR = '<A HREF="%srecombinase/specificity?id=%s&system=%s">'
CLOSE_ANCHOR = '</A>'

def printHeaderHTML():
    #
    # write header to HTML file
    #

    fpHTML.write('</pre>\n')
    fpHTML.write('<H2>MGI Recombinase Alleles Report</H2>')
    fpHTML.write(introBLOG)
    fpHTML.write('<TABLE BORDER=3 WIDTH=100%>')
    fpHTML.write('''
    <th align=left valign=top>Driver</th>
    <th align=left valign=top>Allele Symbol</th>
    <th align=left valign=top>Name</th>
    <th align=left valign=top>Detected in</th>
    <th align=left valign=top>Absent in</th>
    <th align=left valign=top>IMSR Strain</th>
    <th align=left valign=top>Allele ID</th>
    ''')

def printHeaderTAB():
    #
    # write record to tab-delimited file
    #

    fpTAB.write('# MGI Recombinase Alleles Report\n')
    fpTAB.write('# field 1: Driver\n')
    fpTAB.write('# field 2: Allele Symbol\n')
    fpTAB.write('# field 3: Name\n')
    fpTAB.write('# field 4: Detected in (Cre systems with expressed assay results)\n')
    fpTAB.write('# field 5: Absent in (Cre systems with notexpressed assay results)\n')
    fpTAB.write('# field 6: IMSR Strain (list of IMSR lines available (cell lines, cryo materials, and live)\n')
    fpTAB.write('# field 7: Allele ID\n')
    fpTAB.write(CRT*2)

def writeHTML(r, imsrHTML):
    #
    # write record to HTML file
    #

    key = r['_Allele_key']
    accID = r['accID']
    driverNote = str.replace(r['driverNote'], '\n', '')

    # superscript the symbol

    symbol = str.replace(r['symbol'], '<', 'beginss')
    symbol = str.replace(symbol, '>', 'endss')
    symbol = str.replace(symbol, 'beginss', '<sup>')
    symbol = str.replace(symbol, 'endss', '</sup>')

    if r['name'] == r['markerName']:
        name = r['name']
    else:
        name = r['markerName'] + '; ' + r['name']
 
    s = '<tr>' + \
        BEGTD + driverNote + ENDTD + \
        BEGTD + reportlib.create_accession_anchor(r['accID']) + symbol + CLOSE_ANCHOR + ENDTD + \
        BEGTDWRAP + name + ENDTD

    if key in expressedHTML:
        s = s + BEGTD + '%s' % (BREAK.join(expressedHTML[key])) + ENDTD
    else:
        s = s + BLANKFIELD
      
    if key in notexpressedHTML:
        s = s + BEGTD + '%s' % (BREAK.join(notexpressedHTML[key])) + ENDTD
    else:
        s = s + BLANKFIELD
      
    if accID in imsrHTML:
        s = s + BEGTD + '%s' % (BREAK.join(imsrHTML[accID])) + ENDTD
    else:
        s = s + BLANKFIELD
      
    s = s + BEGTD + r['accID'] + ENDTD + CRT

    fpHTML.write(s)

def writeTAB(r, imsrTAB):
    #
    # write record to tab-delimited file
    #

    key = r['_Allele_key']
    accID = r['accID']
    driverNote = str.replace(r['driverNote'], '\n', '')

    if r['name'] == r['markerName']:
        name = r['name']
    else:
        name = r['markerName'] + '; ' + r['name']
 
    fpTAB.write(driverNote + TAB + \
                r['symbol'] + TAB + \
                name + TAB)

    if key in expressedTAB:
      fpTAB.write('|'.join(expressedTAB[key]))
    fpTAB.write(TAB)

    if key in notexpressedTAB:
      fpTAB.write('|'.join(notexpressedTAB[key]))
    fpTAB.write(TAB)

    if accID in imsrTAB:
        fpTAB.write('|'.join(imsrTAB[accID]))
    fpTAB.write(TAB)

    fpTAB.write(r['accID'] + CRT)

#
# Main
#

fpHTML = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None, isHTML = 1)
fpTAB = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 'JAX')

printHeaderHTML()
printHeaderTAB()

#
# select all Cre Alleles
#
# the Marker name was added to this report at the last minue,
# and my preferred method would have been to add the marker key and name
# to the Cre cache table...but at this stage in the process, I'm not going to do that...
# which is why there needs to be another join to the allele/marker tables to grab
# the marker name
#

db.sql('''
       select distinct c._Allele_key, c.symbol, c.name, c.driverNote, 
                c.accID, rtrim(m.name) as markerName
       into temporary table cre
       from ALL_Cre_Cache c, ALL_Allele aa, MRK_Marker m
       where c._Allele_key = aa._Allele_key
       and aa._Marker_key = m._Marker_key
       ''', None)

db.sql('create index cre_idx1 on cre(_Allele_key)', None)

# cre systems of expressed structures

db.sql('''
      select distinct cc._Allele_key, cc.cresystemlabel, c.accID
      into temporary table expressed
      from cre c, ALL_Cre_Cache cc
      where c._Allele_key = cc._Allele_key
      and cc.expressed = 1
      ''', None)
db.sql('create index expressed_idx1 on expressed(_Allele_key)', None)

expressedTAB = {}
expressedHTML = {}
results = db.sql('select * from expressed', 'auto')
for r in results:
    key = r['_Allele_key']
    value = r['cresystemlabel']

    if key not in expressedTAB:
        expressedTAB[key] = []
    expressedTAB[key].append(value)

    value = CRE_ANCHOR % (WI_URL, r['accID'], urllib.parse.quote_plus(r['cresystemlabel'])) + r['cresystemlabel'] + CLOSE_ANCHOR
    if key not in expressedHTML:
        expressedHTML[key] = []
    expressedHTML[key].append(value)

# cre systems of not-expressed structures
# where the system does not exist in the expressed group

results = db.sql('''
      select distinct cc._Allele_key, cc.cresystemlabel, cc.accID
      from cre c, ALL_Cre_Cache cc
      where c._Allele_key = cc._Allele_key
      and cc.expressed = 0
      and not exists (select 1 from expressed e
                        where cc._Allele_key = e._Allele_key
                        and cc.cresystemlabel = e.cresystemlabel)
      ''', 'auto')
notexpressedTAB = {}
notexpressedHTML = {}
for r in results:
    key = r['_Allele_key']
    value = r['cresystemlabel']

    if key not in notexpressedTAB:
        notexpressedTAB[key] = []
    notexpressedTAB[key].append(value)

    value = CRE_ANCHOR % (WI_URL, r['accID'], urllib.parse.quote_plus(r['cresystemlabel'])) + r['cresystemlabel'] + CLOSE_ANCHOR
    if key not in notexpressedHTML:
        notexpressedHTML[key] = []
    notexpressedHTML[key].append(value)

#
# select those alleles that are in IMSR
#

results = db.sql('''
        select distinct c.accID
        from cre c
        ''', 'auto')

# map of imsr allele ID to strain names
imsrStrainMap = createImsrAlleleToStrainDict()

imsrTAB = {}
imsrHTML = {}
for r in results:
    accID = r['accID']
    if accID not in imsrStrainMap:
        continue

    for strain in imsrStrainMap[accID]:

        imsrTAB.setdefault(accID, []).append(strain)

        value = str.replace(strain, '<', 'beginss')
        value = str.replace(value, '>', 'endss')
        value = str.replace(value, 'beginss', '<sup>')
        value = str.replace(value, 'endss', '</sup>')
        value = '%s%s%s' % (reportlib.create_imsrstrain_anchor(strain), value, reportlib.close_accession_anchor())

        imsrHTML.setdefault(accID, []).append(value)

#
# process results
#
results = db.sql('select * from cre order by driverNote', 'auto')

for r in results:
    writeHTML(r, imsrHTML)
    writeTAB(r, imsrTAB)

#
# clean up
#

fpHTML.write('</TABLE>')
fpHTML.write('<pre>')

reportlib.finish_nonps(fpHTML, isHTML = 1)  # non-postscript file
reportlib.finish_nonps(fpTAB)	# non-postscript file

