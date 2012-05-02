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
#	field 1: Driver note
#	field 2: Allele Symbol
#	field 3: Name
#	field 4: Detected In 
#		pipe-deimited list of anatomical system with expressed assay results (expressed = true)
#	field 5: Absent in
#		pipe-deimited list of anatomical system with expressed assay results (expressed = false)
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
 
import sys 
import os
import string
import reportlib

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db


CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

WI_URL = os.environ['WI_URL']
IMSR = os.environ['IMSR_DBNAME']

introBLOG = '''
<p>
This report provides a list of all recombinase-containing alleles in the MGI database.
<p>
Each allele symbol is linked to its respective MGI Allele Detail page, containing phenotypic
and disease model data; each Anatomical System for an allele is linked to its MGI Recombinase Detail page.
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

ALLELE_ANCHOR = '<A HREF="%ssearches/accession_report.cgi?id=%s" target="_blank">'
CRE_ANCHOR = '<A HREF="%srecombinase/specificity?id=%s&systemKey=%s">'
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
    fpTAB.write('# field 4: Detected in (anatomical systems with expressed assay results)\n')
    fpTAB.write('# field 5: Absent in (anatomical systems with notexpressed assay results)\n')
    fpTAB.write('# field 6: IMSR Strain (list of IMSR lines available (cell lines, cryo materials, and live)\n')
    fpTAB.write('# field 7: Allele ID\n')
    fpTAB.write(CRT*2)

def writeHTML(r):
    #
    # write record to HTML file
    #

    key = r['_Allele_key']
    driverNote = string.replace(r['driverNote'], '\n', '')

    # superscript the symbol

    symbol = string.replace(r['symbol'], '<', 'beginss')
    symbol = string.replace(symbol, '>', 'endss')
    symbol = string.replace(symbol, 'beginss', '<sup>')
    symbol = string.replace(symbol, 'endss', '</sup>')

    if r['name'] == r['markerName']:
	name = r['name']
    else:
	name = r['markerName'] + '; ' + r['name']
 
    s = '<tr>' + \
        BEGTD + driverNote + ENDTD + \
        BEGTD + ALLELE_ANCHOR % (WI_URL, r['accID']) + symbol + CLOSE_ANCHOR + ENDTD + \
        BEGTDWRAP + name + ENDTD

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
      
    s = s + BEGTD + r['accID'] + ENDTD + CRT

    fpHTML.write(s)

def writeTAB(r):
    #
    # write record to tab-delimited file
    #

    key = r['_Allele_key']
    driverNote = string.replace(r['driverNote'], '\n', '')

    if r['name'] == r['markerName']:
	name = r['name']
    else:
	name = r['markerName'] + '; ' + r['name']
 
    fpTAB.write(driverNote + TAB + \
                r['symbol'] + TAB + \
		name + TAB)

    if expressedTAB.has_key(key):
      fpTAB.write(string.join(expressedTAB[key], '|'))
    fpTAB.write(TAB)

    if notexpressedTAB.has_key(key):
      fpTAB.write(string.join(notexpressedTAB[key], '|'))
    fpTAB.write(TAB)

    if imsrTAB.has_key(key):
        fpTAB.write(string.join(imsrTAB[key], '|'))
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
		a.accID, rtrim(m.name) as markerName
       into #cre
       from ALL_Cre_Cache c, ACC_Accession a, ALL_Allele aa, MRK_Marker m
       where c._Allele_key = a._Object_key
       and a._MGIType_key = 11
       and a._LogicalDB_key = 1 
       and a.prefixPart = 'MGI:' 
       and a.preferred = 1
       and c._Allele_key = aa._Allele_key
       and aa._Marker_key = m._Marker_key
       ''', None)

db.sql('create index cre_idx1 on #cre(_Allele_key)', None)

# anatomical systems of expressed structures

db.sql('''
      select distinct cc._Allele_key, cc._System_key, cc.system, c.accID
      into #expressed
      from #cre c, ALL_Cre_Cache cc
      where c._Allele_key = cc._Allele_key
      and cc.expressed = 1
      ''', None)
db.sql('create index expressed_idx1 on #expressed(_Allele_key)', None)

expressedTAB = {}
expressedHTML = {}
results = db.sql('select * from #expressed', 'auto')
for r in results:
    key = r['_Allele_key']
    value = r['system']

    if not expressedTAB.has_key(key):
        expressedTAB[key] = []
    expressedTAB[key].append(value)

    #value = CRE_ANCHOR % (WI_URL, key, r['_System_key']) + r['system'] + CLOSE_ANCHOR
    value = CRE_ANCHOR % (WI_URL, r['accID'], r['_System_key']) + r['system'] + CLOSE_ANCHOR
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

    value = string.replace(r['label'], '<', 'beginss')
    value = string.replace(value, '>', 'endss')
    value = string.replace(value, 'beginss', '<sup>')
    value = string.replace(value, 'endss', '</sup>')
    value = '%s%s%s' % (reportlib.create_imsrstrain_anchor(r['label']), value, reportlib.close_accession_anchor())

    if not imsrHTML.has_key(key):
        imsrHTML[key] = []
    if value not in imsrHTML[key]:
        imsrHTML[key].append(value)

#
# process results
#
results = db.sql('select * from #cre order by driverNote', 'auto')

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

