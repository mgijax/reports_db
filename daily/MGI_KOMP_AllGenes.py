#!/usr/local/bin/python

'''
#
# MGI_KOMP_AllGenes.py
#
# Reports:
#
#       TR 8028
#
#	MGI Mouse Markers:
#	. of type "Gene"
#	. associated with at least one:  
#	  NCBI Gene Model
#	  Ensembl Gene Model
#	  VEGA Gene Model
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
# lec	12/08/2006
#	- created, TR 8028
#
'''
 
import sys 
import os
import db
import string
import reportlib

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

introBLOG = '''
This report provides a list of all genes that may be considered potential candidates for being knocked out through the NIH KOMP project. Genes are sorted by chromosome, and alphabetically by gene symbol within a chromosome.
<P>

<UL>
<LI>Each gene symbol is linked to its respective MGI Gene Detail page for additional information.
<P>

<LI>Genes with "KOMP" in the KOMP Target column are in the current list of 4,400+ genes in the KOMP target gene list. The KOMP target gene list consists of those genes in the mouse cCDS uniquely identified as having a consistent one-to-one relationship among genes assigned by MGI, EntrezGene, ENSEMBL, and VEGA and have never been reported as knocked out or gene trapped in the public domain (see 4th bullet point, below).
<P>

<LI>Genes with "Regeneron" or "CSD" in their respective columns indicate genes that currently are in the KOMP pipeline at the Regeneron or CHORI/Sanger/UCDavis KOMP project sites.
<P>

<LI>Genes with "+" in the Other Knockout Source column are those genes that have been reported as knocked out or gene trapped in the public domain (i.e., these genes are listed by the IGTC as having one or more gene traps or are recorded in MGI as having one or more knockouts or gene traps; or are in the pipeline of EUCOMM or the Sanger Institute towards being knocked out).
<P>
</UL>

This All Genes List was created as follows. All genetic markers assigned type GENE in MGI were screened for those that also have at least one of the following types of IDs associated with them: NCBI EntrezGene or ENSEMBL or Vega. Only nuclear encoded genes are included (no mitochondrial encoded genes) and each gene has been placed on the genome (i.e., has a chromosomal location and genome coordinates in NCBI Build 36). 
<P>
'''

BEGTD = '<td><font size="-1">'
ENDTD = '</font></td>'
BLANKFIELD = '%s&nbsp;%s' % (BEGTD, ENDTD)

fpHTML = None
fpTAB = None

kompGeneListFile = os.environ['KOMP_GENELIST']
kompTargetListFile = os.environ['KOMP_TARGETLIST']

KOMPHTMLTAG = BEGTD + 'KOMP' + ENDTD
KOMPTAG = 'KOMP'

kompGeneDict = {}		# dictionary; key = MGI ID, values = dictionary of komp fields
kompTargetList = []		# a list of MGI IDs that appear in the KOMP Target list

ncbimodel = {}		# dictionary of MGI marker key:ncbi gene model ids
ensemblmodel = {}	# dictionary of MGI marker key:ensembl gene model ids
vegamodel = {}		# dictionary of MGI marker key:vega gene model ids

def init():

    global fpHTML, fpTAB, kompGeneDict, kompTargetList

    fpHTML = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'], isHTML = 1)
    fpTAB = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])

    printHeaderHTML()
    printHeaderTAB()

    # KOMP gene list
    # key = MGI ID
    # values = dictionary of key:value pairs

    kompfile = open(kompGeneListFile, 'r')
    for line in kompfile.readlines():
	tokens = string.split(line[:-1], TAB)

	mgiID = tokens[1]

        k = {}
	k['cCDS'] = tokens[4]
	k['regeneron'] = ''
	k['csd'] = ''
	k['eucomm'] = ''
	k['sanger'] = ''
	k['igtc'] = ''
	k['mgi'] = ''

	if len(tokens) >= 8:
	    k['regeneron'] = tokens[7]
	if len(tokens) >= 9:
	    k['csd'] = tokens[8]
	if len(tokens) >= 10:
	    k['eucomm'] = tokens[9]
	if len(tokens) >= 11:
	    k['sanger'] = tokens[10]
	if len(tokens) >= 12:
	    k['igtc'] = tokens[11]
	if len(tokens) >= 13:
	    k['mgi'] = tokens[12]

	kompGeneDict[mgiID] = k

    kompfile.close()

    # KOMP target list

    kompfile = open(kompTargetListFile, 'r')
    for line in kompfile.readlines():
	tokens = string.split(line[:-1], TAB)
	mgiID = tokens[1]
	kompTargetList.append(mgiID)
    kompfile.close()

def writeHTML(r):

    # print marker records in HTML format

    key = r['_Marker_key']

    s = '<tr>' + \
	BEGTD + r['accID'] + '</td>\n' + \
        '%s%s%s%s%s\n' % (BEGTD, reportlib.create_accession_anchor(r['accID']), r['symbol'], reportlib.close_accession_anchor(), ENDTD) + \
	'<td width="25%"><font size="-1">' + r['name'] + ENDTD

    # if MGI ID found in KOMP target file...

    if r['accID'] in kompTargetList:
	s = s + KOMPHTMLTAG
    else:
	s = s + BLANKFIELD
	
    cCDS = ''
    if kompGeneDict.has_key(r['accID']):

	k = kompGeneDict[r['accID']]

	cCDS = k['cCDS']

	if len(k['regeneron']) > 0:
	    s = s + BEGTD + k['regeneron'] + ENDTD
	else:
	    s = s + BLANKFIELD
	
	if len(k['csd']) > 0:
	    s = s + BEGTD + k['csd'] + ENDTD
	else:
	    s = s + BLANKFIELD
	
	if len(k['eucomm']) > 0 or len(k['sanger']) > 0 or len(k['igtc']) > 0 or len(k['mgi']) > 0:
	    s = s + BEGTD + '+' + ENDTD
	else:
	    s = s + BLANKFIELD

    else:
	s = s + BLANKFIELD + BLANKFIELD + BLANKFIELD

    s = s + BEGTD + r['chromosome'] + ENDTD + \
	BEGTD + str(r['startC']) + ENDTD + \
	BEGTD + str(r['endC']) + ENDTD + \
	BEGTD + str(r['strand']) + ENDTD

    # if NCBI gene model exists

    if ncbimodel.has_key(key):
	s = s + '%s%s%s' % (BEGTD, string.join(ncbimodel[key], "<BR>"), ENDTD)
    else:
	s = s + BLANKFIELD

    # if Ensembl gene model exists

    if ensemblmodel.has_key(key):
	s = s + '%s%s%s' % (BEGTD, string.join(ensemblmodel[key], "<BR>"), ENDTD)
    else:
	s = s + BLANKFIELD

    # if VEGA gene model exists

    if vegamodel.has_key(key):
	s = s + '%s%s%s' % (BEGTD, string.join(vegamodel[key], "<BR>"), ENDTD)
    else:
	s = s + BLANKFIELD

    if len(cCDS) > 0:
	s = s + BEGTD + cCDS + ENDTD
    else:
	s = s + BLANKFIELD

    s = s + CRT

    return s

def writeTAB(r):

    # print marker records in tab-delimited format

    key = r['_Marker_key']

    s = r['accID'] + TAB + r['symbol'] + TAB + r['name'] + TAB

    # if MGI ID found in KOMP Target file...

    if r['accID'] in kompTargetList:
	s = s + KOMPTAG
    s = s + TAB
	
    # if MGI ID found in KOMP file...

    cCDS = ''
    if kompGeneDict.has_key(r['accID']):

	k = kompGeneDict[r['accID']]

	cCDS = k['cCDS']
	s = s + k['regeneron'] + TAB
	s = s + k['csd'] + TAB

	if len(k['eucomm']) > 0 or len(k['sanger']) > 0 or len(k['igtc']) > 0 or len(k['mgi']) > 0:
	    s = s + '+'
	s = s + TAB

    else:
	s = s + TAB + TAB + TAB

    s = s + r['chromosome'] + TAB + \
	str(r['startC']) + TAB + \
	str(r['endC']) + TAB + \
	str(r['strand']) + TAB

    # if NCBI gene model exists

    if ncbimodel.has_key(key):
	s = s + '%s' % (string.join(ncbimodel[key], ","))
    s = s + TAB

    # if Ensembl gene model exists

    if ensemblmodel.has_key(key):
	s = s + '%s' % (string.join(ensemblmodel[key], ","))
    s = s + TAB

    # if VEGA gene model exists

    if vegamodel.has_key(key):
	s = s + '%s' % (string.join(vegamodel[key], ","))
    s = s + TAB

    s = s + cCDS + CRT

    return s

def printHeaderHTML():

    fpHTML.write('</pre>\n')
    fpHTML.write('<H2>MGI - All Genes List (for candidates to be knocked out by KOMP)</H2>')
    fpHTML.write(introBLOG)
    fpHTML.write('<TABLE BORDER=3 WIDTH=100%>')
    fpHTML.write('<th align=left valign=top>MGI Gene ID</th>')
    fpHTML.write('<th align=left valign=top>Gene Symbol</th>')
    fpHTML.write('<th align=left valign=top>Gene Name</th>')
    fpHTML.write('<th align=left valign=top>KOMP Target</th>')
    fpHTML.write('<th align=left valign=top>Regeneron</th>')
    fpHTML.write('<th align=left valign=top>CSD</th>')
    fpHTML.write('<th align=left valign=top>Other Knockout Source</th>')
    fpHTML.write('<th align=left valign=top>Chr</th>')
    fpHTML.write('<th align=left valign=top>Start Coordinate</th>')
    fpHTML.write('<th align=left valign=top>End Coordinate</th>')
    fpHTML.write('<th align=left valign=top>Strand</th>')
    fpHTML.write('<th align=left valign=top>NCBI ID</th>')
    fpHTML.write('<th align=left valign=top>Ensembl ID</th>')
    fpHTML.write('<th align=left valign=top>VEGA ID</th>')
    fpHTML.write('<th align=left valign=top>cCDS</th>')

def printHeaderTAB():

    fpTAB.write('# MGI - All Genes List (for candidates to be knocked out by KOMP)' + CRT)
    fpTAB.write('# MGI Gene ID' + TAB)
    fpTAB.write('Gene Symbol' + TAB)
    fpTAB.write('Gene Name' + TAB)
    fpTAB.write('KOMP Target' + TAB)
    fpTAB.write('Regeneron' + TAB)
    fpTAB.write('CSD' + TAB)
    fpTAB.write('Other Knockout Source' + TAB)
    fpTAB.write('Chr' + TAB)
    fpTAB.write('Start Coordinate' + TAB)
    fpTAB.write('End Coordinate' + TAB)
    fpTAB.write('Strand' + TAB)
    fpTAB.write('NCBI ID' + TAB)
    fpTAB.write('Ensembl ID' + TAB)
    fpTAB.write('VEGA ID' + TAB)
    fpTAB.write('cCDS' + CRT)

def process():

    global ncbimodel, ensemblmodel, vegamodel

    #
    # genes where MGI gene is associated with at least one of the following:
    #
    #  NCBI EntrezGene
    #  Ensembl
    #  VEGA
    #

    db.sql('select m._Marker_key, m.symbol, name = substring(m.name,1,75), ma.accID ' + \
	    'into #markers ' + \
	    'from MRK_Marker m, ACC_Accession ma ' + \
	    'where m._Organism_key = 1 ' + \
	    'and m._Marker_Type_key = 1 ' + \
	    'and m._Marker_Status_key in (1,3) ' + \
	    'and m.chromosome != "MT" ' + \
	    'and m._Marker_key = ma._Object_key ' + \
	    'and ma._MGIType_key = 2 ' + \
	    'and ma._LogicalDB_key = 1 ' + \
	    'and ma.prefixPart = "MGI:" ' + \
	    'and ma.preferred = 1 ' + \
	    'and exists (select 1 from SEQ_Marker_Cache c where m._Marker_key = c._Marker_key ' + \
	    'and c._LogicalDB_key in (59, 60, 85))', None)
    db.sql('create index idx1 on #markers(_Marker_key)', None)

    #
    # retrieve location information
    #

    db.sql('select m.*, l.chromosome, l.strand, c.sequenceNum, ' + \
	    'startC = convert(int, l.startCoordinate), ' + \
	    'endC = convert(int, l.endCoordinate) ' + \
	    'into #location ' + \
	    'from #markers m, MRK_Location_Cache l, MRK_Chromosome c ' + \
	    'where m._Marker_key = l._Marker_key ' + \
	    'and l.chromosome = c.chromosome ' + \
	    'and c._Organism_key = 1', None)
    db.sql('create index idx1 on #location(_Marker_key)', None)
    db.sql('create index idx2 on #location(sequenceNum)', None)

    #
    # retrieve gene model information
    
    results = db.sql('select m._Marker_key, c._LogicalDB_key, c.accID ' + \
	    'from #markers m, SEQ_Marker_Cache c ' + \
	    'where m._Marker_key = c._Marker_key ' + \
	    'and c._LogicalDB_key in (59, 60, 85)', 'auto')
    for r in results:
        key = r['_Marker_key']
        value = r['accID']

        provider = r['_LogicalDB_key']

        if provider == 59:
          if not ncbimodel.has_key(key):
	      ncbimodel[key] = []
          ncbimodel[key].append(value)

        elif provider == 60:
          if not ensemblmodel.has_key(key):
	      ensemblmodel[key] = []
          ensemblmodel[key].append(value)

        elif provider == 85:
          if not vegamodel.has_key(key):
	      vegamodel[key] = []
          vegamodel[key].append(value)

    #
    # process results
    #
    results = db.sql('select * from #location order by sequenceNum, symbol', 'auto')

    for r in results:
        fpHTML.write('<td>%s</td>' % (writeHTML(r)))
	fpTAB.write(writeTAB(r))

    fpHTML.write('</TABLE>')
    fpHTML.write('<pre>')

    reportlib.finish_nonps(fpHTML, isHTML = 1)	# non-postscript file
    reportlib.finish_nonps(fpTAB)	# non-postscript file

#
# Main
#

init()
process()

