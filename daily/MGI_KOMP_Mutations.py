#!/usr/local/bin/python

'''
#
# MGI_KOMP_Mutations.py
#
# Reports:
#
#       TR 8028
#
#	Markers for which there is one or more allele 
#	(excluding wild-type (+) alleles and alleles of type=QTL)
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
# lec	12/13/2006
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
This report provides a list of those genes considered to be potential candidates for being knocked out through the NIH KOMP project and have one or more phenotypic mutations recorded in MGI. Note that this list does not include genes that are defined only through their phenotypic mutation (e.g., cub, curly bare) or those that are not assigned genome coordinates in NCBI Build 36.
<P>

<UL>
<LI>Each gene symbol is linked to its respective MGI Gene Detail page for additional information.
<P>
<LI>Counts of phenotypic alleles are linked to listings of all mutant alleles for a gene and to subsets of alleles of particular classes (spontaneous, ENU-induced, targeted knock-out, etc.). Following these links leads to phenotypic information. 
<P>
</UL>
'''

BEGTD = '<td><font size="-1">'
ENDTD = '</font></td>'
BLANKFIELD = '%s&nbsp;%s' % (BEGTD, ENDTD)

ALLELE_ANCHOR = '<A HREF="http://www.informatics.jax.org/searches/allele_report.cgi?_Marker_key=%s&int:_Set_key=%s">'
CLOSE_ALLELE_ANCHOR = '</A>'

fpHTML = None
fpTAB = None

# key = marker key, value = dictionary of allele categories & counts
markers = {}	

# key = allele type key, value = dictionary of categoryKey that contains 'term' and 'count'
alleleCategories = {}

def init():

    global fpHTML, fpTAB

    fpHTML = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'], isHTML = 1)
    fpTAB = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])

    printHeaderHTML()
    printHeaderTAB()

def writeHTML(r):

    # print marker records in HTML format

    key = r['_Marker_key']

    s = '<tr>' + \
	BEGTD + r['accID'] + '</td>\n' + \
        '%s%s%s%s%s\n' % (BEGTD, reportlib.create_accession_anchor(r['accID']), r['symbol'], reportlib.close_accession_anchor(), ENDTD) + \
	'<td width="25%"><font size="-1">' + r['name'] + ENDTD

    s = s + BEGTD

    for ctermkey in markers[key]:
        s = s + ALLELE_ANCHOR % (key, ctermkey) + \
		'%s(%d)' % (markers[key][ctermkey]['term'], markers[key][ctermkey]['count']) + \
		CLOSE_ALLELE_ANCHOR + ','

    # drop trailing comma
    s = s[:-1] + ENDTD + CRT

    return s

def writeTAB(r):

    # print marker records in tab-delimited format

    key = r['_Marker_key']

    s = r['accID'] + TAB + r['symbol'] + TAB + r['name'] + TAB

    for ctermkey in markers[key]:
        s = s + '%s(%d),' % (markers[key][ctermkey]['term'], markers[key][ctermkey]['count'])

    # drop trailing comma
    s = s[:-1] + CRT

    return s

def printHeaderHTML():

    fpHTML.write('</pre>\n')
    fpHTML.write('<H2>MGI - Genes with phenotypic mutations</H2>')
    fpHTML.write(introBLOG)
    fpHTML.write('<TABLE BORDER=3 WIDTH=100%>')
    fpHTML.write('<th align=left valign=top>MGI Gene ID</th>')
    fpHTML.write('<th align=left valign=top>Gene Symbol</th>')
    fpHTML.write('<th align=left valign=top>Gene Name</th>')
    fpHTML.write('<th align=left valign=top>Allele Types(#)</th>')

def printHeaderTAB():

    fpTAB.write('# MGI - Genes with phenotypic mutations' + CRT)
    fpTAB.write('# MGI Gene ID' + TAB)
    fpTAB.write('Gene Symbol' + TAB)
    fpTAB.write('Gene Name' + TAB)
    fpTAB.write('Allele Types(#)' + CRT)

def process():

    global markers, alleleCategories

    #
    # lookup of Marker Detail Allele Categories-to-Allele Type key mappings
    #

    results = db.sql('select categoryKey = vt1._Term_key, vt1.term, alleleTypeKey = vt2._Term_key ' + \
           'from MGI_VocAssociationType mvat, MGI_VocAssociation mva, VOC_Term vt1, VOC_Term vt2 ' + \
           'where mvat._AssociationType_key = 1001 ' + \
	   'and mvat._AssociationType_key = mva._AssociationType_key ' + \
	   'and mva._Term_key_1 = vt1._Term_key ' + \
           'and mva._Term_key_2 = vt2._Term_key', 'auto')
    for r in results:
	alleleCategories[r['alleleTypeKey']] = r

    #
    # select all genes with mutations
    # exclude wild type (+)
    # exclude QTLs
    #

    db.sql('select m._Marker_key, a._Allele_Type_key ' + \
	    'into #markers ' + \
	    'from MRK_Marker m, ALL_Allele a ' + \
	    'where m._Organism_key = 1 ' + \
	    'and m._Marker_Type_key = 1 ' + \
	    'and m._Marker_Status_key in (1,3) ' + \
	    'and m._Marker_key = a._Marker_key ' + \
	    'and a.isWildType = 0 ' + \
	    'and a._Allele_Type_key != 847130 ' + \
	    'order by m.symbol', None)
    db.sql('create index idx1 on #markers(_Marker_key)', None)

    results = db.sql('select * from #markers', 'auto')

    for r in results:

	key = r['_Marker_key']
        ctermkey = alleleCategories[r['_Allele_Type_key']]['categoryKey']
	cterm = alleleCategories[r['_Allele_Type_key']]['term']

	#
	# marker "value" = dictionary where key = term key 
	# and values = dictionary of 'term' and 'count'
	#
	# for example:
	#
	# markers[48132] = 
	#
	# {847156: {'count': 2, 'term': 'Targeted, knock-out'}, 
	#  847157: {'count': 1, 'term': 'Targeted, other'}, 
	#  847158: {'count': 3, 'term': 'Gene trapped'}
	# }
	#
	# where 847156 is the term key, 'term' = 'Targeted, knock-out' 
	# and 'count' = the number of alleles of this type
	#

	if not markers.has_key(key):
	    c = {}
        else:
	    c = markers[key]

	if not c.has_key(ctermkey):
	    c[ctermkey] = {}
	    c[ctermkey]['term'] = cterm
	    c[ctermkey]['count'] = 1
        else:
	    c[ctermkey]['count'] = c[ctermkey]['count'] + 1

	markers[key] = c

    results = db.sql('select distinct m._Marker_key, mm.symbol, name = substring(mm.name,1,75), ma.accID ' + \
	    'from #markers m, MRK_Marker mm, ACC_Accession ma ' + \
	    'where m._Marker_key = mm._Marker_key ' + \
	    'and m._Marker_key = ma._Object_key ' + \
	    'and ma._MGIType_key = 2 ' + \
	    'and ma._LogicalDB_key = 1 ' + \
	    'and ma.prefixPart = "MGI:" ' + \
	    'and ma.preferred = 1 ' + \
	    'order by mm.symbol', 'auto')

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

