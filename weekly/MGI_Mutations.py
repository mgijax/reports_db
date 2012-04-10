#!/usr/local/bin/python

'''
#
# MGI_Mutations.py
#
# Reports:
#
#       TR 8028
#
#       MGI Mouse Markers:
#       . of type "Gene"
#       . associated with at least one:  
#         NCBI Gene Model
#         Ensembl Gene Model
#         VEGA Gene Model
#	. for which there is one or more public allele
#	  (excluding wild-type (+) alleles and alleles of type=QTL)
#
# History:
#
# lec	12/13/2006
#	- created, TR 8028
#
# jmason 27/08/2007
#   - changed name from MGI_KOMP_Mutations to MGI_Mutations
#   - Updated header text
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
TAB = reportlib.TAB

introBLOG = '''
This report lists those genes that have one or more phenotypic mutations recorded in MGI. Note that this list does not include genes that are defined only through their phenotypic mutation (e.g., cub, curly bare) or those that are not assigned genome coordinates in NCBI Build 36.
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
BEGTD2 = '<td width="25%"><font size="-1">'
BLANKFIELD = '%s&nbsp;%s' % (BEGTD, ENDTD)

MARKER_ANCHOR = '<A HREF="http://www.informatics.jax.org/searches/accession_report.cgi?id=%s" target="_blank">'
CLOSE_ANCHOR = '</A>'

ALLELE_ANCHOR1 = '<A HREF="http://www.informatics.jax.org/searches/allele_report.cgi?_Marker_key=%s&int:_Set_key=%s">'
ALLELE_ANCHOR2 = '<A HREF="http://www.informatics.jax.org/searches/allele_report.cgi?_Marker_key=%s">'

fpHTML = None
fpTAB = None

# key = marker key, value = dictionary of allele categories & counts
markers = {}	

# key = allele type key, value = dictionary of categoryKey that contains 'term' and 'count'
alleleCategories = {}

def init():
    #
    # initialize files
    #

    global fpHTML, fpTAB

    fpHTML = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'], isHTML = 1)
    fpTAB = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])

    printHeaderHTML()
    printHeaderTAB()

def writeHTML(r):
    #
    # write record to HTML file
    #

    # print marker records in HTML format

    key = r['_Marker_key']

    s = '<tr>' + \
	BEGTD + r['accID'] + '</td>\n' + \
	BEGTD + MARKER_ANCHOR % (r['accID']) + r['symbol'] + CLOSE_ANCHOR + ENDTD + \
	BEGTD2 + r['name'] + ENDTD

    s = s + BEGTD

    allAlleles = 0
    for ctermkey in markers[key]:
	allAlleles = allAlleles + markers[key][ctermkey]['count']
    s = s + 'All phenotypic alleles (' + ALLELE_ANCHOR2 % (key) + '%d' % (allAlleles) + CLOSE_ANCHOR + '):'

    for ctermkey in markers[key]:
        s = s + '%s(' % (markers[key][ctermkey]['term']) + \
	    ALLELE_ANCHOR1 % (key, ctermkey) + '%d' % (markers[key][ctermkey]['count']) + CLOSE_ANCHOR + ') '

    s = s + ENDTD + CRT

    return s

def writeTAB(r):
    #
    # write record to tab-delimited file
    #

    # print marker records in tab-delimited format

    key = r['_Marker_key']

    s = r['accID'] + TAB + r['symbol'] + TAB + r['name'] + TAB

    for ctermkey in markers[key]:
        s = s + '%s(%d)|' % (markers[key][ctermkey]['term'], markers[key][ctermkey]['count'])

    s = s + CRT

    return s

def printHeaderHTML():
    #
    # write header to HTML file
    #

    fpHTML.write('</pre>\n')
    fpHTML.write('<H2>MGI - Genes with phenotypic mutations</H2>')
    fpHTML.write(introBLOG)
    fpHTML.write('<TABLE BORDER=3 WIDTH=100%>')
    fpHTML.write('<th align=left valign=top>MGI Gene ID</th>')
    fpHTML.write('<th align=left valign=top>Gene Symbol</th>')
    fpHTML.write('<th align=left valign=top>Gene Name</th>')
    fpHTML.write('<th align=left valign=top>Allele Types(#)</th>')

def printHeaderTAB():
    #
    # write record to tab-delimited file
    #

    fpTAB.write('# MGI - Genes with phenotypic mutations' + CRT)
    fpTAB.write('# MGI Gene ID' + TAB)
    fpTAB.write('Gene Symbol' + TAB)
    fpTAB.write('Gene Name' + TAB)
    fpTAB.write('Allele Types(#)' + CRT)

def process():
    #
    # process data
    #

    global markers, alleleCategories

    #
    # lookup of Marker Detail Allele Categories-to-Allele Type key mappings
    #

    results = db.sql('''
	select vt1._Term_key as categoryKey, vt1.term, vt2._Term_key as alleleTypeKey
        from MGI_VocAssociationType mvat, MGI_VocAssociation mva, VOC_Term vt1, VOC_Term vt2 
        where mvat._AssociationType_key = 1001 
	and mvat._AssociationType_key = mva._AssociationType_key 
	and mva._Term_key_1 = vt1._Term_key
        and mva._Term_key_2 = vt2._Term_key
	''', 'auto')
    for r in results:
	alleleCategories[r['alleleTypeKey']] = r

    #
    # same query as MGI_KOMP_AllGenes
    #
    # genes where MGI gene is associated with at least one of the following:
    #
    #  NCBI EntrezGene
    #  Ensembl
    #  VEGA
    #

    db.sql('''
	select m._Marker_key, m.symbol, substring(m.name,1,75) as name, ma.accID 
	into #markers1 
	from MRK_Marker m, ACC_Accession ma 
	where m._Organism_key = 1 
	and m._Marker_Type_key = 1 
	and m._Marker_Status_key in (1,3) 
	and m.chromosome != "MT" 
	and m._Marker_key = ma._Object_key 
	and ma._MGIType_key = 2 
	and ma._LogicalDB_key = 1 
	and ma.prefixPart = "MGI:" 
	and ma.preferred = 1 
	and exists (select 1 from SEQ_Marker_Cache c where m._Marker_key = c._Marker_key 
	and c._LogicalDB_key in (59, 60, 85))
	''', None)
    db.sql('create index markers1_idx1 on #markers1(_Marker_key)', None)

    #
    # select all genes with mutations
    # exclude wild type (+)
    # exclude QTLs
    # include Approved (public) Alleles only
    #

    db.sql('''
	select m.*, a._Allele_Type_key 
	into #markers2 
	from #markers1 m, ALL_Allele a 
	where m._Marker_key = a._Marker_key 
	and a.isWildType = 0 
	and a._Allele_Type_key != 847130 
	and a._Allele_Status_key = 847114 
	order by m.symbol
	''', None)
    db.sql('create index markers2_idx1 on #markers2(_Marker_key)', None)

    results = db.sql('select * from #markers2', 'auto')

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

    results = db.sql('select distinct m._Marker_key, m.symbol, m.name, m.accID from #markers2 m order by m.symbol', 'auto')
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

