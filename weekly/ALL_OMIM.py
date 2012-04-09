#!/usr/local/bin/python
'''
#
# tr9711.py
#
# 
# ORIGINAL REQUEST:
#
# For each OMIM disease annotation to a mouse genotype record in the database, include a row:
#
# 1. OMIM ID (can someone link to our disease detail page using an OMIM ID? 
#    if not, then we'd need to add an additional column to the report for the disease term key - yuck!)
# 2. OMIM disease name
# 3. text of the genotype of the model with html markups for links to MGI allele pages.
#    so this should look like the last column on an allele summary page
# 4. (maybe, if easy) text of the genotype of the model without html markups
#    so as above, but remove html markups, include angle brackets around superscripts 
#    (perhaps this stored in the db already?)
#    the purpose of this is to build text indexes from.
# 5. (maybe empty for now?) - the intent is a link to a page that says where someone can 
#    obtain this model. So for us, this should be a link to IMSR, but how do we link to 
#    just this model? Janan, any thoughts?
#
# Had to make some compromises (see TR progress report section), so here are the actual columns:
#   1. OMIM id
#   2. OMIM term
#   3. database key for OMIM term
#   4. strain (plain text)
#   5. allele pairs (hyperlinks)
#   6. allele pairs (plain text) 
# 
# History:
#
# jer	07/02/2009
#	- created
# jer	09/24/2009
#	- modified so that wild type alleles are NOT hyperlinked in col. 5
# jer	09/25/2009
#	- modified to make it a weekly report
#
'''
 
import sys
import os
import re
import db
import reportlib

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#--------------------------------------------------------------------------------------------
# Most of the work is done with regular expression matching and replacement on values
# found in the genotypeDisplay1 and genotypeDisplay2 columns of the MRK_OMIM_Cache table.
#
# This regex matches mutant alleles, which are coded like:
#	\Allele(MGI:1931521|Fgfr3<tm1Dor>|)
# It captures the MGI id in group 1 and the allele symbol in group 2
#
allele_re1 = re.compile(r'\\Allele\((MGI:[0-9]+)\|([^|]*)\|\)')

#
# This regex matches wild type alleles, which are coded like:
#	\AlleleSymbol(MGI:2178355|0) 
# It captures the MGI id in group 1. 
#
allele_re2 = re.compile(r'\\AlleleSymbol\((MGI:[0-9]+)\|0\)')

#
# This regex matches the allele portion of an allele symbol, i.e., the stuff
# between the angle brackets. That part is captured in group 1.
#
allele_re3 = re.compile(r'<([^>]*)>')

#
# This is a replacement pattern for creating an HTML link from matches against allele_re1
# 
link_repl   = r'<a href="http://www.informatics.jax.org/javawi2/servlet/WIFetch?page=alleleDetail&id=\1">\2</a>' 

#--------------------------------------------------------------------------------------------
fp = None
nrecs = 0
#--------------------------------------------------------------------------------------------

def generate(r):
    global fp
    global nrecs

    # glue the two genotypeDisplay fields together, if the second one has a value
    g = r['genotypeDisplay1']
    if r['genotypeDisplay2'] is not None:
	g+=r['genotypeDisplay2']
    g = g.strip().replace(CRT, SPACE)

    # plain text
    g1 = allele_re1.sub(r'\2',allele_re2.sub('+',g))

    # HTML, including links
    g2 = allele_re1.sub(link_repl, 
	    allele_re2.sub('+',
		allele_re3.sub(r'<sup>\1</sup>',g)))

    orec = [
	r['termID'],
	r['term'],
	str(r['_term_key']),
	r['strain'],
	g2,
	g1,
	]
    fp.write( TAB.join(orec) + CRT )
    nrecs += 1

def main():
    global fp
    global nrecs

    fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

    cmd = '''
	select distinct 
	    moc.term, moc.termID, moc._term_key, moc._genotype_key, 
	    moc.genotypeDisplay1, moc.genotypeDisplay2, moc.strain
	from mrk_omim_cache moc
	where moc._genotype_key is not null
	and moc.genotypeDisplay1 is not null
	and moc._organism_key = 1
	and moc.qualifier is null
	order by moc.term, moc._genotype_key
	'''

    results = db.sql(cmd, 'auto')
    for r in results:
	generate(r)

    reportlib.finish_nonps(fp)

###

main()
