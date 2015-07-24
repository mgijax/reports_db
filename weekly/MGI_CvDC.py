#!/usr/local/bin/python

'''
#
#
# TR11134
# 
# Ellen Devlin (edevlin@NERIscience.com) 
# has requested a custom SQL report of data generated and reported to MGI as part of the 
# Cardiovascular Development Consortium. This is expected to become a weekly run report 
# once we have finalized the content and format.
#
# 1. Data Required:
# 
# A tab-delimited report of allele and phenotype data for CvDC. 
# Currently, alleles associated with this project have been assigned reference 
# number J:175213 and this may be used to extract these records. However, it is 
# anticipated that additional alleles with other reference numbers will be added 
# to this report at a later date. We will consider how best to extract these at that time.
#
# 1. Allele Symbol
# 2. Allele MGI ID
# 3. Allele Name
# 4. Allele Synonyms (command delimited)
# 5. Allele Type
# 6. Gene Symbol
# 7. Gene MGI ID
# 8. Gene Name
# 9. Chromosome (number or UN)
# 10. Summative Diagnosis text 
#	This text is located in the General Notes field.
# 11. Fyler Codes table 
#	This is located in the General Notes field.
# 12. MP IDs and terms associated with genotypes for the allele in column 1 (comma delimited)
# 13. JAX Strain Registry ID
#	if present. This ID found in the strains table; 
#	please list the JAX Registry number for Public strains 
#	(i.e. "Private?" ="No") that are associated with the allele in column 1. 
#
# History
#
# lec   05/12/2015
#       - TR12020/use ALL_Allele._Marker_key
#
# lec	05/01/2014
#	- TR11669/missing gene info
#
# lec	08/02/2012
#	- TR11134
#
'''
 
import sys 
import os
import mgi_utils
import reportlib
import db

db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE()

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

db.useOneConnection(1)

#
# remove old *.html files
#
os.system('rm -rf ${CVDCDIR}/*.html')

fp = reportlib.init(sys.argv[0], outputdir = os.environ['CVDCDIR'], printHeading = None)

fp.write('# 1. Allele Symbol' + CRT)
fp.write('# 2. Allele MGI ID' + CRT)
fp.write('# 3. Allele Name' + CRT)
fp.write('# 4. Allele Synonyms ' + CRT)
fp.write('# 5. Allele Type' + CRT)
fp.write('# 6. Gene Symbol' + CRT)
fp.write('# 7. Gene MGI ID' + CRT)
fp.write('# 8. Gene Name' + CRT)
fp.write('# 9. Chromosome ' + CRT)
fp.write('# 10. MP IDs and terms ' + CRT)
fp.write('# 11. JAX Strain Registry ID' + 2*CRT)

db.sql('''
	select distinct a._Allele_key, a._Marker_key, a.symbol, a.name, aa.accID, t.term as alleletype
	into temporary table allele
	from ALL_Allele a, ACC_Accession aa, VOC_Term t, MGI_Reference_Assoc r
	where a._Allele_Status_key in (847114)
	and a._Allele_Type_key = t._Term_key
	and a._Allele_key = aa._Object_key
	and aa._MGIType_key = 11
	and aa._LogicalDB_key = 1
	and aa.prefixPart = 'MGI:'
	and aa.preferred = 1
	and a._Allele_key = r._Object_key
	and r._MGIType_key = 11
	and r._Refs_key in (176309, 186964)
	''', None)
db.sql('create index idx_allele on allele(_Allele_key)', None)

#
# genes
#

results = db.sql('''
	select a._Allele_key, m.symbol, m.name, m.chromosome, aa.accID
	from allele a, MRK_Marker m, ACC_Accession aa
	where a._Marker_key = m._Marker_key
	and m._Marker_key = aa._Object_key
	and aa._MGIType_key = 2
	and aa._LogicalDB_key = 1
	and aa.prefixPart = 'MGI:'
	and aa.preferred = 1
	''', 'auto')

markers = {}
for r in results:
    key = r['_Allele_key']
    value = r
    if not markers.has_key(key):
	markers[key] = []
    markers[key].append(value)

#
# synonyms
#
results = db.sql('''
	select a._Allele_key, s.synonym 
    	from allele a, MGI_Synonym s
    	where a._Allele_key = s._Object_key 
    	and s._MGIType_key = 11
	''', 'auto')
syns = {}
for r in results:
    key = r['_Allele_key']
    value = r['synonym']
    if not syns.has_key(key):
        syns[key] = []
    syns[key].append(value)

#
# mgi ids & terms
#

results = db.sql('''
        select distinct a._Allele_key, t.term, aa.accID 
        from allele a, GXD_AlleleGenotype ga, VOC_Annot na, VOC_Term t, ACC_Accession aa 
        where a._Allele_key = ga._Allele_key 
        and ga._Genotype_key = na._Object_key 
        and na._AnnotType_key = 1002 
	and na._Term_key = t._Term_key
        and na._Term_key = aa._Object_key 
        and aa._MGIType_key = 13 
        and aa.preferred = 1
	order by aa.accID, t.term
        ''', 'auto')
mpIDs = {}
for r in results:
	key = r['_Allele_key']
	value = r
        if not mpIDs.has_key(key):
                mpIDs[key] = []
        mpIDs[key].append(value)

#
# summative diagnosis text
# fluer codes
# general notes (1020)
#
results = db.sql('''
        select a._Allele_key, nc.note 
        from allele a, MGI_Note n, MGI_NoteChunk nc 
        where a._Allele_key = n._Object_key 
        and n._MGIType_key = 11 
	and n._NoteType_key = 1020
        and n._Note_key = nc._Note_key 
        order by a._Allele_key, nc.sequenceNum
        ''', 'auto')
notes = {}
for r in results:
    key = r['_Allele_key']
    value = r['note']
    if notes.has_key(key):
	notes[key] = notes[key] + value
    else:
        notes[key] = value

#
# jax strain registry id
#
results = db.sql('''
	select a._Allele_key, aa.accID
    	from allele a, PRB_Strain_Marker p, ACC_Accession aa
    	where a._Allele_key = p._Allele_key
	and p._Strain_key = aa._Object_key
	and aa._MGIType_key = 10
	and aa._LogicalDB_key = 22
	''', 'auto')
jaxID = {}
for r in results:
    key = r['_Allele_key']
    value = r['accID']
    if not jaxID.has_key(key):
        jaxID[key] = []
    jaxID[key].append(value)

#
# main select
#
results = db.sql('select * from allele', 'auto')

for r in results:
    key = r['_Allele_key']

# 1. Allele Symbol
# 2. Allele MGI ID
# 3. Allele Name
    fp.write(r['symbol'] + TAB)
    fp.write(r['accID'] + TAB)
    fp.write(r['name'] + TAB)

# 4. Allele Synonyms (command delimited)

    if syns.has_key(key):
	fp.write('|'.join(syns[key]))
    fp.write(TAB)

# 5. Allele Type
    fp.write(r['alleletype'] + TAB)

# 6. Gene Symbol
# 7. Gene MGI ID
# 8. Gene Name
# 9. Chromosome (number or UN)
    if markers.has_key(key):
	for m in markers[key]:
	    fp.write(m['symbol'] + TAB)
	    fp.write(m['accID'] + TAB)
	    fp.write(m['name'] + TAB)
	    fp.write(m['chromosome'] + TAB)
    else:
	fp.write(TAB + TAB + TAB + TAB)

    #
    # Summative Diagnosis
    #
    #if notes.has_key(key):
#	n = notes[key]
#	spattern = n.find('<! start SD')
#	epattern = n.find('<! end SD')
#	fp.write(n[spattern + 12:epattern])
#    fp.write(TAB)

    #
    # Fyler Codes
    #
#    if notes.has_key(key):
#	n = notes[key]
#	fpattern = re.compile('<td>[0-9][0-9][0-9][0-9]</td>')
#	for f in fpattern.findall(n):
#	    fyler = f.replace('<td>', '')
#	    fyler = fyler.replace('</td>', '')
#	    fp.write(fyler + '|')
#    fp.write(TAB)

    if mpIDs.has_key(key):
	for m in mpIDs[key]:
	    fp.write(m['accID'] + '|' + m['term'] + '||')
    fp.write(TAB)

    if jaxID.has_key(key):
	fp.write('|'.join(jaxID[key]))
    fp.write(CRT)

    if notes.has_key(key):
        fphtml = reportlib.init(r['accID'], fileExt = '.html', outputdir = os.environ['CVDCDIR'], printHeading = None)
	fphtml.write(notes[key] + CRT)
        reportlib.finish_nonps(fphtml)	# non-postscript file

reportlib.finish_nonps(fp)	# non-postscript file
db.useOneConnection(0)

