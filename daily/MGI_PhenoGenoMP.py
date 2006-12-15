#!/usr/local/bin/python

'''
#
# MGI_PhenoGeno.py
#
# Report:
#       Phenotype/Genotype MP Annotations
#
#	Tab-delimited report where a row in the report
#	represents a Genotype annotation to a MP term
#	and all of the references (J:) used to annotate that term.
#
#	field 1: MGI Marker ID
#	field 2: Allele Combination
#	field 3: Strain
#	field 4: MP ID
#	field 5: comma-delimited set of J:
#
# Usage:
#       MGI_PhenoGenoMP.py
#
# History:
#
# lec	10/04/2005
#	- TR 5188; GO Qualifier
#
# lec	06/03/2005
#	- created
#
'''
 
import sys 
import os
import db
import re
import string
import reportlib

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# select all MP annotations
#

db.sql('select distinct a._Object_key, a._Term_key, e._Refs_key ' + 
	'into #mp ' + \
	'from VOC_Annot a, VOC_Evidence e, VOC_Term t ' + \
	'where a._AnnotType_key = 1002 ' + \
	'and a._Qualifier_key = t._Term_key ' + \
	'and t.term is null ' + \
	'and a._Annot_key = e._Annot_key', None)
db.sql('create index idx1 on #mp(_Object_key)', None)
db.sql('create index idx2 on #mp(_Term_key)', None)
db.sql('create index idx3 on #mp(_Refs_key)', None)

#
# resolve MP ids
#
results = db.sql('select distinct m._Term_key, a.accID from #mp m, ACC_Accession a ' + \
	'where m._Term_key = a._Object_key ' + \
	'and a._MGIType_key = 13 ' + \
	'and a.preferred = 1', 'auto')
mpID = {}
for r in results:
    key = r['_Term_key']
    value = r['accID']
    mpID[key] = value
    
#
# resolve References
#
results = db.sql('select distinct m._Object_key, m._Term_key, a.accID from #mp m, ACC_Accession a ' + \
	'where m._Refs_key = a._Object_key ' + \
	'and a._LogicalDB_key = 29 ', 'auto')
mpRef = {}
for r in results:
    key = `r['_Object_key']` + `r['_Term_key']`
    value = r['accID']
    if not mpRef.has_key(key):
	mpRef[key] = []
    mpRef[key].append(value)

#
# resolve Genotype Strains
#
results = db.sql('select distinct m._Object_key, s.strain from #mp m, GXD_Genotype g, PRB_Strain s ' + \
	'where m._Object_key = g._Genotype_key ' + \
	'and g._Strain_key = s._Strain_key', 'auto')
mpStrain = {}
for r in results:
    key = r['_Object_key']
    value = r['strain']
    mpStrain[key] = value

#
# resolve Allele Combination display
#
results = db.sql('select distinct m._Object_key, nc.note from #mp m, MGI_Note n, MGI_NoteChunk nc ' + \
	'where m._Object_key = n._Object_key ' + \
	'and n._NoteType_key = 1016 ' + \
	'and n._Note_key = nc._Note_key', 'auto')
mpDisplay = {}
for r in results:
    key = r['_Object_key']
    value = re.sub('\n', ',', string.strip(r['note']))
    mpDisplay[key] = value

#
# resolve Marker ID
#
results = db.sql('select distinct m._Object_key, a.accID from #mp m, GXD_AlleleGenotype g, ACC_Accession a ' + \
	'where m._Object_key = g._Genotype_key ' + \
	'and g._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1', 'auto')
mpMarker = {}
for r in results:
    key = r['_Object_key']
    value = r['accID']
    if not mpMarker.has_key(key):
	mpMarker[key] = []
    mpMarker[key].append(value)

#
# process results
#
results = db.sql('select distinct _Object_key, _Term_key from #mp order by _Object_key, _Term_key', 'auto')

for r in results:

    genotype = r['_Object_key']
    term = r['_Term_key']
    refKey = `genotype` + `term`

    # we only want to list the Genotypes that have Allele Pairs
    if mpDisplay.has_key(genotype):

        fp.write(mpDisplay[genotype] + TAB + \
	    mpStrain[genotype] + TAB + \
	    mpID[term] + TAB)
        if mpRef.has_key(refKey):
	    fp.write(string.join(mpRef[refKey], ','))
	fp.write(TAB + string.join(mpMarker[genotype], ',') + CRT)

reportlib.finish_nonps(fp)	# non-postscript file

