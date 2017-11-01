#!/usr/local/bin/python

'''
#
# MGI_Pheno_Sex.py
#
# Report:
#       Tab-delimited file of Phenotypes where males and females are different
#
# Usage:
#       MGI_Pheno_Sex.py
#
# Used by:
#
#	Peter Robinson/JAX	
#
# Output format:
#
# History:
#
# sc	10/31/2017
#	- TR12632 Report of phenotypes with differences between male and female
#
'''

import sys
import os
import string
import mgi_utils
import reportlib
import db

db.setTrace()

TAB = reportlib.TAB
CRT = reportlib.CRT
#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
results = db.sql('''select a._Object_key as _Genotype_key, a2.accid as genotypeID, a._Term_key, t1.term as mpTerm,  p.value as gender, a1.accid as mpID, t2.term as annotQual, s.strain as bgStrain, nc.note as alleleComp
    from VOC_Annot a, VOC_Evidence e, VOC_Evidence_Property p, VOC_Term t1, VOC_Term t2, ACC_Accession a1, ACC_Accession a2, GXD_Genotype g, PRB_Strain s,  MGI_Note n, MGI_NoteChunk nc
    where a._AnnotType_key = 1002 -- MP
    and a._Annot_key = e._Annot_key
    and e._AnnotEvidence_key = p._AnnotEvidence_key
    and p._PropertyTerm_key = 8836535 --vocabKey 86, MP-Sex-Specificity
    and lower(p.value) in ('m', 'f') -- most upper, a few are lower
    and a._Qualifier_key = t2._Term_key
    and a._Term_key = t1._Term_key
    and a._Term_key = a1._Object_key
    and a1._MGIType_key = 13 -- VOC_Term
    and a1._LogicalDB_key = 34 -- MP
    and a1.preferred = 1
    and a._Object_key = a2._Object_key
    and a2._MGIType_key = 12 -- GXD_Genotype
    and a2._LogicalDB_key = 1 -- MGI
    and a2.prefixPart = 'MGI:'
    and a2.preferred = 1
    and a._Object_key = g._Genotype_key
    and g._Strain_key = s._Strain_key -- background strain
    and g._Genotype_key = n._Object_key
    and n._NOtetype_key = 1016
    and n._Note_key = nc._Note_key
    order by a._Object_key''', 'auto')

fp.write('Genotype ID%sSex%sMP ID%sMP Term%sAllelic Composition%sBackground Strain%sSex-specific Normal Y/N%s' % (TAB, TAB, TAB, TAB, TAB, TAB, CRT))
for r in results:
    annotQual = 'N'
    if r['annotQual'] != None:
	annotQual = 'Y'
    alleleComp = string.strip(r['alleleComp'])
    fp.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (r['genotypeID'], TAB, r['gender'], TAB, r['mpID'], TAB, r['mpTerm'], TAB, alleleComp, TAB, r['bgStrain'], TAB, annotQual, CRT)) 

reportlib.finish_nonps(fp)
db.useOneConnection(0)

