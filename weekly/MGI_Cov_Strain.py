'''
#
# MGI_Cov_Strain.py
#
# Report:
#       TR13459
#       http://wts.informatics.jax.org/searches/tr.detail.cgi?TR_Nr=TR13459
#
# Requirements:
#
# allele-exclude set:
#     allele type = 'Transgenic' AND subtypes only in: 'Reporter', 'Transactivator', 'Recombinase', 'Inducible'
#
# strain by strain-reference associations
#     strain is associated with covid-reference (mgi_reference_assoc)
# 
# for each Strain:
#    -> if no genotype OR
#    -> if genotype
#        -> if allele
#           -> allele is not wild-type
#           -> allele not in alleleExclude
#        -> genotype with MP header OR
#        -> genotype with DO annotation
#
# Output format:
#
# 1/A) Strain name
# 2/B) Strain MGI ID 
# 3/C) Strain attributes (pipe delimited) 
# 4/D) Genotype MGI ID
# 5/E) Allele symbol
# 6/F) Allele MGI ID 
# 7/G) Allele type 
# 8/H) Allele subtypes  (pipe delimited)
# 9/I) DO term
# 10/J) DO ID 
# 11/K) MP term
# 12/L) MP ID
# 13/M) Reference short citation 
# 14/N) Reference JnumID 
# 15/O) Reference-associated COV tags (pipe delimited) 
#
# History:
#
# lec   02/02/2021
#       - TR13459/Coronavirus Reports for Jax Web Resources
#
'''

import sys
import os
import reportlib
import db

db.setTrace()

TAB = reportlib.TAB
CRT = reportlib.CRT

covidTags = {}
strainAttrs = {}
mpGenotypes = {}
doGenotypes = {}
allelesByGenotype = {}
alleleSubtypes = {}

#
# select covid references
# select covid tags
#
def initializeRefSet():
    global covidTags

    #
    # covid references
    #
    cmd = '''
    select distinct r._refs_key, r.jnumid, r.short_citation
    into temp table covidRefs
    from bib_citation_cache r
    where exists (select 1 from bib_workflow_tag tg, voc_term tt
            where r._refs_key = tg._refs_key
            and tg._tag_key = tt._term_key
            and tt._vocab_key = 129
            and tt.term ilike 'COV:%'
            )
    '''
    db.sql(cmd, None)
    db.sql('create index covidKey on covidRefs(_refs_key)', None)

    #
    # covid tags lookup
    #
    cmd = '''
    select distinct r._refs_key, tt.term as covidtag
    from covidRefs r, bib_workflow_tag tg, voc_term tt
    where r._refs_key = tg._refs_key
    and tg._tag_key = tt._term_key
    and tt._vocab_key = 129
    and tt.term ilike 'COV:%'
    '''
    results = db.sql(cmd, 'auto')
    for r in results:
        key = r['_refs_key']
        value = r['covidtag']
        if key not in covidTags:
            covidTags[key] = []
        covidTags[key].append(value)
    #print(covidTags)

#
# allele exclusions
#     allele type = 'Transgenic' AND subtypes only in: 'Reporter', 'Transactivator', 'Recombinase', 'Inducible'
#
# genotype exclusion
#     genotype that has > 1 marker
#
def initializeAlleleExclude():

    cmd = '''
    select distinct a._allele_key
    into temp table alleleExclude
    from all_allele a
    where a._allele_type_key = 847126
    and exists (select 1 from voc_annot v
            where a._allele_key = v._object_key
            and v._annottype_key = 1014
            and v._term_key in (11025588,11025589,13289567,11025592)
            )
    and not exists (select 1 from voc_annot v
            where a._allele_key = v._object_key
            and v._annottype_key = 1014
            and v._term_key not in (11025588,11025589,13289567,11025592)
            )
    '''
    db.sql(cmd, None)
    db.sql('create index alleleKey1 on alleleExclude(_allele_key)', None)

    cmd = '''
    select distinct g._genotype_key, g._marker_key
    into temp table genotypeExclude
    from gxd_allelegenotype g
    group by _genotype_key, _marker_key having count(*) > 1
    '''
    db.sql(cmd, None)
    db.sql('create index genotypeKey1 on genotypeExclude(_genotype_key)', None)

#
# strain by strain-reference associations
#
def initializeStrain():
    global allelesByGenotype

    #
    # strain set
    #
    cmd = '''
    select distinct s._strain_key, s.strain, a.accid as strainid, cv.*
    into temp table strainSet
    from prb_strain s, acc_accession a, covidRefs cv, mgi_reference_assoc r
    where s._strain_key = a._object_key
    and a._mgitype_key = 10
    and a._logicaldb_key = 1
    and a.preferred = 1
    and s._strain_key = r._Object_key
    and cv._refs_key = r._refs_key
    and r._mgitype_key = 10
    ''' 

    db.sql(cmd, None)
    db.sql('create index strainKey on strainSet(_strain_key)', None)
    db.sql('create index refsKey1 on strainSet(_refs_key)', None)

    #
    # alleles by genotype
    #
    allelesByGenotype = {}
    cmd = '''
    select distinct g._genotype_key, a.accid as alleleid, aa._allele_key, aa.symbol, t.term as alleletype
    into temp table alleleSet
    from strainSet s, gxd_genotype g, gxd_allelegenotype ag, acc_accession a, all_allele aa, voc_term t
    where s._strain_key = g._strain_key
    and g._genotype_key = ag._genotype_key
    and ag._allele_key = a._object_key
    and a._mgitype_key = 11
    and a._logicaldb_key = 1
    and a.preferred = 1
    and ag._allele_key = aa._allele_key
    and aa.iswildtype = 0
    and aa._allele_type_key = t._term_key
    and not exists (select 1 from alleleExclude e where aa._allele_key = e._allele_key)
    '''
    db.sql(cmd, None)
    db.sql('create index alleleKey2 on alleleSet(_allele_key)', None)
    results = db.sql('select * from alleleSet', 'auto')
    for r in results:
        key = r['_genotype_key']
        value = r
        if key not in allelesByGenotype:
                allelesByGenotype[key] = []
        allelesByGenotype[key].append(value)
    #print(allelesByGenotype)

    #
    # alleleSubtype lookup
    #
    cmd = '''
    select distinct ag._allele_key, t.term
    from strainSet s, gxd_genotype g, gxd_allelegenotype ag, alleleSet aa, voc_annot v, voc_term t
    where s._strain_key = g._strain_key
    and g._genotype_key = ag._genotype_key
    and ag._allele_key = aa._allele_key
    and ag._allele_key = v._object_key
    and v._annottype_key = 1014
    and v._term_key = t._term_key
    and not exists (select 1 from alleleExclude e where aa._allele_key = e._allele_key)
    '''
    results = db.sql(cmd, 'auto')
    for r in results:
        key = r['_allele_key']
        value = r['term']
        if key not in alleleSubtypes:
                alleleSubtypes[key] = []
        alleleSubtypes[key].append(value)
    #print(alleleSubtypes)

    initializeStrainGeneral()

#
# strain attributes
# mpGenotypes with mp annotations with covid reference
# doGenotypes with DO annotations with at least 1 covid reference
#
def initializeStrainGeneral():
    global strainAttrs
    global mpGenotypes
    global doGenotypes
    global alleleSubtypes

    #
    # strain attributes
    #
    cmd = '''
    select distinct s._strain_key, t.term as strainattr
    from strainSet s, voc_annot va, voc_term t
    where s._strain_key = va._object_key
    and va._annotType_key = 1009
    and va._term_key not in (481372, 481373)
    and va._term_key = t._term_key
    '''
    results = db.sql(cmd, 'auto')
    for r in results:
        key = r['_strain_key']
        value = r['strainattr']
        if key not in strainAttrs:
            strainAttrs[key] = []
        strainAttrs[key].append(value)
    #print(strainAttrs)

    #
    # mpGenotypes with mp annotations with covid reference
    #
    cmd = '''
    select s._strain_key, s._refs_key, g._genotype_key, a.accid as genotypeid, 
        t.label, t._object_key, tt.term as otherLabel, mp.accid as mpid
    from strainSet s, gxd_genotype g, acc_accession a,
            voc_annotheader h, mgi_setmember t, acc_accession mp, voc_term tt
    where s._strain_key = g._strain_key
    and g._genotype_key = a._object_key
    and a._mgitype_key = 12
    and a._logicaldb_key = 1
    and a.preferred = 1
    and g._genotype_key = h._object_key
    and h._annottype_key = 1002
    and h._term_key = t._object_key
    and t._set_key = 1051
    and t._object_key = mp._object_key
    and mp._mgitype_key = 13
    and mp._logicaldb_key = 34
    and mp.preferred = 1
    and t._object_key = tt._term_key
    and exists (select 1 from voc_annot v, voc_evidence e
            where h._object_key = v._object_key
            and v._annottype_key = 1002
            and v._annot_key = e._annot_key
            and e._refs_key = s._refs_key
            )
    order by s._strain_key, s._refs_key, genotypeid,  t.label
    '''
    results = db.sql(cmd, 'auto')
    for r in results:
        if r['label'] == None:
            r['label'] = r['otherLabel']
        key = r['_strain_key'] + r['_refs_key']
        value = r
        if key not in mpGenotypes:
            mpGenotypes[key] = []
        mpGenotypes[key].append(value)
    #print(mpGenotypes)

    #
    # doGenotypes with DO annotations with at least 1 covid reference
    #
    cmd = '''
    select s._strain_key, s._refs_key, g._genotype_key, a.accid as genotypeid, 
        tt.term, aa.accid as doid
    from strainSet s, gxd_genotype g, acc_accession a, acc_accession aa, voc_annot vt, voc_term tt
    where s._strain_key = g._strain_key
    and g._genotype_key = a._object_key
    and a._mgitype_key = 12
    and a._logicaldb_key = 1
    and a.preferred = 1
    and g._genotype_key = vt._object_key
    and vt._annottype_key = 1020
    and vt._term_key = tt._term_key
    and vt._term_key = aa._object_key
    and aa._mgitype_key = 13
    and aa._logicaldb_key = 191
    and aa.preferred = 1
    and exists (select 1 from voc_annot v, voc_evidence e
            where g._genotype_key = v._object_key
            and v._annottype_key = 1020
            and v._annot_key = e._annot_key
            and e._refs_key = s._refs_key
            )
    order by s._strain_key, s._refs_key, genotypeid, tt.term
    '''
    results = db.sql(cmd, 'auto')
    for r in results:
        key = r['_strain_key'] + r['_refs_key']
        value = r
        if key not in doGenotypes:
            doGenotypes[key] = []
        doGenotypes[key].append(value)
    #print(doGenotypes)

#
# mp annotations
#
def processMP(r, key):

    strain = r['strain']
    strainKey = r['_strain_key']
    refsKey = r['_refs_key']

    for g in mpGenotypes[key]:

        gKey = g['_genotype_key']

        # no alleles
        if gKey not in allelesByGenotype:

            fp.write(strain + TAB)
            fp.write(r['strainid'] + TAB)

            if strainKey in strainAttrs:
                fp.write('|'.join(strainAttrs[strainKey]) + TAB)
            else:
                fp.write(TAB)

            fp.write(g['genotypeid'] + TAB)
            fp.write(TAB*6)
            fp.write(g['label'] + TAB)
            fp.write(g['mpid'] + TAB)
            fp.write(r['short_citation'] + TAB)
            fp.write(r['jnumid'] + TAB)
            fp.write('|'.join(covidTags[refsKey]) + CRT)

            continue
                
        for a in allelesByGenotype[gKey]:
            aKey = a['_allele_key']

            fp.write(strain + TAB)
            fp.write(r['strainid'] + TAB)

            if strainKey in strainAttrs:
                fp.write('|'.join(strainAttrs[strainKey]) + TAB)
            else:
                fp.write(TAB)

            fp.write(g['genotypeid'] + TAB)
            fp.write(a['symbol'] + TAB)
            fp.write(a['alleleid'] + TAB)
            fp.write(a['alleletype'] + TAB)

            if aKey in alleleSubtypes:
                fp.write('|'.join(alleleSubtypes[aKey]) + TAB)
            else:
                fp.write(TAB)

            fp.write(TAB*2)
            fp.write(g['label'] + TAB)
            fp.write(g['mpid'] + TAB)
            fp.write(r['short_citation'] + TAB)
            fp.write(r['jnumid'] + TAB)
            fp.write('|'.join(covidTags[refsKey]) + CRT)

#
# do annotations
#
def processDO(r, key):

    strain = r['strain']
    strainKey = r['_strain_key']
    refsKey = r['_refs_key']

    for g in doGenotypes[key]:

        gKey = g['_genotype_key']

        # no alleles
        if gKey not in allelesByGenotype:
            fp.write(strain + TAB)
            fp.write(r['strainid'] + TAB)

            if strainKey in strainAttrs:
                fp.write('|'.join(strainAttrs[strainKey]) + TAB)
            else:
                fp.write(TAB)

            fp.write(g['genotypeid'] + TAB)
            fp.write(TAB*4)
            fp.write(g['term'] + TAB)
            fp.write(g['doid'] + TAB)
            fp.write(TAB*2)
            fp.write(r['short_citation'] + TAB)
            fp.write(r['jnumid'] + TAB)
            fp.write('|'.join(covidTags[refsKey]) + CRT)
            continue
                
        for a in allelesByGenotype[gKey]:
            aKey = a['_allele_key']

            fp.write(strain + TAB)
            fp.write(r['strainid'] + TAB)

            if strainKey in strainAttrs:
                fp.write('|'.join(strainAttrs[strainKey]) + TAB)
            else:
                fp.write(TAB)

            fp.write(g['genotypeid'] + TAB)
            fp.write(a['symbol'] + TAB)
            fp.write(a['alleleid'] + TAB)
            fp.write(a['alleletype'] + TAB)

            if aKey in alleleSubtypes:
                fp.write('|'.join(alleleSubtypes[aKey]) + TAB)
            else:
                fp.write(TAB)

            fp.write(g['term'] + TAB)
            fp.write(g['doid'] + TAB)
            fp.write(TAB*2)
            fp.write(r['short_citation'] + TAB)
            fp.write(r['jnumid'] + TAB)
            fp.write('|'.join(covidTags[refsKey]) + CRT)

#
# process strain by allele-reference associations
#
def processStrain():

    #
    # strainSet by strain
    #
    results = db.sql('select * from strainSet order by strain, short_citation', 'auto')
    for r in results:

        strain = r['strain']
        strainKey = r['_strain_key']
        refsKey = r['_refs_key']
        key = strainKey + refsKey

        # no genotypes with mp or do annotations
        if key not in mpGenotypes and key not in doGenotypes:
            fp.write(r['strain'] + TAB)
            fp.write(r['strainid'] + TAB)

            if strainKey in strainAttrs:
                fp.write('|'.join(strainAttrs[strainKey]) + TAB)
            else:
                fp.write(TAB)

            fp.write(TAB*9)
            fp.write(r['short_citation'] + TAB)
            fp.write(r['jnumid'] + TAB)
            fp.write('|'.join(covidTags[refsKey]) + CRT)
            continue

        # MP annotations
        if key in mpGenotypes:
           processMP(r, key)

        # DO annotations
        if key in doGenotypes:
           processDO(r, key)

#
# main
#

db.useOneConnection(1)

fp = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])

fp.write('#\n')
fp.write('#1/A) Strain name\n')
fp.write('#2/B) Strain MGI ID \n')
fp.write('#3/C) Strain attributes (pipe delimited) \n')
fp.write('#4/D) Genotype MGI ID\n')
fp.write('#5/E) Allele symbol\n')
fp.write('#6/F) Allele MGI ID \n')
fp.write('#7/G) Allele type \n')
fp.write('#8/H) Allele subtypes (pipe delimited) \n')
fp.write('#9/I) DO term\n')
fp.write('#10/J) DO ID \n')
fp.write('#11/K) MP term\n')
fp.write('#12/L) MP ID\n')
fp.write('#13/M) Reference short citation \n')
fp.write('#14/N) Reference JnumID \n')
fp.write('#15/O) Reference-associated COV tags (pipe delimited) \n')
fp.write('#\n')

initializeRefSet()
initializeAlleleExclude()
initializeStrain()
processStrain()

reportlib.finish_nonps(fp)
db.useOneConnection(0)

