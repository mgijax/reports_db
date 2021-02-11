'''
#
# MGI_Cov_StrainAllele.py
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
# allele by allele-reference associations
#     allele is associated with covid-reference (mgi_reference_assoc) and reference is not in by-strain set
#     allele is associated with a genotype (gxd_allelegenotype)
#     allele is not wild-type
#     allele not in alleleExclude
#     genotype with 1 marker
#
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
def initializeByStrain():
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
# allele by allele-reference associations
#
def initializeByAllele():
    global allelesByGenotype
    global alleleSubtype

    #
    # allele is associated with covid-reference (mgi_reference_assoc)
    # allele is associated with a genotype (gxd_allelegenotype)
    # allele is not wild-type
    # allele not in alleleExclude
    #
    db.sql('drop table if exists alleleSet', None)
    cmd = '''
    select distinct aa._allele_key, aa.symbol, a.accid as alleleid, t.term as alleletype, cv.*
    into temp table alleleSet
    from all_allele aa, acc_accession a, covidRefs cv, mgi_reference_assoc r, voc_term t
    where aa._allele_key = a._object_key
    and a._mgitype_key = 11
    and a._logicaldb_key = 1
    and a.preferred = 1
    and aa._allele_key = r._Object_key
    and cv._refs_key = r._refs_key
    and r._mgitype_key = 11
    and aa._allele_type_key = t._term_key
    and aa.iswildtype = 0
    and not exists (select 1 from alleleExclude e where aa._allele_key = e._allele_key)
    and aa.symbol in ('Tg(K18-ACE2)2Prlmn', 'Ccl2<tm1Rol>')
    '''
    db.sql(cmd, None)
    db.sql('create index alleleKey2 on alleleSet(_allele_key)', None)
    db.sql('create index refsKey2 on alleleSet(_refs_key)', None)

    # genotype not in genotypeExclude
    # genotype has MP or DO annotation for covidRefs
    allelesByGenotype = {}

    cmd = '''
    select distinct aa.*, g._genotype_key, ga.accid as genotypeid
    into temp table genotypeSet
    from alleleSet aa, gxd_allelegenotype ag, gxd_genotype g, acc_accession ga
    where aa._allele_key = ag._allele_key
    and ag._genotype_key = g._genotype_key
    and g._genotype_key = ga._object_key
    and ga._mgitype_key = 12
    and ga._logicaldb_key = 1
    and ga.preferred = 1
    and not exists (select 1 from genotypeExclude e where g._genotype_key = e._genotype_key)
    and (exists (select 1 from voc_annot v, voc_evidence e
            where g._genotype_key = v._object_key
            and v._annottype_key = 1002
            and v._annot_key = e._annot_key
            and e._refs_key = aa._refs_key
            )
        or exists (select 1 from voc_annot v, voc_evidence e
            where g._genotype_key = v._object_key
            and v._annottype_key = 1020
            and v._annot_key = e._annot_key
            and e._refs_key = aa._refs_key
            )
        )
    '''
    db.sql(cmd, None)
    db.sql('create index genotypeKey3 on genotypeSet(_genotype_key)', None)
    db.sql('create index alleleKey3 on genotypeSet(_allele_key)', None)
    db.sql('create index refsKey3 on genotypeSet(_refs_key)', None)
    results = db.sql('select * from genotypeSet', 'auto')
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
    select distinct aa._allele_key, t.term
    from alleleSet aa, voc_annot v, voc_term t
    where aa._allele_key = v._object_key
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

    # mpGenotypes with mp annotations with covid reference
    #
    cmd = '''
    select s._allele_key, s._refs_key, g._genotype_key, a.accid as genotypeid, 
        t.label, t._object_key, tt.term as otherLabel, mp.accid as mpid
    from genotypeSet s, gxd_genotype g, acc_accession a,
            voc_annotheader h, mgi_setmember t, acc_accession mp, voc_term tt
    where s._genotype_key = g._genotype_key
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
    order by s._allele_key, s._refs_key, genotypeid, t.label
    '''
    results = db.sql(cmd, 'auto')
    for r in results:
        if r['label'] == None:
            r['label'] = r['otherLabel']
        key = r['_allele_key'] + r['_refs_key']
        value = r
        if key not in mpGenotypes:
            mpGenotypes[key] = []
        mpGenotypes[key].append(value)
    #print(mpGenotypes)

    #
    # doGenotypes with DO annotations with at least 1 covid reference
    #
    cmd = '''
    select s._allele_key, s._refs_key, g._genotype_key, a.accid as genotypeid, 
        tt.term, aa.accid as doid
    from genotypeSet s, gxd_genotype g, acc_accession a, acc_accession aa, voc_annot vt, voc_term tt
    where s._genotype_key = g._genotype_key
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
    order by s._allele_key, s._refs_key, genotypeid, tt.term
    '''
    results = db.sql(cmd, 'auto')
    for r in results:
        key = r['_allele_key'] + r['_refs_key']
        value = r
        if key not in doGenotypes:
            doGenotypes[key] = []
        doGenotypes[key].append(value)
    #print(doGenotypes)

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
def processMP(mode, fp, r, key):

    if mode == 'strain':
        strain = r['strain']
        strainKey = r['_strain_key']

    refsKey = r['_refs_key']

    for g in mpGenotypes[key]:

        gKey = g['_genotype_key']

        # no alleles
        if mode == 'strain' and gKey not in allelesByGenotype:

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

            if mode == 'strain':

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

            else:

                fp.write(a['symbol'] + TAB)
                fp.write(a['alleleid'] + TAB)
                fp.write(a['alleletype'] + TAB)
                fp.write(g['genotypeid'] + TAB)

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
def processDO(mode, fp, r, key):

    if mode == 'strain':
        strain = r['strain']
        strainKey = r['_strain_key']

    refsKey = r['_refs_key']

    for g in doGenotypes[key]:

        gKey = g['_genotype_key']

        # no alleles
        if mode == 'strain' and gKey not in allelesByGenotype:
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

            if mode == 'strain':

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

            else:

                fp.write(a['symbol'] + TAB)
                fp.write(a['alleleid'] + TAB)
                fp.write(a['alleletype'] + TAB)
                fp.write(g['genotypeid'] + TAB)

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
def processByStrain():

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
            fpStrain.write(r['strain'] + TAB)
            fpStrain.write(r['strainid'] + TAB)

            if strainKey in strainAttrs:
                fpStrain.write('|'.join(strainAttrs[strainKey]) + TAB)
            else:
                fpStrain.write(TAB)

            fpStrain.write(TAB*9)
            fpStrain.write(r['short_citation'] + TAB)
            fpStrain.write(r['jnumid'] + TAB)
            fpStrain.write('|'.join(covidTags[refsKey]) + CRT)
            continue

        # MP annotations
        if key in mpGenotypes:
           processMP('strain', fpStrain, r, key)

        # DO annotations
        if key in doGenotypes:
           processDO('strain', fpStrain, r, key)

#
# process strain by allele-reference associations
#
def processByAllele():

    #
    # alleleSet
    #
    results = db.sql('select * from alleleSet order by jnumid, symbol, short_citation', 'auto')
    for r in results:

        alleleKey = r['_allele_key']
        refsKey = r['_refs_key']
        key = alleleKey + refsKey

        # genotypes but no mp or do annotations
        if key not in mpGenotypes and key not in doGenotypes:
            fpAllele.write(r['symbol'] + TAB)
            fpAllele.write(r['alleleid'] + TAB)
            fpAllele.write(r['alleletype'] + TAB)
            fpAllele.write(TAB)

            if alleleKey in alleleSubtypes:
                fpAllele.write('|'.join(alleleSubtypes[alleleKey]) + TAB)
            else:
                fpAllele.write(TAB)

            fpAllele.write(TAB*4)
            fpAllele.write(r['short_citation'] + TAB)
            fpAllele.write(r['jnumid'] + TAB)
            fpAllele.write('|'.join(covidTags[refsKey]) + CRT)

            continue

        # MP annotations
        if key in mpGenotypes:
           processMP('allele', fpAllele, r, key)

        # DO annotations
        if key in doGenotypes:
           processDO('allele', fpAllele, r, key)

#
# main
#

db.useOneConnection(1)

fpStrain = reportlib.init('MGI_Cov_Strain', printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])
fpAllele = reportlib.init('MGI_Cov_Allele', printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])

fpStrain.write('#\n')
fpStrain.write('#1/A) Strain name\n')
fpStrain.write('#2/B) Strain MGI ID \n')
fpStrain.write('#3/C) Strain attributes (pipe delimited) \n')
fpStrain.write('#4/D) Genotype MGI ID\n')
fpStrain.write('#5/E) Allele symbol\n')
fpStrain.write('#6/F) Allele MGI ID \n')
fpStrain.write('#7/G) Allele type \n')
fpStrain.write('#8/H) Allele subtypes (pipe delimited) \n')
fpStrain.write('#9/I) DO term\n')
fpStrain.write('#10/J) DO ID \n')
fpStrain.write('#11/K) MP term\n')
fpStrain.write('#12/L) MP ID\n')
fpStrain.write('#13/M) Reference short citation \n')
fpStrain.write('#14/N) Reference JnumID \n')
fpStrain.write('#15/O) Reference-associated COV tags (pipe delimited) \n')
fpStrain.write('#\n')

fpAllele.write('#\n')
fpAllele.write('#1/A) Allele symbol\n')
fpAllele.write('#2/B) Allele MGI ID \n')
fpAllele.write('#3/C) Allele type \n')
fpAllele.write('#4/D) Genotype MGI ID\n')
fpAllele.write('#5/E) Allele subtypes (pipe delimited) \n')
fpAllele.write('#6/F) DO term\n')
fpAllele.write('#7/G) DO ID \n')
fpAllele.write('#8/H) MP term\n')
fpAllele.write('#9/I) MP ID\n')
fpAllele.write('#10/J) Reference short citation \n')
fpAllele.write('#11/K) Reference JnumID \n')
fpAllele.write('#12/L) Reference-associated COV tags (pipe delimited) \n')
fpAllele.write('#\n')

initializeRefSet()
initializeAlleleExclude()

#initializeByStrain()
#processByStrain()

initializeByAllele()
processByAllele()

reportlib.finish_nonps(fpStrain)
reportlib.finish_nonps(fpAllele)
db.useOneConnection(0)

