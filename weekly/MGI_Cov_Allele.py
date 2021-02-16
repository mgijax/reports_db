'''
#
# MGI_Cov_Allele.py
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
# alleles covid-reference associations
#     allele is not wild-type
#     allele not in alleleExclude
#
# for alleles that contain genotype (gxd_allelegenotype)
#     genotype with 1 marker
#     genotype must contain MP or DO annotation using covid-reference
#
# Output format:
#
# 1/A) Allele symbol
# 2/B) Allele MGI ID 
# 3/C) Allele type 
# 4/D) Allele subtypes  (pipe delimited)
# 5/E) Genotype MGI ID
# 6/F) DO term
# 7/G) DO ID 
# 8/H) MP term
# 9/I) MP ID
# 10/J) Reference short citation 
# 11/K) Reference JnumID 
# 12/L) Reference-associated COV tags (pipe delimited) 
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
mpGenotypes = {}
doGenotypes = {}
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
# allele by allele-reference associations
#
def initializeAllele():
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
    --and aa.symbol in ('Tg(K18-ACE2)2Prlmn')
    '''
    db.sql(cmd, None)
    db.sql('create index alleleKey2 on alleleSet(_allele_key)', None)
    db.sql('create index refsKey2 on alleleSet(_refs_key)', None)

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
    select distinct s._allele_key, s._refs_key, g._genotype_key, a.accid as genotypeid, 
        t.label, t._object_key, tt.term as otherLabel, mp.accid as mpid
    from alleleSet s, gxd_allelegenotype g, acc_accession a,
            voc_annotheader h, mgi_setmember t, acc_accession mp, voc_term tt, voc_annot v, voc_evidence e
    where s._allele_key = g._allele_key
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
    and h._object_key = v._object_key
    and v._annottype_key = 1002
    and v._annot_key = e._annot_key
    and e._refs_key = s._refs_key
    and not exists (select 1 from genotypeExclude e where g._genotype_key = e._genotype_key)
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
    # doGenotypes with DO annotations with covid reference
    #
    cmd = '''
    select distinct s._allele_key, s._refs_key, g._genotype_key, a.accid as genotypeid, tt.term, aa.accid as doid
    from alleleSet s, gxd_allelegenotype g, acc_accession a, acc_accession aa, voc_annot vt, voc_term tt, voc_annot v, voc_evidence e
    where s._allele_key = g._allele_key
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
    and g._genotype_key = v._object_key
    and v._annottype_key = 1002
    and v._annot_key = e._annot_key
    and e._refs_key = s._refs_key
    and not exists (select 1 from genotypeExclude e where g._genotype_key = e._genotype_key)
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
# mp annotations
#
def processMP(r, key):

    alleleKey = r['_allele_key']
    refsKey = r['_refs_key']

    for g in mpGenotypes[key]:

        fp.write(r['symbol'] + TAB)
        fp.write(r['alleleid'] + TAB)
        fp.write(r['alleletype'] + TAB)

        if alleleKey in alleleSubtypes:
            fp.write('|'.join(alleleSubtypes[alleleKey]) + TAB)
        else:
            fp.write(TAB)

        fp.write(g['genotypeid'] + TAB)
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

    alleleKey = r['_allele_key']
    refsKey = r['_refs_key']

    for g in doGenotypes[key]:

        fp.write(r['symbol'] + TAB)
        fp.write(r['alleleid'] + TAB)
        fp.write(r['alleletype'] + TAB)

        if alleleKey in alleleSubtypes:
            fp.write('|'.join(alleleSubtypes[alleleKey]) + TAB)
        else:
            fp.write(TAB)

        fp.write(g['genotypeid'] + TAB)
        fp.write(g['term'] + TAB)
        fp.write(g['doid'] + TAB)
        fp.write(TAB*2)
        fp.write(r['short_citation'] + TAB)
        fp.write(r['jnumid'] + TAB)
        fp.write('|'.join(covidTags[refsKey]) + CRT)

#
# process allele-reference associations
#
def processAllele():

    #
    # alleleSet
    #
    results = db.sql('select * from alleleSet order by symbol, jnumid, short_citation', 'auto')
    for r in results:

        alleleKey = r['_allele_key']
        refsKey = r['_refs_key']
        key = alleleKey + refsKey

        # allele without genotype/mp or genotype/do annotations
        if key not in mpGenotypes and key not in doGenotypes:
            fp.write(r['symbol'] + TAB)
            fp.write(r['alleleid'] + TAB)
            fp.write(r['alleletype'] + TAB)

            if alleleKey in alleleSubtypes:
                fp.write('|'.join(alleleSubtypes[alleleKey]) + TAB)
            else:
                fp.write(TAB)

            fp.write(TAB*5)
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
fp.write('# This report contains mouse alleles used in coronavirus research references in MGI.\n')
fp.write('# When available Human Disease Ontology (DO) and system-level Mammalian Phenotype ontology (MP) annotations for the allele are included.\n')
fp.write('# A row consists of an allele-reference-DO/MP term triplet.\n')
fp.write('# If no MP or DO terms are associated with an allele-reference pair then just the strain and reference are present in the row.\n')
fp.write('#\n')
fp.write('#1/A) Allele symbol\n')
fp.write('#2/B) Allele MGI ID \n')
fp.write('#3/C) Allele type \n')
fp.write('#4/D) Allele subtypes (pipe delimited) \n')
fp.write('#5/E) Genotype MGI ID\n')
fp.write('#6/F) DO term\n')
fp.write('#7/G) DO ID \n')
fp.write('#8/H) MP term\n')
fp.write('#9/I) MP ID\n')
fp.write('#10/J) Reference short citation \n')
fp.write('#11/K) Reference JnumID \n')
fp.write('#12/L) Reference-associated COV tags (pipe delimited) \n')
fp.write('#\n')

initializeRefSet()
initializeAlleleExclude()
initializeAllele()
processAllele()

reportlib.finish_nonps(fp)
db.useOneConnection(0)

