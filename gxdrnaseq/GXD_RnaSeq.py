'''
#
# GXD_RnaSeq.py
#
# Generate 1 report per experiment
#
# 1:  MGI Gene ID
# 2:  Ensembl ID
# 3:  Gene Symbol
# 4:  Gene Name
# 5:  Experiment ID
# 6:  Anatomical Structure [gxd_htsample._emapa_key]
# 7:  Theiler Stage [gxd_htsample._stage_key]
# 8:  Age [gxd_htsample.age]
# 9:  Sex [gxd_htsample._sex_key]
# 10: Strain [gxd_htsample._genotype_key -> gxd_genotype -> prb_strain]
# 11: Mutant Allele Pair(s) (comma delimited, if multiple)
# 12: Notes (RNA-Seq) [htsample note]
# 13: Sample ID (maybe combine the bioreplicate samples - comma delimited)
# 14: Number of Biological Replicates
# 15: Bioreplicate Set Label : as used in the heat map
#   emapa term + "_" + experiiment ID + "_" + gxd_htsample_rnaseqset._rnaseqset_key
#   example: hippocampus_E-GEOD-76567_323)
# 16: Detected : if TPM Level == 'Below Cutoff' then No, else Yes
# 17: avg_TPM [gxd_htsample_rnaseq..averagetpm]
# 18: qnTPM [gxd_htsample_rnaseq..quantilenormalizedtpm]
# 19: avg_qnTPM [gxd_htsample_rnaseqcombined..averagequantilenormalizedtpm]
# 20: TPM Level
#
# Sort
# 1) Gene Symbol
# 2) Experiment ID
# 3) Anatomical Structure
# 4) Theiler Stage
# 5) Age
# 6) Sex
# 7) Strain
# 8) Mutant Allele Pair(s) (if multiple, comma delimited and alpha sorted)
# 9) Sample Note
# 
# History:
#
# lec	09/11/2024
#	- wts2-495/e4g2/Public Report(s) of Comprehensive GXD RNA-Seq Data
#
'''
 
import sys 
import os
import reportlib
import db

db.setTrace()

CRT = reportlib.CRT
#TAB = reportlib.TAB
TAB = "|"

# distinct experiments
db.sql('''
select distinct s._experiment_key, a.accid as exptId
into temp table experiments
from GXD_HTSample s, ACC_Accession a
where exists (select 1 from GXD_HTSample_RNASeqSetMember rna where s._sample_key = rna._sample_key)
and s._experiment_key = a._object_key
and a._mgitype_key = 42
and a._logicaldb_key in (189)
--and a.accid in ('E-MTAB-9192','E-MTAB-7279')
and a.accid in ('E-GEOD-22131')
order by a.accid
''', None)
db.sql('create index eidx1 on experiments (_experiment_key);', None)

# sampleNotes by sample key
sampleNotes = {}
results = db.sql('''
select distinct n._object_key, n.note
from experiments e, GXD_HTSample s, MGI_Note n
where e._experiment_key = s._experiment_key
and s._sample_key = n._object_key
and n._notetype_key = 1048
and exists (select 1 from GXD_HTSample_RNASeqSetMember rna where s._sample_key = rna._sample_key)
''', 'auto')
for r in results:
    key = r['_object_key']
    value = r['note']
    sampleNotes[key] = value
#print(sampleNotes)

# alleles by genotype
alleles = {}
results = db.sql('''
select distinct s._genotype_key, n.note as alleles
from experiments e, GXD_HTSample_RNASeqSet s, MGI_Note n
where e._experiment_key = s._experiment_key
and s._genotype_key = n._object_key
and n._mgitype_key = 12
and n._notetype_key = 1016
and exists (select 1 from GXD_HTSample_RNASeqSetMember rna where s._rnaseqset_key = rna._rnaseqset_key)
order by alleles
''', 'auto')
for r in results:
    key = r['_genotype_key']
    value = r['alleles'].replace('\n',',')
    if value.endswith(","):
        value = value[:-1]
    alleles[key] = value
#print(alleles)

# distinct markers used in RNASeqCombined
markers = {}
db.sql('''
WITH marker AS (
select distinct m._marker_key, m.symbol, m.name
from GXD_HTSample_RNASeqCombined rna, MRK_Marker m
where rna._marker_key = m._marker_key
)
select m.*, a1.accid as mgiId, array_to_string(array_agg(distinct a2.accid),',') as ensId
into temp table markers
from marker m, ACC_Accession a1, ACC_Accession a2
where m._marker_key = a1._object_key
and a1._mgitype_key = 2
and a1._logicaldb_key = 1
and a1.preferred = 1
and m._marker_key = a2._object_key
and a2._mgitype_key = 2
and a2._logicaldb_key = 60
and a2.preferred = 1
group by 1,2,3,4
''', None)
db.sql('create index midx1 on markers (_marker_key);', None)
results = db.sql('select * from markers', 'auto')
for r in results:
    # skip marker if > 1 ensId
    if "," in r['ensId']:
        continue
    key = r['_marker_key']
    value = r
    markers[key] = value
#print(markers)

# iterate thru one experimen at a time
counter=1
eresults = db.sql('select * from experiments order by _experiment_key', 'auto')
for e in eresults:

    # create 1 report per experiment
    print(str(counter) + ':' + e['exptId'])
    counter += 1

    fp = reportlib.init(e['exptId'], outputdir = os.environ['GXDRNASEQDIR'], printHeading = None)
    fp.write('MGI Gene ID' + TAB)
    fp.write('Ensembl ID' + TAB)
    fp.write('Gene Symbol' + TAB)
    fp.write('Gene Name' + TAB)
    fp.write('Experiment ID' + TAB)
    fp.write('Anatomical Structure' + TAB)
    fp.write('Theiler Stage' + TAB)
    fp.write('Age' + TAB)
    fp.write('Sex' + TAB)
    fp.write('Strain' + TAB)
    fp.write('Mutant Allele Pair(s)' + TAB)
    fp.write('Notes' + TAB)
    fp.write('Sample ID' + TAB)
    fp.write('Number of Biological Replicates' + TAB)
    fp.write('Bioreplicate Set Label' + TAB)
    fp.write('Detected' + TAB)
    fp.write('avg_TPM' + TAB)
    fp.write('qnTPM' + TAB)
    fp.write('avg_qnTPM' + TAB)
    fp.write('TPM Level\n')

    eKey = e['_experiment_key']

    # sample info of given experiment
    sampleByExpt = {}
    results = db.sql('''
    select distinct e.exptId, s._sample_key, s.name, ss.*,
        s1.term as termStruct, s2.term as termSex, gs.strain
    from experiments e, GXD_HTSample s, GXD_HTSample_RNASeqSet ss, GXD_HTSample_RNASeqSetMember sm,
         voc_term s1, voc_term s2, gxd_genotype g, prb_strain gs
    where e._experiment_key = %s
    and e._experiment_key = s._experiment_key
    and s._sample_key = sm._sample_key
    and sm._rnaseqset_key = ss._rnaseqset_key
    and ss._emapa_key = s1._term_key
    and ss._sex_key = s2._term_key
    and ss._genotype_key = g._genotype_key
    and g._strain_key = gs._strain_key
    order by exptId, termStruct, _stage_key, age, termSex, strain
    ''' % (eKey), 'auto')
    for r in results:
        key = r['_sample_key']
        value = r
        sampleByExpt[key] = value
    #print(sampleByExpt)

    for mKey in markers:

        #print('mKey: ', mKey)

        # sample info of given experiment/marker
        results = db.sql('''
           select distinct rna._sample_key,
                rna.averagetpm, rna.quantilenormalizedtpm,
                rnaC._level_key, rnaC.numberofbiologicalreplicates, rnaC.averagequantilenormalizedtpm,
                rnaC._rnaseqset_key,
                s1.term as tpmLevel
            from experiments e, GXD_HTSample s, GXD_HTSample_RNASeq rna, GXD_HTSample_RNASeqCombined rnaC, voc_term s1
            where e._experiment_key = %s
            and e._experiment_key = s._experiment_key
            and s._sample_key = rna._sample_key
            and rna._marker_key = %s
            and rna._rnaseqcombined_key = rnaC._rnaseqcombined_key
            and rnaC._level_key = s1._term_key
            and rnaC._createdby_key = 1613
            union
            select distinct rna._sample_key,
                -1, -1,
                rnaC._level_key, rnaC.numberofbiologicalreplicates, rnaC.averagequantilenormalizedtpm,
                rnaC._rnaseqset_key,
                s1.term as tpmLevel
            from experiments e, GXD_HTSample s, GXD_HTSample_RNASeqSetMember rna, GXD_HTSample_RNASeqCombined rnaC, voc_term s1
            where e._experiment_key = %s
            and e._experiment_key = s._experiment_key
            and s._sample_key = rna._sample_key
            and rna._rnaseqset_key = rnaC._rnaseqset_key
            and rnaC._level_key = s1._term_key
            and rnaC._marker_key = %s
            and rnaC._createdby_key = 1673
        ''' % (e['_experiment_key'], mKey, e['_experiment_key'], mKey), 'auto')

        # iterate thru each experiment/sample result
        for r in results:

            sKey = r['_sample_key']
            #print('sKey: ', sKey)

            # 1:  MGI Gene ID
            # 2:  Ensembl ID
            # 3:  Gene Symbol
            # 4:  Gene Name
            # 5:  Experiment ID
            fp.write(markers[mKey]['mgiId'] + TAB)
            fp.write(markers[mKey]['ensId'] + TAB)
            fp.write(markers[mKey]['symbol'] + TAB)
            fp.write(markers[mKey]['name'] + TAB)
            fp.write(sampleByExpt[sKey]['exptId'] + TAB)

            # 6:  Anatomical Structure [gxd_htsample._emapa_key]
            # 7:  Theiler Stage [gxd_htsample._stage_key]
            # 8:  Age [gxd_htsample.age]
            # 9:  Sex [gxd_htsample._sex_key]
            # 10: Strain [gxd_htsample._genotype_key -> gxd_genotype -> prb_strain]
            fp.write(sampleByExpt[sKey]['termStruct'] + TAB)
            fp.write(str(sampleByExpt[sKey]['_stage_key']) + TAB)
            fp.write(sampleByExpt[sKey]['age'] + TAB)
            fp.write(sampleByExpt[sKey]['termSex'] + TAB)
            fp.write(sampleByExpt[sKey]['strain'] + TAB)

            # may/may not exist
            # 11: Mutant Allele Pair(s) (comma delimited, if multiple)
            gKey = sampleByExpt[sKey]['_genotype_key']
            if gKey in alleles:
                fp.write(alleles[gKey])
            fp.write(TAB)

            # 12: Notes (RNA-Seq) [htsample note _notetype_key = 1048]
            if sKey in sampleNotes:
                fp.write(sampleNotes[sKey])
            fp.write(TAB)

            # 13: Sample ID (name)
            fp.write(sampleByExpt[sKey]['name'] + TAB)

            # 14: Number of Biological Replicates
            fp.write(str(r['numberofbiologicalreplicates']) + TAB)

            # 15: Bioreplicate Set Label
            fp.write(sampleByExpt[sKey]['termStruct'] + '_' + sampleByExpt[sKey]['exptId'] + '_' + str(r['_rnaseqset_key']) + TAB)

            # 16: Detected
            if r['tpmLevel'] == 'Below Cutoff':
                fp.write('No' + TAB)
            else:
                fp.write('Yes' + TAB)

            # 17: avg_TPM
            # 18: qnTPM
            # 19: avg_qnTPM
            # 20: TPM Level

            if r['averagetpm'] >= 0:
                fp.write(str(r['averagetpm']))
            fp.write(TAB)

            if r['quantilenormalizedtpm'] >= 0:
                fp.write(str(r['quantilenormalizedtpm']))
            fp.write(TAB)

            fp.write(str(r['averagequantilenormalizedtpm']) + TAB)
            fp.write(r['tpmLevel'] + CRT)
    
    reportlib.finish_nonps(fp)	# non-postscript file
    sys.stdout.flush()

# end for e in eresults:

