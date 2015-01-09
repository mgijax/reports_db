#!/usr/local/bin/python

'''
#
# MGI_iphone_app.py
#
# Report:
#       Tab-delimited file for iphone app
#	counts and accession ids by mouse marker
#	see TR11046 for more details
#
# Used by:
#	Jill Recla (jillrecla@jax.org), Carol Bult (carol.bult@jax.org)
#
# Output format:
#
# iphone-genes
#
#	Markers of all types
#	Markers of status official or interum (_marker_status_key in (1,3)
#
#	1:  MGI Marker ID
#	2:  Marker key
#	3:  Symbol
#	4:  Name
#
#	5: MGI Ref ID (MGI:xxx||MGI:xxx...)
#
#	6: MGI Allele ID (MGI:xxx|MGI:xxx|...)
#	7: GO ID (C group): (GO:xxxx||GO:xxxx||...)
#
#	8: GO ID (F group): (GO:xxxx||GO:xxxx||...)
#
#	9: GO ID (P group): (GO:xxxx||GO:xxxx||...)
#
#	10: MP ID: (MP:xxxx||MP:xxxx||...)
#
#	11: OMIM ID: (xxxx||xxxx||...)
#
#	12: OMIM ID: (xxxx||xxxx||...)
#
#####
#TR11732 Changes
#Remove the following 23 columns:
#       5:  URL prefix (by marker id)
#       6.  Genetic position (ex. 'ChrX 36.02 cM')
#       7.  Genomic position (ex. 'ChrX:67639052-67642592 (+)')
#       8:  # of References
#       10: URL prefix (by reference (J:))
#       11:  # of Alleles
#       13: URL prefix (by allele id)
#       14: # of GO annotations (C group)
#       15: GO ID (C group): (GO:xxxx|term||GO:xxxx|term||...)
#       16: URL prefix (by GO id)
#       17: # of GO annotations (F group)
#       19: URL prefix (by GO id)
#       20: # of GO annotations (P group)
#       22: URL prefix (by GO id)
#       23: # of MP annotations
#       25: URL prefix (by MP id)
#       26: # of OMIM annotations (via genotype annotations)
#       28: URL prefix (by OMIM id)
#       29: # of OMIM annotations (via human disease via mouse ortholog)
#       31: URL prefix (by OMIM id)
#       32: # of Nomenclature Events (via marker history)
#       33: Nomenclature sequence number and event:  
#       34: URL prefix (by marker)
#
#Remove the terms from the following columns (leave only IDs):
#       15 - GO IDs (cellular component group)
#       18 - GO IDs (molecular function group)
#       21 - GO IDs (biological process group)
#       24 - MP IDs
#       27 - OMIM IDs
#       30 - OMIM IDs
##################################
# iphone-mp
#
#	MP annotations
#
#	1: MP ID
#	2: MP Term
#	3: MP Definition
#
#	4: MGI Ref ID (MGI:xxx||MGI:xxx||...)
#
#	5: MGI Genotype ID (MGI:xxx|MGI:xxx|...)
#
#       6: MGI Marker ID (MGI:xxx|MGI:xxx|...)
#
#	7: MGI Allele ID (MGI:xxx|MGI:xxx|...)
######
#TR11732 Changes
#Remove the following 8 columns:
#       4 - # of references
#       6 - url prefix
#       7 - # of genotypes
#       9 - url prefix
#       10 - # of markers
#       12 - url prefix
#       13 - # of alleles
#       15 - url prefix
#
#Remove the J#s from the following columns (leave only IDs):
#       5 - MGI IDs of annotated references
#
###################################
# iphone-omim
#
#	1: OMIM ID
#	2: OMIM Definition
#	3: URL prefix (by omim id: omim page)
#
#	for (_annottype_key = 1005)
#	4: # of References
#	5: MGI Ref ID (MGI:xxx|J:xxx||MGI:xxx|J:xxx||...)
#	6: URL prefix (by reference (J:))
#
#	7: # of Genotypes
#	8: MGI Genotype ID (MGI:xxx|MGI:xxx|...)
#	9: URL prefix (by genotype)
#
#	10: # of Markers (_AnnotType_key = 1016)
#	11: MGI Marker ID (MGI:xxx|MGI:xxx|...)
#	12: URL prefix (by marker)
#
#	13: # of Alleles
#	14: MGI Allele ID (MGI:xxx|MGI:xxx|...)
#	15: URL prefix (by allele)
#
#	for (_annottype_key = 1006)
#	16: # of References
#	17: MGI Ref ID (MGI:xxx|J:xxx||MGI:xxx|J:xxx||...)
#	18: URL prefix (by reference (J:))
#
#	19: # of Human Gene Annotations
#	20: Human EntrezGene ID (xxx|xxx|...)
#	21: URL prefix (by omim id: human disease page)
#
######
#TR11732 Changes
#Remove the following 13 columns:
#       3 - url prefix
#       4 - # of references
#       6 - url prefix
#       7 - # of genotypes
#       9 - url prefix
#       10 - # of markers
#       12 - url prefix
#       13 - # of alleles
#       15 - url prefix
#       16 - # of references
#       18 - url prefix
#       19 - # of human gene annotations
#       21 - url prefix
#
#Remove the J#s from the following columns (leave only IDs):
#       5 - MGI IDs of annotated references
#       17 - MGI IDs of annotated references
#
####
##TR11732 Changes
# replaced the double pipes separator, in each report, with single pipe  
# and got rid of the trailing double pipes 
#####################################
#
# History:
#
# lnh  08/07/2014
#      - TR11732 : trim down iphone app reports
#
# lnh  07/23/2014
#      - TR11725 : store previous week's iPhone app reports
#
# sc	04/18/2013 
#	- N2MO; update to use new MRK_Cluster* tables
#
# lec	05/17/2012
#	- TR11046; new
#
'''

import sys
import os
import subprocess
import string
import mgi_utils
import reportlib
import pg_db
db = pg_db
db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE()

TAB = reportlib.TAB
CRT = reportlib.CRT

# constants

MP_MARKER_ANNOT_TYPE = 1015
OMIM_MARKER_ANNOT_TYPE = 1016

NOT_QUALIFIER = 1614157
NORMAL_QUALIFIER = 2181424

#
# URL prefxies
#
# markers & references:  use reportlib.create_accession_anchor()
#
allele_urlprefix = 'http://www.informatics.jax.org/allele/'
go_urlprefix = 'http://www.informatics.jax.org/searches/GO.cgi?id='
mp_urlprefix = 'http://www.informatics.jax.org/searches/Phat.cgi?id='
geno_urlprefix = 'http://www.informatics.jax.org/javawi2/servlet/WIFetch?page=mpAnnotSummary&id='
omim_urlprefix = 'http://www.omim.org/entry/'
omimdisease_urlprefix = 'http://www.informatics.jax.org/disease/'

#
# report directory, archive, etc.
#
reportDir = os.environ['REPORTOUTPUTDIR']
archiveDir = os.environ['IPHONEARCHIVE']
currentDate = mgi_utils.date('%m-%d-%Y')

#
# little ditty to remove the '<A HREF="' and '">' from the accession anchor
# returned by reportlib.create_accession_anchor()
#
def print_url_prefix(url):

    url = url.replace('<A HREF="', '')
    url = url.replace('">', '')
    return url

#
# report lives in repoort directory (not the FTP directory)
# archive lives in the FTP directory (reports/archive/...)
#
def init_report(reportName):

    reportWithDate = '%s-%s' % (reportName, currentDate)
    currentReport = '%s-current.rpt' % (reportName)
    #
    # addons to support TR 11725 start here
    #
    oldReportLink='%s-old.rpt' % (reportName)
    tempDir='%s/iphone-temp' % (reportDir)
    currentReportFile='%s/%s-current_file' % (tempDir,reportName)
    currentReportName=""
    currentDir=os.getcwd()
    if os.path.isdir('%s' % (tempDir) ):
       os.system('rm -rf %s' % (tempDir))
    if os.path.isfile('%s/%s' % (reportDir, currentReport)):
       #create a temp directory to temporary store current report info 
       os.system('mkdir %s' % (tempDir))
       os.system('readlink %s/%s > %s'% (reportDir,currentReport,currentReportFile))
       if os.path.isfile('%s' % (currentReportFile)):
          f = open(currentReportFile)
          lines = f.readlines()     
          f.close()
          for report_name in lines:
              report_name= report_name.strip()
              if report_name:
                 os.system('cp -p %s/%s %s/' % (reportDir,report_name,tempDir))
                 currentReportName='%s' % (report_name)
       #Remove old report symbolic link if exists
       if os.path.isfile('%s/%s' % (reportDir,oldReportLink)):
          os.remove('%s/%s' % (reportDir, oldReportLink))
       # remove current report link if it exists
       if os.path.isfile('%s/%s' % (reportDir, currentReport)):
          os.remove('%s/%s' % (reportDir, currentReport))
    # move existing reports to the archive
    os.system('mv -f %s/%s*.rpt %s' % (reportDir, reportName, archiveDir))
    
    #restore old report - and create the "old" symbolic link to point to previous week report 
    if os.path.isfile('%s/%s' % (tempDir,currentReportName)):
       os.system('mv %s/%s %s'%(tempDir,currentReportName,reportDir))     
       os.chdir(reportDir)
       os.symlink(currentReportName, oldReportLink)
       os.chdir(currentDir)

    if os.path.isdir('%s' % (tempDir) ):
       os.system('rm -rf %s' % (tempDir)) 

    fp = reportlib.init(reportWithDate, outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

    return fp, reportWithDate, currentReport

#
# iphone_genes
#
def iphone_genes():

    fp, reportWithDate, currentReport = init_report('iphone_app_snapshot_genes')

    #
    # for testing
    #	and m.chromosome = '5'
    #	and m.symbol = 'Kit'
    #	and m._marker_type_key = 1
    #
    db.sql('''
	    select m._marker_key, m.symbol, m.name, m.chromosome, o.offset, a.accid, a.numericPart
	    into #markers
	    from MRK_Marker m, ACC_Accession a, MRK_Offset o
	    where m._organism_key = 1
	    and m._marker_status_key in (1,3)
	    and m._marker_key = o._marker_key
	    and o.source = 0
	    and m._marker_key = a._object_key
	    and a._mgitype_key = 2
	    and a._logicaldb_key = 1
	    and a.prefixpart = 'MGI:'
	    and a.preferred = 1
	    ''', None)
    db.sql('create index marker_idx on #markers(_marker_key)', None)
    #
    # refs
    #
    results = db.sql('''
            select distinct m._marker_key, r.mgiid
            from #markers m, MRK_Reference r
            where m._marker_key = r._marker_key
            ''', 'auto')
    refs = {}
    for r in results:
        key = r['_marker_key']
        value = r
        if not refs.has_key(key):
            refs[key] = []
        refs[key].append(value)

    #
    # alleles
    #
    results = db.sql('''
            select distinct m._marker_key, a.accid
            from #markers m, ALL_Marker_Assoc aa, ACC_Accession a
            where m._marker_key = aa._marker_key
	    and aa._allele_key = a._object_key
	    and a._mgitype_key = 11
	    and a._logicaldb_key = 1
	    and a.prefixpart = 'MGI:'
	    and a.preferred = 1
            ''', 'auto')
    alleles = {}
    for r in results:
        key = r['_marker_key']
        value = r['accid']
        if not alleles.has_key(key):
	    alleles[key] = []
        alleles[key].append(value)

    #
    # GO annotations
    #
    results = db.sql('''
            select distinct m._marker_key, a.accid,d.dag
            from #markers m, VOC_Annot aa, ACC_Accession a,VOC_term t, DAG_Node_View d
            where m._marker_key = aa._object_key
	    and aa._annottype_key = 1000
	    and aa._term_key = a._object_key
	    and a._mgitype_key = 13
	    and a.preferred = 1
	    and aa._Term_key = t._Term_key
	    and t._Term_key = d._Object_key 
	    and t._Vocab_key = d._Vocab_key 
            ''', 'auto')
    goCannots = {}
    goFannots = {}
    goPannots = {}
    for r in results:
        key = r['_marker_key']
        value = r
    
        if r['dag'] == 'Cellular Component':
	    goannots = goCannots
        elif r['dag'] == 'Molecular Function':
	    goannots = goFannots
        elif r['dag'] == 'Biological Process':
	    goannots = goPannots

        if not goannots.has_key(key):
	    goannots[key] = []
        goannots[key].append(value)

    #
    # Phenotype Annotations (now using MP/Marker annotations derived by the
    # rollupload product)
    #
    results = db.sql('''
            select distinct m._marker_key, a.accid
            from #markers m, VOC_Annot aa, ACC_Accession a
            where m._marker_key = aa._object_key
	    and aa._annottype_key = %d
	    and aa._term_key = a._object_key
	    and a._mgitype_key = 13
	    and a.preferred = 1
	    and aa._Qualifier_key != %d
            ''' % (MP_MARKER_ANNOT_TYPE, NORMAL_QUALIFIER), 'auto')
    phenoannots = {}
    for r in results:
        key = r['_marker_key']
        value = r
    
        if not phenoannots.has_key(key):
	    phenoannots[key] = []
        phenoannots[key].append(value)

    #
    # OMIM annotations rolled up to markers (now using OMIM/Marker annotations
    # derived by the rollupload)
    #
    results = db.sql('''
	    select distinct m._marker_key, a.accid
	    from #markers m, VOC_Annot aa, ACC_Accession a
	    where m._marker_key = aa._object_key
	    and aa._annottype_key = %d
	    and aa._term_key = a._object_key
	    and a._mgitype_key = 13
	    and a.preferred = 1
	    and aa._Qualifier_key != %d
            ''' % (OMIM_MARKER_ANNOT_TYPE, NOT_QUALIFIER), 'auto')
    omimgenotype = {}
    for r in results:
        key = r['_marker_key']
        value = r
    
        if not omimgenotype.has_key(key):
	    omimgenotype[key] = []
        omimgenotype[key].append(value)

    #
    # OMIM human disease (_annottype_key = 1006)
    #
    db.sql('''select c.clusterID, cm.*, m._Organism_key
        into #mouse
        from MRK_Cluster c, MRK_ClusterMember cm, MRK_Marker m
        where c._ClusterType_key = 9272150
        and c._ClusterSource_key = 9272151
        and c._Cluster_key = cm._Cluster_key
        and cm._Marker_key = m._Marker_key
        and m._Organism_key = 1''', None)

    db.sql('create index mouse_idx on #mouse(clusterID)', None)

    db.sql('''select c.clusterID, cm.*, m._Organism_key
        into #human
        from MRK_Cluster c, MRK_ClusterMember cm, MRK_Marker m
        where c._ClusterType_key = 9272150
        and c._ClusterSource_key = 9272151
        and c._Cluster_key = cm._Cluster_key
        and cm._Marker_key = m._Marker_key
        and m._Organism_key = 2''', None)

    db.sql('create index human_idx on #human(clusterID)', None)

    results = db.sql('''select distinct m._marker_key, a.accid
            from #markers m, #mouse hm, #human hh,
	    VOC_Annot aa, ACC_Accession a
            where m._marker_key = hm._marker_key
            and hm.clusterID = hh.clusterID
            and hh._Marker_key = aa._object_key
            and aa._annottype_key = 1006
            and aa._term_key = a._object_key
            and a._mgitype_key = 13
            and a.preferred = 1''', 'auto')
    omimhuman = {}
    for r in results:
        key = r['_marker_key']
        value = r
    
        if not omimhuman.has_key(key):
	    omimhuman[key] = []
        omimhuman[key].append(value)

    #
    # report
    #
    results = db.sql('select * from #markers order by numericPart', 'auto')
    
    for r in results:
    
        key = r['_marker_key']
    #	MGI Marker ID
    #	Marker key
    #	Symbol
    #	Name
        fp.write(r['accid'] + TAB)
        fp.write(str(r['_marker_key']) + TAB)
        fp.write(r['symbol'] + TAB)
        fp.write(r['name'] + TAB)
    
    #	MGI Ref ID (MGI:xxx|MGI:xxx|...)
        if refs.has_key(key):
            i=0
	    for n in refs[key]:
                if i>0:
	           fp.write('|'+str(n['mgiid']))
                else:
                   fp.write(str(n['mgiid']))
                   i=1;
            fp.write(TAB)
        else:
            fp.write('0' + TAB)
    
    #	MGI Allele ID (MGI:xxx|MGI:xxx|...)
        if alleles.has_key(key):
	    fp.write(string.join(alleles[key], '|') + TAB)
        else:
            fp.write('0' + TAB)
    #	GO ID (C group): (GO:xxxx|GO:xxxx|...)
        if goCannots.has_key(key):
            i=0
	    for n in goCannots[key]:
                if i>0:
	           fp.write('|'+str(n['accid']))
                else:
	           fp.write(str(n['accid']))
                   i=1
	    fp.write(TAB)
        else:
            fp.write('0' + TAB)
    #	GO ID (F group): (GO:xxxx|GO:xxxx|...)
        if goFannots.has_key(key):
            i=0
	    for n in goFannots[key]:
                if i>0:
                   fp.write('|'+str(n['accid']))
                else:
                   fp.write(str(n['accid']))
                   i=1
	    fp.write(TAB)
        else:
            fp.write('0' +TAB)
    
    #	GO ID (P group): (GO:xxxx|GO:xxxx|...)
        if goPannots.has_key(key):
            i=0
	    for n in goPannots[key]:
                if i>0:
                   fp.write('|'+str(n['accid']))
                else:
                   fp.write(str(n['accid']))
                   i=1
	    fp.write(TAB)
        else:
            fp.write('0' +TAB)
    
    #	MP ID: (MP:xxxx|MP:xxxx|...)
        if phenoannots.has_key(key):
            i=0
	    for n in phenoannots[key]:
                if i>0:
                   fp.write('|'+str(n['accid']))
                else:
                   fp.write(str(n['accid']))
                   i=1
	    fp.write(TAB)
        else:
            fp.write('0' +TAB)

    #	OMIM ID: (xxxx|xxxx|...)
        if omimgenotype.has_key(key):
            i=0
	    for n in omimgenotype[key]:
                if i>0:
                   fp.write('|'+str(n['accid']))
                else:
                   fp.write(str(n['accid']))
                   i=1
	    fp.write(TAB)
        else:
            fp.write('0' + TAB)
    
    #	OMIM ID: (xxxx|xxxx|...)
        if omimhuman.has_key(key):
            i=0
	    for n in omimhuman[key]:
                if i>0:
                   fp.write('|'+str(n['accid']))
                else:
                   fp.write(str(n['accid']))
                   i=1
	    #fp.write(TAB)
        else:
            fp.write('0')
        fp.write(CRT)
    reportlib.finish_nonps(fp)
    
    # re-create a symbolic link between the new file and the current file
    os.chdir(reportDir)
    os.symlink(reportWithDate + '.rpt', currentReport)

#
# iphone_mp
#
# MP terms -> genotype, genes, alleles
#
def iphone_mp():

    fp, reportWithDate, currentReport = init_report('iphone_app_snapshot_mp')

    #
    # MP terms
    #
    db.sql('''
            select t._term_key, t.term, a.accid
	    into #mp
            from VOC_Term t, ACC_Accession a
            where t._Vocab_key = 5
	    and t._term_key = a._object_key
	    and a._mgitype_key = 13
	    and a.preferred = 1
	    order by a.accid
            ''', None)

    db.sql('create index mp_idx on #mp(_term_key)', None)

    #
    # definitions
    #
    results = db.sql('''    
        select m._term_key, x.note
        from #mp m, VOC_Text x
        where m._term_key = x._term_key
	order by m._term_key, x.sequencenum
            ''', 'auto')
    notes = {}
    for r in results:
        key = r['_term_key']
	value = r['note']
	if notes.has_key(key):
            value = notes[key] + r['note']
        notes[key] = value

    #
    # Mammalian Phenotype/References
    #
    results = db.sql('''
            select distinct m._term_key, r.mgiid
            from #mp m, VOC_Annot aa, VOC_Evidence e, BIB_Citation_Cache r
            where m._term_key = aa._term_key
	    and aa._annottype_key = 1002
	    and aa._annot_key = e._annot_key
	    and e._refs_key = r._refs_key
            ''', 'auto')
    refs = {}
    for r in results:
        key = r['_term_key']
        value = r
        if not refs.has_key(key):
	    refs[key] = []
        refs[key].append(value)

    #
    # Mammalian Phenotype/Genotype by Genotype
    #
    results = db.sql('''
            select distinct m._term_key, a.accid
            from #mp m, VOC_Annot aa, ACC_Accession a
            where m._term_key = aa._term_key
	    and aa._annottype_key = 1002
	    and aa._object_key = a._object_key
	    and a._mgitype_key = 12
	    and a._logicaldb_key = 1
	    and a.prefixpart = 'MGI:'
	    and a.preferred = 1
            ''', 'auto')
    genoannots = {}
    for r in results:
        key = r['_term_key']
        value = r['accid']
    
        if not genoannots.has_key(key):
	    genoannots[key] = []
        genoannots[key].append(value)

    #
    # Mammalian Phenotype/Genotype by Marker (updated to use the pre-computed
    # MP/Marker annotations, computed by rollupload)
    #
    results = db.sql('''
            select distinct m._term_key, a.accid
            from #mp m, VOC_Annot aa, ACC_Accession a
            where m._term_key = aa._term_key
	    and aa._annottype_key = %d
	    and aa._object_key = a._object_key
	    and a._mgitype_key = 2
	    and a._logicaldb_key = 1
	    and a.prefixpart = 'MGI:'
	    and a.preferred = 1
	    and aa._Qualifier_key != %d
            ''' % (MP_MARKER_ANNOT_TYPE, NORMAL_QUALIFIER), 'auto')
    markerannots = {}
    for r in results:
        key = r['_term_key']
        value = r['accid']
    
        if not markerannots.has_key(key):
	    markerannots[key] = []
        markerannots[key].append(value)

    #
    # Mammalian Phenotype/Genotype by Allele
    #
    results = db.sql('''
            select distinct m._term_key, a.accid
            from #mp m, VOC_Annot aa, GXD_AlleleGenotype g, ACC_Accession a
            where m._term_key = aa._term_key
	    and aa._annottype_key = 1002
	    and aa._object_key = g._genotype_key
	    and g._allele_key = a._object_key
	    and a._mgitype_key = 11
	    and a._logicaldb_key = 1
	    and a.prefixpart = 'MGI:'
	    and a.preferred = 1
            ''', 'auto')
    alleleannots = {}
    for r in results:
        key = r['_term_key']
        value = r['accid']
        if not alleleannots.has_key(key):
	    alleleannots[key] = []
        alleleannots[key].append(value)

    #
    # report
    #
    results = db.sql('select * from #mp', 'auto')
    
    for r in results:
    
        key = r['_term_key']

	fp.write(r['accid'] + TAB)
	fp.write(r['term'] + TAB)

	if notes.has_key(key):
	    fp.write(notes[key])
        else:
            fp.write('0')
	fp.write(TAB)
    #	MGI Ref ID (MGI:xxx|MGI:xxx|...)
        if refs.has_key(key):
            i=0
	    for n in refs[key]:
                if i>0:
                   fp.write('|'+str(n['mgiid']))
                else:
                   fp.write(str(n['mgiid']))
                   i=1
            fp.write(TAB)
        else:
            fp.write('0' +TAB)
    #	MGI Genotype ID (MGI:xxx|MGI:xxx|...)
        if genoannots.has_key(key):
	    fp.write(string.join(genoannots[key], '|') + TAB)
        else:
            fp.write('0' +TAB)
    #	MGI Marker ID (MGI:xxx|MGI:xxx|...)
        if markerannots.has_key(key):
	    fp.write(string.join(markerannots[key], '|') + TAB)
        else:
            fp.write('0' +TAB)
    #	MGI Allele ID (MGI:xxx|MGI:xxx|...)
        if alleleannots.has_key(key):
	    fp.write(string.join(alleleannots[key], '|'))
        else:
            fp.write('0')
        fp.write(CRT)

    reportlib.finish_nonps(fp)
    
    # re-create a symbolic link between the new file and the current file
    os.chdir(reportDir)
    os.symlink(reportWithDate + '.rpt', currentReport)

#
# iphone_omim
#
# OMIM terms -> genotype, genes, alleles
#
def iphone_omim():

    fp, reportWithDate, currentReport = init_report('iphone_app_snapshot_omim')

    #
    # OMIM terms
    #
    db.sql('''
            select t._term_key, t.term, a.accid
	    into #omim
            from VOC_Term t, ACC_Accession a
            where t._Vocab_key = 44
	    and t._term_key = a._object_key
	    and a._mgitype_key = 13
	    and a.preferred = 1
	    order by a.accid
            ''', None)

    db.sql('create index omim_idx on #omim(_term_key)', None)

    #
    # OMIM/Genotype/References
    #
    results = db.sql('''
            select distinct m._term_key, r.mgiid
            from #omim m, VOC_Annot aa, VOC_Evidence e, BIB_Citation_Cache r
            where m._term_key = aa._term_key
	    and aa._annottype_key = 1005
	    and aa._annot_key = e._annot_key
	    and e._refs_key = r._refs_key
            ''', 'auto')
    refs1 = {}
    for r in results:
        key = r['_term_key']
        value = r
        if not refs1.has_key(key):
	    refs1[key] = []
        refs1[key].append(value)

    #
    # OMIM/Genotype by Genotype
    #
    results = db.sql('''
            select distinct m._term_key, a.accid
            from #omim m, VOC_Annot aa, ACC_Accession a
            where m._term_key = aa._term_key
	    and aa._annottype_key = 1005
	    and aa._object_key = a._object_key
	    and a._mgitype_key = 12
	    and a._logicaldb_key = 1
	    and a.prefixpart = 'MGI:'
	    and a.preferred = 1
            ''', 'auto')
    genoannots1 = {}
    for r in results:
        key = r['_term_key']
        value = r['accid']
    
        if not genoannots1.has_key(key):
	    genoannots1[key] = []
        genoannots1[key].append(value)

    #
    # OMIM/Marker pairs (now updated to use ones pre-computed by rollupload)
    #
    results = db.sql('''
            select distinct m._term_key, a.accid
            from #omim m, VOC_Annot aa, ACC_Accession a
            where m._term_key = aa._term_key
	    and aa._annottype_key = %d
	    and aa._object_key = a._object_key
	    and a._mgitype_key = 2
	    and a._logicaldb_key = 1
	    and a.prefixpart = 'MGI:'
	    and a.preferred = 1
	    and aa._Qualifier_key != %d
            ''' % (OMIM_MARKER_ANNOT_TYPE, NOT_QUALIFIER), 'auto')
    markerannots1 = {}
    for r in results:
        key = r['_term_key']
        value = r['accid']
    
        if not markerannots1.has_key(key):
	    markerannots1[key] = []
        markerannots1[key].append(value)

    #
    # OMIM/Genotype by Allele
    #
    results = db.sql('''
            select distinct m._term_key, a.accid
            from #omim m, VOC_Annot aa, GXD_AlleleGenotype g, ACC_Accession a
            where m._term_key = aa._term_key
	    and aa._annottype_key = 1005
	    and aa._object_key = g._genotype_key
	    and g._allele_key = a._object_key
	    and a._mgitype_key = 11
	    and a._logicaldb_key = 1
	    and a.prefixpart = 'MGI:'
	    and a.preferred = 1
            ''', 'auto')
    alleleannots1 = {}
    for r in results:
        key = r['_term_key']
        value = r['accid']
    
        if not alleleannots1.has_key(key):
	    alleleannots1[key] = []
        alleleannots1[key].append(value)

    #
    # OMIM/Human Marker/References
    #
    results = db.sql('''
            select distinct m._term_key, r.mgiid
            from #omim m, VOC_Annot aa, VOC_Evidence e, BIB_Citation_Cache r
            where m._term_key = aa._term_key
	    and aa._annottype_key = 1006
	    and aa._annot_key = e._annot_key
	    and e._refs_key = r._refs_key
            ''', 'auto')
    refs2 = {}
    for r in results:
        key = r['_term_key']
        value = r
        if not refs2.has_key(key):
	    refs2[key] = []
        refs2[key].append(value)

    #
    # OMIM/Human Marker by Marker
    #
    results = db.sql('''
            select distinct m._term_key, a.accid
            from #omim m, VOC_Annot aa, ACC_Accession a
            where m._term_key = aa._term_key
	    and aa._annottype_key = 1006
	    and aa._object_key = a._object_key
	    and a._mgitype_key = 2
	    and a._logicaldb_key = 55
	    and a.preferred = 1
            ''', 'auto')
    markerannots2 = {}
    for r in results:
        key = r['_term_key']
        value = r['accid']
    
        if not markerannots2.has_key(key):
	    markerannots2[key] = []
        markerannots2[key].append(value)

    #
    # report
    #
    results = db.sql('select * from #omim', 'auto')
    
    for r in results:
    
        key = r['_term_key']

	fp.write(r['accid'] + TAB)
	fp.write(r['term'] + TAB)

    #   OMIM (_annottype_key = 1006)
    #	MGI Ref ID (MGI:xxx|MGI:xxx|...)
        if refs1.has_key(key):
            i=0
	    for n in refs1[key]:
                if i>0:
                   fp.write('|'+str(n['mgiid']))
                else:
                   fp.write(str(n['mgiid']))
                   i=1
            fp.write(TAB)
        else:
            fp.write('0' + TAB)

    #   OMIM (_annottype_key = 1006)
    #	MGI Genotype ID (MGI:xxx|MGI:xxx|...)
        if genoannots1.has_key(key):
	    fp.write(string.join(genoannots1[key], '|') + TAB)
        else:
            fp.write('0' + TAB)

    #   OMIM (_annottype_key = 1006)
    #	MGI Marker ID (MGI:xxx|MGI:xxx|...)
        if markerannots1.has_key(key):
	    fp.write(string.join(markerannots1[key], '|') + TAB)
        else:
            fp.write('0' +TAB)

    #   OMIM (_annottype_key = 1006)
    #	MGI Allele ID (MGI:xxx|MGI:xxx|...)
        if alleleannots1.has_key(key):
	    fp.write(string.join(alleleannots1[key], '|') + TAB)
        else:
            fp.write('0' + TAB)

    #   OMIM (_annottype_key = 1006)
    #	MGI Ref ID (MGI:xxx|MGI:xxx|...)
        if refs2.has_key(key):
            i=0
	    for n in refs2[key]:
                if i>0:
                   fp.write('|'+str(n['mgiid']))
                else:
                   fp.write(str(n['mgiid']))
                   i=1
            fp.write(TAB)
        else:
            fp.write('0' + TAB)

    #   OMIM (_annottype_key = 1006)
    #	EntrezGene ID (xxx|xxx|...)
        if markerannots2.has_key(key):
	    fp.write(string.join(markerannots2[key], '|'))
        else:
            fp.write('0')
        fp.write(CRT)
    reportlib.finish_nonps(fp)
    
    # re-create a symbolic link between the new file and the current file
    os.chdir(reportDir)
    os.symlink(reportWithDate + '.rpt', currentReport)

#
# main
#

iphone_genes()
iphone_mp()
iphone_omim()

