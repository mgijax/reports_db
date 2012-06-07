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
#	5:  URL prefix (by marker id)
#	6:  Start Coordinate
#	7:  End Coordinate
#	8:  Strand
#
#	9:  # of References
#	10: MGI Ref ID (MGI:xxx|J:xxx||MGI:xxx|J:xxx||...)
#	11: URL prefix (by reference (J:))
#
#	12:  # of Alleles
#	13: MGI Allele ID (MGI:xxx|MGI:xxx|...)
#	14: URL prefix (by allele id)
#
#	15: # of GO annotations (C group)
#	16: GO ID (C group): (GO:xxxx|term||GO:xxxx|term||...)
#	17: URL prefix (by GO id)
#
#	18: # of GO annotations (F group)
#	19: GO ID (F group): (GO:xxxx|term||GO:xxxx|term||...)
#	20: URL prefix (by GO id)
#
#	21: # of GO annotations (P group)
#	22: GO ID (P group): (GO:xxxx|term||GO:xxxx|term||...)
#	23: URL prefix (by GO id)
#
#	24: # of MP annotations
#	25: MP ID: (MP:xxxx|term||MP:xxxx|term||...)
#	26: URL prefix (by MP id)
#
#	27: # of OMIM annotations (via genotype annotations)
#	28: OMIM ID: (xxxx|disease term||xxxx|disease term||...)
#	29: URL prefix (by OMIM id)
#
#	30: # of OMIM annotations (via human disease via mouse ortholog)
#	31: OMIM ID: (xxxx|disease term||xxxx|disease term||...)
#	32: URL prefix (by OMIM id)
#
#	33: # of Nomenclature Events (via marker history)
#	34: Nomenclature sequence number and event:  
#	35: URL prefix (by marker)
#			1|assigned||2|rename||3|split||4|deletion
# iphone-mp
#
#	MP annotations
#
#	1: MP ID
#	2: MP Term
#	3: MP Definition
#
#	4: # of References
#	5: MGI Ref ID (MGI:xxx|J:xxx||MGI:xxx|J:xxx||...)
#	6: URL prefix (by reference (J:))
#
#	7: # of Genotypes
#	8: MGI Genotype ID (MGI:xxx|MGI:xxx|...)
#	9: URL prefix (by genotype)
#
#	10: # of Markers
#	11: MGI Marker ID (MGI:xxx|MGI:xxx|...)
#	12: URL prefix (by marker)
#
#	13: # of Alleles
#	14: MGI Allele ID (MGI:xxx|MGI:xxx|...)
#	15: URL prefix (by allele)
#
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
#	10: # of Markers
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
# History:
#
# lec	05/17/2012
#	- TR11046; new
#
'''

import sys
import os
import string
import mgi_utils
import reportlib

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
        db.setAutoTranslate(False)
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db

TAB = reportlib.TAB
CRT = reportlib.CRT

#
# URL prefxies
#
# markers & references:  use reportlib.create_accession_anchor()
#
allele_urlprefix = 'http://www.informatics.jax.org/javawi2/servlet/WIFetch?page=alleleDetail&id='
go_urlprefix = 'http://www.informatics.jax.org/searches/GO.cgi?id='
mp_urlprefix = 'http://www.informatics.jax.org/searches/Phat.cgi?id='
geno_urlprefix = 'http://www.informatics.jax.org/javawi2/servlet/WIFetch?page=mpAnnotSummary&id='
omim_urlprefix = 'http://www.omim.org/entry/'
omimdisease_urlprefix = 'http://www.informatics.jax.org/javawi2/servlet/WIFetch?page=humanDisease&id='

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

    # remove current report link if it exists
    if os.path.isfile('%s/%s' % (reportDir, currentReport)):
            os.remove('%s/%s' % (reportDir, currentReport))

    # move existing reports to the archive
    os.system('mv -f %s/%s*.rpt %s' % (reportDir, reportName, archiveDir))

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
	    select m._marker_key, m.symbol, m.name, a.accid
	    into #markers
	    from MRK_Marker m, ACC_Accession a
	    where m._organism_key = 1
	    and m._marker_status_key in (1,3)
	    and m._marker_key = a._object_key
	    and a._mgitype_key = 2
	    and a._logicaldb_key = 1
	    and a.prefixpart = 'mgi:'
	    and a.preferred = 1
	    ''', None)
    db.sql('create index marker_idx on #markers(_marker_key)', None)

    #
    # coordinates
    #
    results = db.sql('''    
        select m._marker_key,
               c.strand, 
               convert(int, c.startcoordinate) as startc,
               convert(int, c.endcoordinate) as endc
        from #markers m, MRK_Location_Cache c
        where m._marker_key = c._marker_key
            ''', 'auto')
    coords = {}
    for r in results:
        key = r['_marker_key']
        value = r
        if not coords.has_key(key):
            coords[key] = []
        coords[key].append(value)

    #
    # refs
    #
    results = db.sql('''
            select distinct m._marker_key, r.mgiid, r.jnumid
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
	    and a.prefixpart = 'mgi:'
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
            select distinct m._marker_key, a.accid, t.term, d.dag
            from #markers m, VOC_Annot aa, ACC_Accession a, VOC_Term t, DAG_Node_View d
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
    # Phenotype Annotations
    #
    results = db.sql('''
            select distinct m._marker_key, a.accid, t.term
            from #markers m, VOC_Annot aa, ACC_Accession a, VOC_Term t
            where m._marker_key = aa._object_key
	    and aa._annottype_key = 1002
	    and aa._term_key = a._object_key
	    and a._mgitype_key = 13
	    and a.preferred = 1
	    and aa._Term_key = t._Term_key
            ''', 'auto')
    phenoannots = {}
    for r in results:
        key = r['_marker_key']
        value = r
    
        if not phenoannots.has_key(key):
	    phenoannots[key] = []
        phenoannots[key].append(value)

    #
    # OMIM by genotype annotations (_annottype_key = 1005)
    #
    results = db.sql('''
	    select distinct m._marker_key, a.accid, t.term
	    from #markers m, GXD_AlleleGenotype g, VOC_Annot aa, ACC_Accession a, VOC_Term t
	    where m._marker_key = g._marker_key
	    and g._genotype_key = aa._object_key
	    and aa._annottype_key = 1005
	    and aa._term_key = a._object_key
	    and a._mgitype_key = 13
	    and a.preferred = 1
	    and aa._Term_key = t._Term_key
            ''', 'auto')
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
    results = db.sql('''
	    select distinct m._marker_key, a.accid, t.term
	    from #markers m, MRK_Homology_Cache h1, MRK_Homology_Cache h2, 
	         VOC_Annot aa, ACC_Accession a, VOC_Term t
	    where m._marker_key = h1._marker_key
	    and h1._Organism_key = 1 
	    and h1._Class_key = h2._Class_key 
	    and h2._Organism_key = 2
	    and h2._Marker_key = aa._object_key
	    and aa._annottype_key = 1006
	    and aa._term_key = a._object_key
	    and a._mgitype_key = 13
	    and a.preferred = 1
	    and aa._Term_key = t._Term_key
            ''', 'auto')
    omimhuman = {}
    for r in results:
        key = r['_marker_key']
        value = r
    
        if not omimhuman.has_key(key):
	    omimhuman[key] = []
        omimhuman[key].append(value)

    #
    # Nomenclature History
    #
    results = db.sql('''
            select distinct m._marker_key, h.sequencenum, e.event
            from #markers m, MRK_History h, MRK_Event e
            where m._marker_key = h._marker_key
	    and h._marker_event_key = e._marker_event_key
            ''', 'auto')
    nomen = {}
    for r in results:
        key = r['_marker_key']
        value = r
        if not nomen.has_key(key):
	    nomen[key] = []
        nomen[key].append(value)

    #
    # report
    #
    results = db.sql('select * from #markers', 'auto')
    
    for r in results:
    
        key = r['_marker_key']
    
    #	MGI Marker ID
    #	Marker key
    #	Symbol
    #	Name
    #	URL prefix (by marker)
    
        fp.write(r['accid'] + TAB)
        fp.write(str(r['_marker_key']) + TAB)
        fp.write(r['symbol'] + TAB)
        fp.write(r['name'] + TAB)
        fp.write(print_url_prefix(reportlib.create_accession_anchor(r['accid'], 'marker')) + TAB)
    
    #	Start Coordinate
    #	End Coordinate
    #	Strand
    
        if coords.has_key(key):
            fp.write(mgi_utils.prvalue(coords[key][0]['startc']) + TAB)
            fp.write(mgi_utils.prvalue(coords[key][0]['endc']) + TAB)
            fp.write(mgi_utils.prvalue(coords[key][0]['strand']) + TAB)
        else:
            fp.write(TAB + TAB + TAB)
    
    #	# of References
    #	MGI Ref ID (MGI:xxx|MGI:xxx|...)
    #	URL prefix (by reference)
    	
        if refs.has_key(key):
	    fp.write(str(len(refs[key])) + TAB)
	    for n in refs[key]:
	        fp.write(str(n['mgiid']) + '|' + n['jnumid'] + '||')
            fp.write(TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(print_url_prefix(reportlib.create_accession_anchor('', 'reference')) + TAB)
    
    #	# of Alleles
    #	MGI Allele ID (MGI:xxx|MGI:xxx|...)
    #	URL prefix (by allele)
    
        if alleles.has_key(key):
	    fp.write(str(len(alleles[key])) + TAB)
	    fp.write(string.join(alleles[key], '|') + TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(allele_urlprefix + TAB)

    #	# of GO annotations (C group)
    #	GO ID (C group): (GO:xxxx|GO:xxxx|...)
    #	URL prefix (by marker)
    
        if goCannots.has_key(key):
	    fp.write(str(len(goCannots[key])) + TAB)
	    for n in goCannots[key]:
	        fp.write(str(n['accid']) + '|' + n['term'] + '||')
	    fp.write(TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(go_urlprefix + TAB)
    
    #	# of GO annotations (F group)
    #	GO ID (F group): (GO:xxxx|GO:xxxx|...)
    #	URL prefix (by marker)
    
        if goFannots.has_key(key):
	    fp.write(str(len(goFannots[key])) + TAB)
	    for n in goFannots[key]:
	        fp.write(str(n['accid']) + '|' + n['term'] + '||')
	    fp.write(TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(go_urlprefix + TAB)
    
    #	# of GO annotations (P group)
    #	GO ID (P group): (GO:xxxx|GO:xxxx|...)
    #	URL prefix (by marker)
    
        if goPannots.has_key(key):
	    fp.write(str(len(goPannots[key])) + TAB)
	    for n in goPannots[key]:
	        fp.write(str(n['accid']) + '|' + n['term'] + '||')
	    fp.write(TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(go_urlprefix + TAB)
    
    #	# of MP annotations
    #	MP ID: (MP:xxxx|MP:xxxx|...)
    #	URL prefix (by marker)
    
        if phenoannots.has_key(key):
	    fp.write(str(len(phenoannots[key])) + TAB)
	    for n in phenoannots[key]:
	        fp.write(str(n['accid']) + '|' + n['term'] + '||')
	    fp.write(TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(mp_urlprefix + TAB)

    #	# of OMIM annotations (via genotype)
    #	OMIM ID: (xxxx|xxxx|...)
    #	URL prefix (by marker)
    
        if omimgenotype.has_key(key):
	    fp.write(str(len(omimgenotype[key])) + TAB)
	    for n in omimgenotype[key]:
	        fp.write(str(n['accid']) + '|' + n['term'] + '||')
	    fp.write(TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(omim_urlprefix + TAB)
    
    #	# of OMIM annotations (via human disease)
    #	OMIM ID: (xxxx|xxxx|...)
    #	URL prefix (by marker)
    
        if omimhuman.has_key(key):
	    fp.write(str(len(omimhuman[key])) + TAB)
	    for n in omimhuman[key]:
	        fp.write(str(n['accid']) + '|' + n['term'] + '||')
	    fp.write(TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(omimdisease_urlprefix + TAB)

    #	# of Nomenclature Events (via marker history)
    #	Nomenclature sequence number:event:  
    #	URL prefix (by marker)
    
        if nomen.has_key(key):
	    fp.write(str(len(nomen[key])) + TAB)
	    for n in nomen[key]:
	        fp.write(str(n['sequencenum']) + '|' + n['event'] + '||')
            fp.write(TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(print_url_prefix(reportlib.create_accession_anchor(r['accid'], 'marker')) + CRT)
    
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
            select distinct m._term_key, r.mgiid, r.jnumid
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
	    and a.prefixpart = 'mgi:'
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
    # Mammalian Phenotype/Genotype by Marker
    #
    results = db.sql('''
            select distinct m._term_key, a.accid
            from #mp m, VOC_Annot aa, GXD_AlleleGenotype g, ACC_Accession a
            where m._term_key = aa._term_key
	    and aa._annottype_key = 1002
	    and aa._object_key = g._genotype_key
	    and g._marker_key = a._object_key
	    and a._mgitype_key = 2
	    and a._logicaldb_key = 1
	    and a.prefixpart = 'mgi:'
	    and a.preferred = 1
            ''', 'auto')
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
	    and a.prefixpart = 'mgi:'
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
	fp.write(TAB)

    #	# of References
    #	MGI Ref ID (MGI:xxx|MGI:xxx|...)
    #	URL prefix (by reference)
    	
        if refs.has_key(key):
	    fp.write(str(len(refs[key])) + TAB)
	    for n in refs[key]:
	        fp.write(str(n['mgiid']) + '|' + n['jnumid'] + '||')
            fp.write(TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(print_url_prefix(reportlib.create_accession_anchor('', 'reference')) + TAB)

    #	# of Genotypes
    #	MGI Genotype ID (MGI:xxx|MGI:xxx|...)
    #	URL prefix (by genotypes)
    
        if genoannots.has_key(key):
	    fp.write(str(len(genoannots[key])) + TAB)
	    fp.write(string.join(genoannots[key], '|') + TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(geno_urlprefix + TAB)

    #	# of Markers
    #	MGI Marker ID (MGI:xxx|MGI:xxx|...)
    #	URL prefix (by marker)
    
        if markerannots.has_key(key):
	    fp.write(str(len(markerannots[key])) + TAB)
	    fp.write(string.join(markerannots[key], '|') + TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(print_url_prefix(reportlib.create_accession_anchor('', 'marker')) + TAB)

    #	# of Alleles
    #	MGI Allele ID (MGI:xxx|MGI:xxx|...)
    #	URL prefix (by allele)
    
        if alleleannots.has_key(key):
	    fp.write(str(len(alleleannots[key])) + TAB)
	    fp.write(string.join(alleleannots[key], '|') + TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(allele_urlprefix + CRT)

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
            select distinct m._term_key, r.mgiid, r.jnumid
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
	    and a.prefixpart = 'mgi:'
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
    # OMIM/Genotype by Marker
    #
    results = db.sql('''
            select distinct m._term_key, a.accid
            from #omim m, VOC_Annot aa, GXD_AlleleGenotype g, ACC_Accession a
            where m._term_key = aa._term_key
	    and aa._annottype_key = 1005
	    and aa._object_key = g._genotype_key
	    and g._marker_key = a._object_key
	    and a._mgitype_key = 2
	    and a._logicaldb_key = 1
	    and a.prefixpart = 'mgi:'
	    and a.preferred = 1
            ''', 'auto')
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
	    and a.prefixpart = 'mgi:'
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
            select distinct m._term_key, r.mgiid, r.jnumid
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
	fp.write(omim_urlprefix + TAB)

    #   OMIM (_annottype_key = 1006)
    #	# of References
    #	MGI Ref ID (MGI:xxx|MGI:xxx|...)
    #	URL prefix (by reference)
    	
        if refs1.has_key(key):
	    fp.write(str(len(refs1[key])) + TAB)
	    for n in refs1[key]:
	        fp.write(str(n['mgiid']) + '|' + n['jnumid'] + '||')
            fp.write(TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(print_url_prefix(reportlib.create_accession_anchor('', 'reference')) + TAB)

    #   OMIM (_annottype_key = 1006)
    #	# of Genotypes
    #	MGI Genotype ID (MGI:xxx|MGI:xxx|...)
    #	URL prefix (by genotypes)
    
        if genoannots1.has_key(key):
	    fp.write(str(len(genoannots1[key])) + TAB)
	    fp.write(string.join(genoannots1[key], '|') + TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(geno_urlprefix + TAB)

    #   OMIM (_annottype_key = 1006)
    #	# of Markers
    #	MGI Marker ID (MGI:xxx|MGI:xxx|...)
    #	URL prefix (by marker)
    
        if markerannots1.has_key(key):
	    fp.write(str(len(markerannots1[key])) + TAB)
	    fp.write(string.join(markerannots1[key], '|') + TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(print_url_prefix(reportlib.create_accession_anchor('', 'marker')) + TAB)

    #   OMIM (_annottype_key = 1006)
    #	# of Alleles
    #	MGI Allele ID (MGI:xxx|MGI:xxx|...)
    #	URL prefix (by allele)
    
        if alleleannots1.has_key(key):
	    fp.write(str(len(alleleannots1[key])) + TAB)
	    fp.write(string.join(alleleannots1[key], '|') + TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(allele_urlprefix + TAB)

    #   OMIM (_annottype_key = 1006)
    #	# of References
    #	MGI Ref ID (MGI:xxx|MGI:xxx|...)
    #	URL prefix (by reference)
    	
        if refs2.has_key(key):
	    fp.write(str(len(refs2[key])) + TAB)
	    for n in refs2[key]:
	        fp.write(str(n['mgiid']) + '|' + n['jnumid'] + '||')
            fp.write(TAB)
        else:
            fp.write('0' + TAB + TAB)
        fp.write(print_url_prefix(reportlib.create_accession_anchor('', 'reference')) + TAB)

    #   OMIM (_annottype_key = 1006)
    #	# of Human Markers
    #	EntrezGene ID (xxx|xxx|...)
    #	URL prefix (by human/marker)
    
        if markerannots2.has_key(key):
	    fp.write(str(len(markerannots2[key])) + TAB)
	    fp.write(string.join(markerannots2[key], '|') + TAB)
        else:
            fp.write('0' + TAB + TAB)
	fp.write(omimdisease_urlprefix + CRT)

    reportlib.finish_nonps(fp)
    
    # re-create a symbolic link between the new file and the current file
    os.chdir(reportDir)
    os.symlink(reportWithDate + '.rpt', currentReport)

#
# main
#

#iphone_genes()
#iphone_mp()
iphone_omim()

