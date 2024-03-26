
'''
#
# GO_gp2protein.py (TR 4877, originally TR 3659)
#
# Report:
#       Tab-delimited file
#
# Usage:
#       GO_gp2protein.py
#
# Output
#
#	A tab-delimited file in this format:
#	field 1: MGI Marker ID
#       field 2: protein or transcript or neither (blank)
#
#	gp2protein.mgi : genes witih protein sequences (fp1)
#	gp2rna.mgi : any genes with/without a transcript (fp2)
#	gp_unlocalized.mgi : all genes not in gp2protein or gp2rna (fp3)
#
# History:
#
# 01/30/2015	sc
#	- TR11879  Fix: non-protein coding genes in gp2protein file
#	- exclude RNA genes associated with uniprot in the gp2protein file
#
# 10/31/2013	lec
#	- TR11510/re-organize reports to include 'feature type' rules
#
# 02/13/2013	lec
#	- TR11272/split protein coding genes into protein coding/rna/neither reports
#
# 12/14/2012	lec
#	- TR11242/inclulde _Marker_Status_key = 1
#
# 05/22/2012	lec
#	- TR11034/UniProt trumps protein; add uniprotkbDict
#
# 05/01/2012	lec
#	- make 'idx1' names unique
#
# 05/01/2012	sc
#	- TR11057
#	- miRNA Transcript IDs were not being written to file as their ldb
#	  was not being considered
#	- added 'else: continue' in the block of code which assigns prefixes
#	- added transcript ldbs to same block of code
#	- Harold wants all Genes with rep genomic sequence to appear in this rpt
#	  if protein, use that
#	  else if transcript,  use that
#	  else just print out marker accID
#
# lec	03/30/2011
#   - TR10652/change 'NCBI:' to 'RefSeq:'
#
# 08/18/2010	lec
#       - TR6839/marker types
#       - marker type 11 (microRNA) moved to marker type 1 (gene)
#
# 08/17/2010	lec
#	- TR 10321/fix prefixes for Ensembl and Vega
#
# 10/21/2008	lec
#	- TR 9325; add microRNA marker types with transcript
#
# 05/02/2008	jer
#	- TR 8994; yet another rewrite; changes to selection logic
#	as well as to output formatting. Want all coding genes, whether
#	or not they have representative protein sequence (leave col 2
#	blank if not). To do this: first get everything that has a
#	representative protein seq; then add in everything else that
#	has a Ensembl, or NCBI gene model. Output formatting:
#	change "MGI:12345" to "MGI:MGI:12345"; and change "SP:" and
#	"TR:" to "UniProtKB:".
#
# 12/24/2007    dbm
#       - TR 8697; complete re-write; Show all protein coding genes and
#         their representative protein sequences, regardless of annotations.
#
# 6/6/2003    lec
#	- Created for TR3659.  
#	Report the MGI ID, SwissProt ID(s), and GO ID(s) for all mouse
#	markers that have associated SwissProt sequences.

'''
 
import sys
import os
import reportlib
import db

db.setTrace()

#
# Main
#

TAB = reportlib.TAB
CRT = reportlib.CRT

fp1 = reportlib.init('gp2protein', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp2 = reportlib.init('gp2rna', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp3 = reportlib.init('gp_unlocalized', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# all markers of type 'gene'
# with feature type (see below)
#
db.sql('''
    select distinct mm._Marker_key, a.accID, tdc.term
    into temporary table allgenes
    from MRK_Marker mm, ACC_Accession a, VOC_Annot_View tdc
    where mm._Marker_Status_key = 1
    and mm._Marker_Type_key = 1
    and mm._Organism_key = 1
    and mm._Marker_key = a._Object_key
    and a._MGIType_key = 2
    and a._LogicalDB_key = 1
    and a.preferred = 1
    and mm._Marker_key = tdc._Object_key
    and tdc._AnnotType_key = 1011    
    and tdc._LogicalDB_key = 146    
    and (tdc.term in ('protein coding gene') or tdc.term like '%RNA gene')
    ''', None)
db.sql('create index allgenes_idx1 on allgenes(_Marker_key)', None)

#
# representative transcripts by marker
# _Qualifier_key = 615420 = 'transcript'
#
db.sql('''
    select distinct mm._Marker_key, mm.accID, mc.accID as seqID, mc._LogicalDB_key
    into temporary table transcripts
    from SEQ_Marker_Cache mc, allgenes mm
    where mc._Marker_key = mm._Marker_key
    and mc._Qualifier_key = 615420
    ''', None)
db.sql('create index transcripts_idx1 on transcripts(_Marker_key)', None)

results = db.sql('''select * from transcripts''', 'auto')
transcriptDict = {}
for r in results:
    transcriptDict [r['_Marker_key']] = [ r['seqID'],r['_LogicalDB_key'] ]

#
# representative proteins by marker
# _Qualifier_key = 615421 = 'polypeptide'
#
db.sql('''
    select distinct mm._Marker_key, mm.accID, mc.accID as seqID, mc._LogicalDB_key
    into temporary table proteins
    from SEQ_Marker_Cache mc, allgenes mm
    where mc._Marker_key = mm._Marker_key
    and mc._Qualifier_key = 615421
    ''', None)
db.sql('create index proteins_idx1 on proteins(_Marker_key)', None)

results = db.sql('''select * from proteins''', 'auto')
proteinDict = {}
for r in results:
    proteinDict[r['_Marker_key']] = [ r['seqID'],r['_LogicalDB_key'] ]

#
# uniprotKB by marker (ldb = 13 SwissProt only)
#
db.sql('''
    select distinct mm._Marker_key, mm.accID, mc.accID as seqID, mc._LogicalDB_key
    into temporary table uniprotkb
    from SEQ_Marker_Cache mc, allgenes mm
    where mc._Marker_key = mm._Marker_key
    and mc._LogicalDB_key in (13)
    ''', None)
db.sql('create index uniprotkb_idx1 on proteins(_Marker_key)', None)

results = db.sql('''select * from uniprotkb''', 'auto')
uniprotkbDict = {}
for r in results:
    uniprotkbDict[r['_Marker_key']] = [ r['seqID'],r['_LogicalDB_key'] ]

#
# all mouse genes
# include: markers of type 'gene'
# include: marker with 'genomic' representative sequence from Ensembl, NCBI gene model
# excludd: feature type in 'protein coding gene'
# excludd: feature type in '%RNA gene'
# exclude: markers in transcriptDict or proteinDict list
#
# _Qualifier_key = 615419 = 'genomic'
#
db.sql('''
    select distinct a.accID
    into temporary table noTransProt
    from SEQ_Marker_Cache mc, MRK_Marker mm, ACC_Accession a, VOC_Annot_View tdc
    where mm._Marker_Status_key = 1
    and mm._Marker_Type_key = 1
    and mm._Organism_key = 1
    and mm._Marker_key = a._Object_key
    and a._MGIType_key = 2
    and a._LogicalDB_key = 1
    and a.preferred = 1
    and mm._Marker_key = mc._Marker_key
    and mc._Qualifier_key = 615419
    and mc._LogicalDB_key in (59,60)
    and mm._Marker_key not in (select _Marker_key from transcripts)
    and mm._Marker_key not in (select _Marker_key from proteins)
    and mm._Marker_key = tdc._Object_key
    and tdc._AnnotType_key = 1011    
    and tdc._LogicalDB_key = 146    
    and tdc.term not in ('protein coding gene') 
    and tdc.term not like '%RNA gene'
    ''', None)

#
# Write a record to the report for each marker/sequence in the results set.
#
# if uniprotKB exists, use it
# else if protein exists, use it
# else use transcript
#
# if the sequence accession id for the given marker is not of interest, then skip
# of interest list:  13,41,9,27,134
#

results = db.sql('select * from allgenes', 'auto')
for r in results:

    accID = "MGI:" + r['accID']
    tdcTerm = r['term']
    markerKey = r['_Marker_key']
    featureType = r['term']
    isfp1 = 0

    # fp1 (protein)
    if markerKey in uniprotkbDict: 
        if tdcTerm == 'protein coding gene':
            l = uniprotkbDict[markerKey]
            seqID = l[0]
            logicalDB = l[1]
            isfp1 = 1
        #else:
        #    print(('RNA w/uniprot seq assoc: %s' % accID))
    elif markerKey in proteinDict:
        if tdcTerm == 'protein coding gene':
            l = proteinDict[markerKey]
            seqID = l[0]
            logicalDB = l[1]
            isfp1 = 1
        #else:
        #    print(('RNA w/rep protein seq assoc: %s' % accID))

    # fp2 (rna)
    elif markerKey in transcriptDict:
        l = transcriptDict[markerKey]
        seqID = l[0]
        logicalDB = l[1]

    # if there is no protein or transcript 
    # 	AND the gene is *not* featureType = RNA
    # then skip it
    elif featureType == 'protein coding gene':
        continue

    # else...there is no protein or transcript AND the gene *is* featureType = RNA
    # so we keep it
    else:
        logicalDB = ''

    #
    # Apply the proper prefix to the seq ID based on the logical DB.
    #

    if logicalDB in [13,41]:
        seqID = 'UniProtKB:' + seqID
    elif logicalDB in [9]:
        seqID = 'EMBL:' + seqID
    elif logicalDB in [27]:
        seqID = 'RefSeq:' + seqID
    elif logicalDB in [133, 134]:
        seqID = 'ENSEMBL:' + seqID
    else:
        seqID = ''

    if isfp1:
        fp1.write(accID + TAB + seqID + CRT)
    else:
        fp2.write(accID + TAB + seqID + CRT)

#
# markers that have neither transcript nor protein
#
results = db.sql('select * from noTransProt', 'auto')
for r in results:
    accID = "MGI:" + r['accID']
    fp3.write(accID + TAB + CRT)

reportlib.finish_nonps(fp1)
reportlib.finish_nonps(fp2)
reportlib.finish_nonps(fp3)

