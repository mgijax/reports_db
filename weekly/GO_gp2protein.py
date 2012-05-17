#!/usr/local/bin/python

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
#       field 2: SwissProt or TrEMBL sequences ('UniPrpt:' + seqID)
#                RefSeq sequence ('RefSeq:' + seqID)
#                GenBank sequence ('EMBL:' + seqID)
#
# Used by:
#
# History:
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
#	  else just print out marker mgiID
#
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
#	has a VEGA, Ensembl, or NCBI gene model. Output formatting:
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
import string
import os
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


#
# Main
#

TAB = reportlib.TAB
CRT = reportlib.CRT

fp = reportlib.init('gp2protein', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

db.useOneConnection(1)

#
# all mouse genes with representative protein sequence ids
# all mouse genes that are microRNA's
#
db.sql('''
    select distinct mc._Marker_key, mc._Marker_Type_key, mm.mgiID, mc.accID as seqID, mc._LogicalDB_key, mc._Qualifier_key
    into #transcripts
    from SEQ_Marker_Cache mc, MRK_Mouse_View mm
    where mc._Marker_key = mm._Marker_key
    and mc._Marker_Type_key = 1
    and mc._Qualifier_key = 615420''', None)

results = db.sql('''select * from #transcripts''', 'auto')
# rep transcript by markerKey
transcriptDict = {}
for r in results:
    transcriptDict [r['_Marker_key']] = [ r['seqID'],r['_LogicalDB_key'] ]

db.sql('create index transcripts_idx1 on #transcripts(_Marker_key)', None)

db.sql('''
    select distinct mc._Marker_key, mc._Marker_Type_key, mm.mgiID, mc.accID as seqID, mc._LogicalDB_key, mc._Qualifier_key
    into #proteins
    from SEQ_Marker_Cache mc, MRK_Mouse_View mm
    where mc._Marker_key = mm._Marker_key
    and mc._Marker_Type_key = 1
    and mc._Qualifier_key = 615421''', None)

results = db.sql('''select * from #proteins''', 'auto')

# rep protein by markerKey (not all markers are protein coding)
proteinDict = {}
for r in results:
    proteinDict[r['_Marker_key']] = [ r['seqID'],r['_LogicalDB_key'] ]

db.sql('create index proteins_idx1 on #proteins(_Marker_key)', None)

#
# all mouse genes not in the first two group that have an Ensembl,
# NCBI, or VEGA gene model as the representative genomic sequence
#
db.sql('''
    select distinct mm.mgiID
    into #noTransProt
    from SEQ_Marker_Cache mc, MRK_Mouse_View mm
    where mc._Marker_key = mm._Marker_key
    and mc._Marker_Type_key = 1
    and mc._Qualifier_key = 615419
    and mc._LogicalDB_key in (59,60,85)
    and mm._Marker_key not in (select _Marker_key from #transcripts)
    and mm._Marker_key not in (select _Marker_key from #proteins)
    ''', None)

#
# Write a record to the report for each marker/sequence in the results set.
# If there is a protein,  use it, if not use the transcript
results = db.sql('select * from #transcripts', 'auto')
for r in results:
    mgiID = "MGI:"+r['mgiID']
    markerKey = r['_Marker_key']
    if proteinDict.has_key(markerKey):
	l = proteinDict[markerKey]
	seqID = l[0]
	logicalDB = l[1]
    else:
	l = transcriptDict[markerKey]
	seqID = l[0]
        logicalDB = l[1]
    #
    # Apply the proper prefix to the seq ID based on the logical DB.
    #
    if logicalDB in [13,41]:
        seqID = 'UniProtKB:' + seqID
    elif logicalDB in [9]:
        seqID = 'EMBL:' + seqID
    elif logicalDB in [27]:
        seqID = 'RefSeq:' + seqID
    elif logicalDB in [131,132]:
        seqID = 'VEGA:' + seqID
    elif logicalDB in [133, 134]:
        seqID = 'ENSEMBL:' + seqID
    else:
	continue
    fp.write(mgiID + TAB + seqID + CRT)

#
#
#
results = db.sql('select * from #noTransProt', 'auto')
for r in results:
    mgiID = "MGI:"+r['mgiID']
    fp.write(mgiID + TAB + CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(1)

