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
    select distinct mc._Marker_key, mc._Marker_Type_key, mm.mgiID, mc.accID as seqID, mc._LogicalDB_key
    into #results1
    from SEQ_Marker_Cache mc, MRK_Mouse_View mm
    where mc._Marker_key = mm._Marker_key
    and mc._Marker_Type_key = 1
    and mc._Qualifier_key = 615421
    union
    select distinct mc._Marker_key, mc._Marker_Type_key, mm.mgiID, mc.accID as seqID, mc._LogicalDB_key
    from SEQ_Marker_Cache mc, MRK_Mouse_View mm, MRK_MCV_Cache mcv
    where mc._Marker_key = mm._Marker_key
    and mc._Marker_Type_key = 1
    and mc._Qualifier_key = 615420
    and mc._Marker_key = mcv._Marker_key
    and mcv.term = "miRNA Gene"
    ''', None)

db.sql('create index ix1 on #results1(_Marker_key)', None)

#
# all mouse genes not in the first group that have an Ensembl,
# NCBI, or VEGA gene model as the representative genomic sequence
#
db.sql('''
    select distinct mm.mgiID
    into #results2
    from SEQ_Marker_Cache mc, MRK_Mouse_View mm
    where mc._Marker_key = mm._Marker_key
    and mc._Marker_Type_key = 1
    and mc._Qualifier_key = 615419
    and mc._LogicalDB_key in (59,60,85)
    and mm._Marker_key not in (select _Marker_key from #results1)
    ''', None)

#
# Write a record to the report for each marker/sequence in the results set.
#
results = db.sql('select * from #results1', 'auto')
for r in results:
    mgiID = "MGI:"+r['mgiID']
    logicalDB = r['_LogicalDB_key']

    #
    # Apply the proper prefix to the seq ID based on the logical DB.
    #
    if logicalDB in [13,41]:
        seqID = 'UniProtKB:' + r['seqID']
    elif logicalDB in [9]:
        seqID = 'EMBL:' + r['seqID']
    elif logicalDB in [27]:
        seqID = 'RefSeq:' + r['seqID']
    elif logicalDB in [132]:
        seqID = 'VEGA:' + r['seqID']
    elif logicalDB in [134]:
        seqID = 'ENSEMBL:' + r['seqID']

    fp.write(mgiID + TAB + seqID + CRT)

#
#
#
results = db.sql('select * from #results2', 'auto')
for r in results:
    mgiID = "MGI:"+r['mgiID']
    fp.write(mgiID + TAB + CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(1)

