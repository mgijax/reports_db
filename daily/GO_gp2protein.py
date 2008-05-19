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
#       field 2: Swiss-Prot sequence ('SWP:' + seqID) OR
#                TrEMBL sequence ('TR:' + seqID) OR
#                RefSeq sequence ('NCBI:' + seqID) OR
#
# Used by:
#
# History:
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
import db
import reportlib

#
# Main
#

TAB = reportlib.TAB
CRT = reportlib.CRT

fp = reportlib.init('gp2protein', fileExt = '.mgi', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

cmds = []

#
# all mouse genes with representative protein sequence ids
#
cmds.append('''
    select distinct mm._Marker_key, mm.mgiID, seqID=mc.accID, mc._LogicalDB_key
    into #results1
    from SEQ_Marker_Cache mc, MRK_Mouse_View mm
    where mc._Marker_key = mm._Marker_key
    and mm._Marker_Type_key = 1
    and mc._Qualifier_key = 615421
    ''')

cmds.append('''
    create index ix1 on #results1(_Marker_key)
    ''')

#
# all mouse genes not in the first group that have an Ensembl,
# NCBI, or VEGA gene model as the representative genomic sequence
#
cmds.append('''
    select distinct mm.mgiID
    into #results2
    from SEQ_Marker_Cache mc, MRK_Mouse_View mm
    where mc._Marker_key = mm._Marker_key
    and mm._Marker_Type_key = 1
    and mc._Qualifier_key = 615419
    and mc._LogicalDB_key in (59,60,85)
    and mm._Marker_key not in (
	select _Marker_key
	from #results1
	)
    ''')

cmds.append('''
    select * 
    from #results1
    ''')

cmds.append('''
    select * 
    from #results2
    ''')

#
results = db.sql(cmds, 'auto')

#
# Write a record to the report for each marker/sequence in the results set.
#
for r in results[3]:
    mgiID = "MGI:"+r['mgiID']
    logicalDB = r['_LogicalDB_key']

    #
    # Apply the proper prefix to the seq ID based on the logical DB.
    #
    if logicalDB in [13,41]:
        seqID = 'UniProtKB:' + r['seqID']
    else:
        seqID = 'NCBI:' + r['seqID']

    fp.write(mgiID + TAB + seqID + CRT)

#
#
#
for r in results[4]:
    mgiID = "MGI:"+r['mgiID']
    fp.write(mgiID + TAB + CRT)

reportlib.finish_nonps(fp)
