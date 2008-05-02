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
#	- TR 8994; several small changes; (1) add MGI: prefix to MGI ids
#	(e.g. MGI:MGI:12345); include all genes with a DNA/RNA 
#	sequence, even if no protein seqs are available
#	(leave col 2 blank in those cases); change "SWP:" and "TR:" prefixes
#	to "UniProtKB:".
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
# Get the representative sequence (_Qualifier_key = 615421) for each marker
# if the sequence is:
# Swiss-Prot (logical DB = 13)
# TrEMBL (logical DB = 41
# RefSeq (logical DB = 27 and accID starts with 'NP_' or 'XP_'
#
cmds.append('select mc._Marker_key, ' + \
                   'mc.accID, ' + \
                   'mc._LogicalDB_key ' + \
            'into #markerseq ' + \
            'from ACC_Accession a, SEQ_Marker_Cache mc ' + \
            'where a._MGIType_key = 2 and ' + \
                  '(a._LogicalDB_key in (13,41) or ' + \
                  '(a._LogicalDB_key = 27 and ' + \
                   'a.prefixPart in ("NP_","XP_"))) and ' + \
                  'a._Object_key = mc._Marker_key and ' + \
                  'a.accID = mc.accID and ' + \
                  'mc._Qualifier_key = 615421 and ' + \
                  'mc._Organism_key = 1')

#
# Get the MGI ID for each of the markers.
#
cmds.append('select a.accID "mgiID", ' + \
                   'ms.accID "seqID", ' + \
                   'ms._LogicalDB_key ' + \
            'from ACC_Accession a, #markerseq ms ' + \
            'where a._Object_key = ms._Marker_key and ' + \
                  'a._MGIType_key = 2 and ' + \
                  'a._LogicalDB_key = 1 and ' + \
                  'a.prefixPart = "MGI:" and ' + \
                  'a.preferred = 1 ' + \
            'order by a.accID')

#
# Get the genes that have at least one nucleic acid sequence, but no
# protein sequence. 
#
cmds.append('''
	select mv.mgiID
	from MRK_Mouse_View mv
	where mv._Marker_key not in (
		select ms._Marker_key
		from #markerseq ms
		)
	and mv._Marker_key in (
		select distinct mc._Marker_key
		from SEQ_Marker_Cache mc
		where mc._Organism_key = 1
		and mc._SequenceType_key in (316373,316346)
		)
	and mv.markerType = 'Gene'
	and mv.status = 'official'
	''')

#
# Get the results set.
#
results = db.sql(cmds, 'auto')

#
# Write a record to the report for each marker/sequence in the results set.
#
for r in results[1]:
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
for r in results[2]:
    mgiID = "MGI:"+r['mgiID']
    fp.write(mgiID + TAB + CRT)

reportlib.finish_nonps(fp)
