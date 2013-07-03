#!/usr/local/bin/python

'''
#
# MGI_InParanoid_FASTA.py
#
# Report:
#	TR 6500
#	Requested by InParanoid
#
#	Takes the FASTA files REFSEQFASTA and UNIPROTFASTA and generates
#	a new FASTA file that contains:
#
#	1. protein ids that are designated as
#           "representative polypeptide" for a given Marker
#	   in MGI.
#
#	2. a new FASTA header:
#
#	   MGI:#### source=MGI; version=MM/DD/YYYY; symbol=####; uniprot=####
#	   MGI:#### source=MGI; version=MM/DD/YYYY; symbol=####; refseq=####
#
# Usage:
#       MGI_InParanoid_FASTA.py
#
# History:
#
# 11/16/2007	dbm
#	- added new "aaseq.1090.fasta" file
#
# 12/12/2006	lec
#	- exclude markers annotated to "deleted" sequences
#	- don't exclude sequences annotated to > 1 marker
#
# 12/05/2006	lec
#	- added refseq proteins
#
# lec   01/25/2005
#       - created
#
'''

import sys
import os
import string
import mgi_utils
import reportlib
import pg_db
db = pg_db
db.setTrace()
db.setAutoTranslateBE()

uniprotHeader = '>%s source=MGI; version=%s; symbol=%s, uniprot=%s\n'
refseqHeader = '>%s source=MGI; version=%s; symbol=%s, refseq=%s\n'
fastaHeader = '>%s%s MGI:%s "%s"\n'
uniprotFileName = os.environ['UNIPROTFASTA']
refseqFileName = os.environ['REFSEQFASTA']
uniprotLabel = 'UniProt:'
refseqLabel = 'NCBI:'

reportNameA = 'Mus-musculus_MGI_' + mgi_utils.date('%m%d%Y') + '_protein-reps'
reportNameB = 'Mus-musculus_MGI_' + mgi_utils.date('%m%d%Y') + '_protein-all'
reportNameC = 'aaseq'
fileExtension = '.fa'
fpA = None
fpB = None
fpC = None

rep = {}
seqs = {}
markers = {}
mgiIDs = {}

def initialize():

    global fpA, fpB, fpC
    global rep, seqs, markers, mgiIDs

    fpA = reportlib.init(reportNameA, outputdir = os.environ['INPARANOIDDIR'], printHeading = None, fileExt = fileExtension)
    fpB = reportlib.init(reportNameB, outputdir = os.environ['INPARANOIDDIR'], printHeading = None, fileExt = fileExtension)
    fpC = reportlib.init(reportNameC, outputdir = os.environ['INPARANOIDDIR'], printHeading = None, fileExt = '.10090.fasta')

    # deleted sequences
    db.sql('''select s._Sequence_key into #deletedsequences 
	    from SEQ_Sequence s 
	    where s._SequenceStatus_key = 316343''', None)
    db.sql('create index deletedsequences_idx1 on #deletedsequences(_Sequence_key)', None)

    # markers with deleted sequences
    db.sql('''select distinct m._Marker_key 
	    into #excludemarkers 
	    from #deletedsequences d, SEQ_Marker_Cache m 
	    where d._Sequence_key = m._Sequence_key
	    and m._Organism_key = 1''', None)
    db.sql('create index excludedmarkers_idx1 on #excludemarkers(_Marker_key)', None)

    #
    # marker/polypeptide sequences that are not excluded
    # markers of type gene only
    #
    db.sql('''select s.accID, s._Marker_key, s._Qualifier_key into #sequences 
	    from SEQ_Marker_Cache s, MRK_Marker m 
            where s._SequenceType_key = 316348 
	    and s._Organism_key = 1 
	    and s._Marker_key = m._Marker_key
	    and m._Marker_Type_key = 1
	    and m._Marker_Status_key in (1,3) 
	    and not exists (select 1 from #excludemarkers x where 
		s._Marker_key = x._Marker_key)''', None)
    db.sql('create index sequences_idx1 on #sequences(accID)', None)
    db.sql('create index sequences_idx2 on #sequences(_Marker_key)', None)

    #
    # cache the representative polypeptide
    #

    results = db.sql('''select accID, _Marker_key 
	from #sequences 
	where _Qualifier_key = 615421''', 'auto')
    for r in results:
        key = r['accID']
        value = r['_Marker_key']
        rep[key] = value

    #
    # cache the seq id:marker key
    #

    results = db.sql('select accID, _Marker_key from #sequences', 'auto')
    for r in results:
        key = r['accID']
        value = r['_Marker_key']
	if not seqs.has_key(key):
	    seqs[key] = []
        seqs[key].append(value)

    #
    # cache the marker key:symbol
    #

    results = db.sql('''select s._Marker_key, m.symbol 
	    from #sequences s, MRK_Marker m 
	    where s._Marker_key = m._Marker_key''', 'auto')
    for r in results:
        key = r['_Marker_key']
        value = r['symbol']
        markers[key] = value

    #
    # cache the marker key:mgi id
    #

    results = db.sql('''select s._Marker_key, a.accID 
	    from #sequences s, ACC_Accession a 
	    where s._Marker_key = a._Object_key
	    and a._MGIType_key = 2 
	    and a._LogicalDB_key = 1 
	    and a.prefixPart = 'MGI:' 
	    and a.preferred = 1''', 'auto')
    for r in results:
        key = r['_Marker_key']
        value = r['accID']
        mgiIDs[key] = value

def process(inFileName, column, header, label):

    inFile = open(inFileName, 'r')
    skipRecord = 0

    for line in inFile.readlines():

        # a header line is one the starts with ">"

        if line[0] == '>':

            skipRecord = 0
	    tokens1 = string.split(line[:-1], '|')
	    tokens2 = string.split(tokens1[column], '.')
	    seqID = tokens2[0]

	    # if the seq ID in the FASTA file is a representative polypeptide, then process id

	    if seqs.has_key(seqID):
		for markerKey in seqs[seqID]:
	            symbol = markers[markerKey]
	            mgiID = mgiIDs[markerKey]
                    newLine = header % (mgiID, mgi_utils.date('%m/%d/%Y'), symbol, seqID)
	            if rep.has_key(seqID):
		        fpA.write(newLine)
		        fpC.write(fastaHeader % (label, seqID, mgiID, symbol))
	            fpB.write(newLine)

	    # else, skip until we find the next header record

            else:
	        skipRecord = 1

        # if the record ain't skipped, print out the sequence

        else:
	    if not skipRecord:
	        if rep.has_key(seqID):
                    fpA.write(line)
                    fpC.write(line)
                fpB.write(line)

    inFile.close()

#
# main
#

initialize()
process(uniprotFileName, 1, uniprotHeader, uniprotLabel)
process(refseqFileName, 3, refseqHeader, refseqLabel)
reportlib.finish_nonps(fpA)
reportlib.finish_nonps(fpB)
reportlib.finish_nonps(fpC)
