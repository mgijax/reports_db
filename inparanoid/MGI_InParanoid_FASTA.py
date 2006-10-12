#!/usr/local/bin/python

'''
#
# MGI_InParanoid_FASTA.py
#
# Report:
#	TR 6500
#	Requested by InParanoid
#
#	Takes the FASTA file ${INPARANOIDFASTA} and generates
#	a new FASTA file that contains:
#
#		1. uniprot ids that are designated as
#	           "representative polypeptide" for a given Marker
#		   in MGI.
#	
#		2. a new FASTA header:
#
#		   MGI:#### source=MGI; version=MM/DD/YYYY; symbol=####; uniprot=####
#
# Usage:
#       MGI_InParanoid_FASTA.py
#
# History:
#
# lec   01/25/2005
#       - created
#
'''

import sys
import os
import string
import db
import mgi_utils
import reportlib

newHeader = '>%s source=MGI; version=%s; symbol=%s, uniprot=%s\n'
inFileName = os.environ['INPARANOIDFASTA']

reportNameA = 'Mus-musculus_MGI_' + mgi_utils.date('%m%d%Y') + '_protein-reps'
reportNameB = 'Mus-musculus_MGI_' + mgi_utils.date('%m%d%Y') + '_protein-all'
fileExtension = '.fa'
fpA = reportlib.init(reportNameA, outputdir = os.environ['INPARANOIDDIR'], printHeading = None, fileExt = fileExtension)
fpB = reportlib.init(reportNameB, outputdir = os.environ['INPARANOIDDIR'], printHeading = None, fileExt = fileExtension)

#
# select all representative polypeptide sequences
#

db.sql('select a.accID, s._Marker_key, s._Qualifier_key into #allsequences ' + \
	'from SEQ_Marker_Cache s, ACC_Accession a ' + \
        'where s._Sequence_key = a._Object_key ' + \
        'and a._MGIType_key = 19 ' + \
        'and a.preferred = 1 ' + \
        'and a._LogicalDB_key in (13, 41)', None)
db.sql('create index idx1 on #allsequences(accID)', None)

#
# eliminate those sequences annotated to > 1 marker
#

db.sql('select accID, _Marker_key, _Qualifier_key into #sequences from #allsequences group by accID having count(*) = 1', None)
db.sql('create index idx1 on #sequences(_Marker_key)', None)

#
# cache the representative polypeptide
#

rep = {}
results = db.sql('select accID, _Marker_key from #sequences where _Qualifier_key = 615421', 'auto')
for r in results:
    key = r['accID']
    value = r['_Marker_key']
    rep[key] = value

#
# cache the seq id:marker key
#

seqs = {}
results = db.sql('select accID, _Marker_key from #sequences', 'auto')
for r in results:
    key = r['accID']
    value = r['_Marker_key']
    seqs[key] = value

#
# cache the marker key:symbol
#

markers = {}
results = db.sql('select s._Marker_key, m.symbol ' + \
	'from #sequences s, MRK_Marker m ' + \
	'where s._Marker_key = m._Marker_key', 'auto')
for r in results:
    key = r['_Marker_key']
    value = r['symbol']
    markers[key] = value

#
# cache the marker key:mgi id
#

mgiIDs = {}
results = db.sql('select s._Marker_key, a.accID ' + \
	'from #sequences s, ACC_Accession a ' + \
	'where s._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1', 'auto')
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    mgiIDs[key] = value

inFile = open(inFileName, 'r')
skipRecord = 0

for line in inFile.readlines():

    # a header line is one the starts with ">"

    if line[0] == '>':
	tokens = string.split(line[:-1], '|')
	seqID = tokens[1]

	# if the seq ID in the FASTA file is a representative polypeptide, then process id

	if seqs.has_key(seqID):
	    markerKey = seqs[seqID]
	    symbol = markers[markerKey]
	    mgiID = mgiIDs[markerKey]
            newLine = newHeader % (mgiID, mgi_utils.date('%m/%d/%Y'), symbol, seqID)
	    if rep.has_key(seqID):
		fpA.write(newLine)
	    fpB.write(newLine)
	    skipRecord = 0	# re-set skip with every new header

	# else, skip until we find the next header record

        else:
	    skipRecord = 1

    # if the record ain't skipped, print out the sequence

    else:
	if not skipRecord:
	    if rep.has_key(seqID):
                fpA.write(line)
            fpB.write(line)

inFile.close()

reportlib.finish_nonps(fpA)
reportlib.finish_nonps(fpB)

