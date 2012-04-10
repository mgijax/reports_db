#!/usr/local/bin/python

'''
#
# TR 4970 (custom SQL TR 4338)
#
# Report:
#       Produce a tab-delimited report with the following output fields:
#
#       MGD ID
#       Gene Symbol
#       Gene Name
#       Chromosome
#       cM Position
#       Ensembl Gene ID
#	Ensembl Transcript IDs - space delimited
#	Ensembl Protein IDs - space delimited
# Usage:
#       MRK_ENSEMBL.py
#
# History:
#
# sc    03/12/2010
#       - TR9774 Add Ensembl and VEGA transcripts

# dbm	09/24/2008
#	- removed marker type restriction that was only looking for "genes"
#
# lec	08/26/2003
#	- added to reports_db product
#	- changed order of cM position and chromosome to be consistent
#	  with other Sequence reports.
#	- added outputdir, import os
#
# dbm   12/13/2002
#       - created
#
'''
 
import sys 
import os
import string
import db
import reportlib

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# Transcript ID lookup by Genomic ID
db.sql('select sa._Sequence_key_1 as transcriptKey, ' + \
        'sa._Sequence_key_2 as genomicKey ' + \
        'into #transGen ' + \
        'from SEQ_Sequence_Assoc sa ' + \
        'where sa._Qualifier_key = 5445464', None)  # transcribed from qualifier
db.sql('create index idx1 on ' + \
        '#transGen(transcriptKey)', None)
db.sql('create index idx2  on ' + \
        '#transGen(genomicKey)', None)
results = db.sql('select a1.accID as genomicID, a2.accID as transcriptID ' + \
        'from  #transGen t, ACC_Accession a1, ACC_Accession a2 ' + \
        'where t.genomicKey = a1._Object_key ' + \
        'and a1._MGIType_key = 19 ' + \
        'and a1.preferred = 1 ' + \
        'and t.transcriptKey = a2._Object_key ' + \
        'and a2._MGIType_key = 19 ' + \
        'and a2.preferred = 1 ' + \
        'order by a1.accID', 'auto')

genomicToTranscript = {}
for r in results:
    #print 'r:%s' % r
    key = r['genomicID']
    value = r['transcriptID']
    if not genomicToTranscript.has_key(key):
        genomicToTranscript[key] = []
    genomicToTranscript[key].append(value)

# Protein ID lookup by Genomic ID
db.sql('select tg.genomicKey, ' + \
    'sa._Sequence_key_1 as proteinKey ' + \
    'into #protGen ' + \
    'from #transGen tg, SEQ_Sequence_Assoc sa ' + \
    'where sa._Qualifier_key = 5445465 ' + \
    'and tg.transcriptKey =  sa._Sequence_key_2', None)

results = db.sql('select a1.accID as genomicID, a2.accID as proteinID ' + \
        'from  #protGen t, ACC_Accession a1, ACC_Accession a2 ' + \
        'where t.genomicKey = a1._Object_key ' + \
        'and a1._MGIType_key = 19 ' + \
        'and a1.preferred = 1 ' + \
        'and t.proteinKey = a2._Object_key ' + \
        'and a2._MGIType_key = 19 ' + \
        'and a2.preferred = 1 ' + \
        'order by a1.accID', 'auto')
genomicToProtein = {}
for r in results:
    key = r['genomicID']
    value = r['proteinID']
    if not genomicToProtein.has_key(key):
        genomicToProtein[key] = []
    genomicToProtein[key].append(value)

cmd = 'select a1.accID "MGI", m.symbol, m.name, m.chromosome, o.offset, a2.accID "Ensembl" ' + \
      'from ACC_Accession a1, ACC_Accession a2, MRK_Marker m, MRK_Offset o ' + \
      'where a1._Object_key = a2._Object_key and ' + \
            'a1._Object_key = m._Marker_key and ' + \
            'm._Marker_key = o._Marker_key and ' + \
            'a1._LogicalDB_key = 1 and ' + \
            'a1._MGIType_key = 2 and ' + \
            'a1.prefixPart = "MGI:" and ' + \
            'a1.preferred = 1 and ' + \
            'a2._LogicalDB_key = 60 and ' + \
            'a2._MGIType_key = 2 and ' + \
            'o.source = 0 ' + \
      'order by m.chromosome, m.symbol'

results = db.sql(cmd, 'auto')

for r in results:
    genomicID = r['Ensembl']
    fp.write(r['MGI'] + TAB + 
	    r['symbol'] + TAB + 
	    r['name'] + TAB +
	    str(r['offset']) + TAB +
            r['chromosome'] + TAB + 
            genomicID + TAB)
    if genomicToTranscript.has_key(genomicID):
        fp.write(string.join(genomicToTranscript[genomicID], ' '))
    fp.write(TAB)
    if genomicToProtein.has_key(genomicID):
        fp.write(string.join(genomicToProtein[genomicID], ' '))
    fp.write(CRT)
reportlib.finish_nonps(fp)
