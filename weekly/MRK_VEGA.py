#!/usr/local/bin/python

'''
#
# TR 8055
#
# Report:
#       Produce a tab-delimited report with the following output fields:
#
#       MGD ID
#       Gene Symbol
#       Gene Name
#       Chromosome
#       cM Position
#       VEGA Gene ID
#	VEGA Transcript IDs - space delimited
#	VEGA Protein IDs - space delimited
#	Feature Type(s)
#	Coordinate Start
#	Coordinate End
#	Stand
#	BioType(s)
#
# Usage:
#       MRK_VEGA.py
#
# History:
#
# lec	04/23/2012
#	- TR11035/postgres cleanup; added feature type, coordinates & biotype
#
# sc    03/12/2010
#       - TR9774 Add Ensembl and VEGA transcripts
#
# lec   12/13/2006
#       - created
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
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db


CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# Transcript ID lookup by Genomic ID
#
db.sql('''
        select sa._Sequence_key_1 as transcriptKey, 
        sa._Sequence_key_2 as genomicKey 
        into #transGen 
        from SEQ_Sequence_Assoc sa 
        where sa._Qualifier_key = 5445464
	''', None)  # transcribed from qualifier
db.sql('create index idx1 on #transGen(transcriptKey)', None)
db.sql('create index idx2  on #transGen(genomicKey)', None)

results = db.sql('''
	select a1.accID as genomicID, a2.accID as transcriptID 
	from #transGen t, ACC_Accession a1, ACC_Accession a2 
	where t.genomicKey = a1._Object_key 
	and a1._MGIType_key = 19 
	and a1.preferred = 1 
	and t.transcriptKey = a2._Object_key 
        and a2._MGIType_key = 19 
        and a2.preferred = 1 
	order by a1.accID
	''', 'auto')

genomicToTranscript = {}
for r in results:
    #print 'r:%s' % r
    key = r['genomicID']
    value = r['transcriptID']
    if not genomicToTranscript.has_key(key):
        genomicToTranscript[key] = []
    genomicToTranscript[key].append(value)

#
# Protein ID lookup by Genomic ID
#
db.sql('''
    select tg.genomicKey, 
    sa._Sequence_key_1 as proteinKey 
    into #protGen 
    from #transGen tg, SEQ_Sequence_Assoc sa 
    where sa._Qualifier_key = 5445465 
    and tg.transcriptKey =  sa._Sequence_key_2
    ''', None)

results = db.sql('''
	select a1.accID as genomicID, a2.accID as proteinID 
        from #protGen t, ACC_Accession a1, ACC_Accession a2 
        where t.genomicKey = a1._Object_key 
        and a1._MGIType_key = 19 
        and a1.preferred = 1 
        and t.proteinKey = a2._Object_key 
        and a2._MGIType_key = 19 
        and a2.preferred = 1 
        order by a1.accID
	''', 'auto')
genomicToProtein = {}
for r in results:
    key = r['genomicID']
    value = r['proteinID']
    if not genomicToProtein.has_key(key):
        genomicToProtein[key] = []
    genomicToProtein[key].append(value)

#
# markers
#
db.sql('''
      select a1.accID as mgiID, 
	     m._Marker_key, 
	     m.symbol, 
	     m.name, 
	     m.chromosome, 
	     o.offset, 
	     a2.accID as vegaID
      into #markers
      from ACC_Accession a1, ACC_Accession a2, MRK_Marker m, MRK_Offset o 
      where a1._Object_key = a2._Object_key 
            and a1._Object_key = m._Marker_key 
            and m._Marker_key = o._Marker_key 
            and a1._LogicalDB_key = 1 
            and a1._MGIType_key = 2 
            and a1.prefixPart = 'MGI:' 
            and a1.preferred = 1 
            and a2._LogicalDB_key = 85 
            and a2._MGIType_key = 2 
            and o.source = 0 
      ''', None)
db.sql('create index marker_idx1 on #markers(_Marker_key)', None)

#
# feature types
#
results = db.sql('''
	select m._marker_key, s.term 
        from #markers m, MRK_MCV_Cache s 
        where m._marker_key = s._marker_key 
        and s.qualifier = 'D'
	''', 'auto')
featureTypes = {}
for r in results:
        key = r['_marker_key']
        value = r['term']
        if not featureTypes.has_key(key):
                featureTypes[key] = []
        featureTypes[key].append(value)

#
# coordinates
#
results = db.sql('''	
    select m._marker_key,
           c.strand, 
	   convert(int, c.startCoordinate) as startC,
	   convert(int, c.endCoordinate) as endC
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
# biotype
#
results = db.sql('''
	select distinct m._Marker_key, s.rawbiotype 
	from #markers m, SEQ_Marker_Cache s 
	where m._Marker_key = s._Marker_key 
	and s.rawbiotype is not null
	''', 'auto')
bioTypes = {}
for r in results:
	key = r['_Marker_key']
	value = r['rawbiotype']
	if not bioTypes.has_key(key):
		bioTypes[key] = []
	bioTypes[key].append(value)

#
# final select
#
results = db.sql('select * from #markers order by chromosome, symbol', 'auto')

for r in results:

    genomicID = r['vegaID']

    # column 1: mgi id
    # column 2: symbol
    # column 3: name
    # column 4: cM position
    # column 5: chromosome
    # column 6: genomic id

    fp.write(r['mgiID'] + TAB + 
	     r['symbol'] + TAB + 
	     r['name'] + TAB +
	     str(r['offset']) + TAB +
             r['chromosome'] + TAB + 
             genomicID + TAB)

    # column 7
    if genomicToTranscript.has_key(genomicID):
	fp.write(string.join(genomicToTranscript[genomicID], ' '))
    fp.write(TAB)

    # column 8
    if genomicToProtein.has_key(genomicID):
	fp.write(string.join(genomicToProtein[genomicID], ' '))
    fp.write(TAB)

    # column 9: feature types
    if featureTypes.has_key(r['_Marker_key']):	
	fp.write(string.join(featureTypes[r['_Marker_key']], '|'))
    fp.write(TAB)

    # column 10-11-12
    if coords.has_key(key):
        fp.write(mgi_utils.prvalue(coords[r['_Marker_key']][0]['startC']) + TAB)
        fp.write(mgi_utils.prvalue(coords[r['_Marker_key']][0]['endC']) + TAB)
        fp.write(mgi_utils.prvalue(coords[r['_Marker_key']][0]['strand']) + TAB)
    else:
        fp.write(TAB + TAB + TAB)

    # column 13: biotypes
    if bioTypes.has_key(r['_Marker_key']):	
	fp.write(string.join(bioTypes[r['_Marker_key']], '|'))
    fp.write(CRT)

reportlib.finish_nonps(fp)
