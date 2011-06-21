#!/usr/local/bin/python

'''
#
# Report:
#       Tab-delimited file of MGI coordinates for gene trap/gene unification pipeline
#	excludes DNA segments and QTLs
#	(TR 8120)
#
#   1. chromosome: "chr##"
#   2. source of feature: "MGI"
#   3. gene feature: marker type or "GeneModel"
#   4. start coordinate
#   5. end coordinate
#   6. empty: "."
#   7. strand
#   8. empty "."
#   9. ID=MGI ID;Name=marker symbol;Note=marker feature" or blank
#
# Usage:
#       MGI_GTGUP.py
#
# History:
#
# 07/28/2010 lec
#	- TR 6839/added marker features to field 9 as a "note"
#
# 01/23/2007 lec
#	- created
#
'''
 
import sys
import os
import string
import db
import reportlib
import mgi_utils

CRT = reportlib.CRT
TAB = reportlib.TAB

field1 = 'chr%s'
field2 = 'MGI'
field3 = 'GeneModel'
field6 = '.'
field8 = '.'
field9 = 'ID=%s;Name=%s'

idField = 'ID='

# Note=feature1,feature2,feature3,...
# concatenate notes by ','
noteField = ';Note='

# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], fileExt = '.gff', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# lookup of Marker Featers (MRK_MCV_Cache)
# direct features only
#

results = db.sql('''select _Marker_key, term 
		    from MRK_MCV_Cache where qualifier = 'D'
		 ''', 'auto')

mcvLookup = {}
for r in results:
    key = r['_Marker_key']
    value = r['term']

    if not mcvLookup.has_key(key):
	mcvLookup[key] = []
    mcvLookup[key].append(value)

# mapped features that have marker associations
# exclude DNA Segments and QTLs

results = db.sql('''select l.chromosome, l.strand, m.symbol, m._Marker_key,
	markerType = mt.name, a.accID, 
        startC = convert(int, l.startCoordinate), 
        endC = convert(int, l.endCoordinate) 
        from MRK_Location_Cache l, MRK_Marker m,
	     MRK_Types mt, ACC_Accession a, MRK_Chromosome c
        where m._Organism_key = 1 
	and l.startCoordinate is not null 
	and l._Marker_key = m._Marker_key 
	and m._Marker_Status_key in (1,3) 
	and l._Marker_Type_key not in (2,6) 
	and l._Marker_Type_key = mt._Marker_Type_key 
	and l._Marker_key = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 1 
	and a.prefixPart = "MGI:" 
	and a.preferred = 1 
	and m.chromosome = c.chromosome 
	and m._Organism_key = c._Organism_key 
	order by c.sequenceNum, l.startCoordinate, l.endCoordinate
	''', 'auto')

for r in results:

    key = r['_Marker_key']

    fp.write(field1 % (r['chromosome']) + TAB)
    fp.write(field2 + TAB)
    fp.write(r['markerType'] + TAB)
    fp.write(str(r['startC']) + TAB)
    fp.write(str(r['endC']) + TAB)
    fp.write(field6 + TAB)
    fp.write(str(r['strand']) + TAB)
    fp.write(field8 + TAB)
    fp.write(field9 % (r['accID'], r['symbol']))

    # add "note" if marker features exist
    if mcvLookup.has_key(key):
	fp.write(noteField + string.join(mcvLookup[key], ','))

    fp.write(CRT)

# mapped sequence features that don't have marker associations

results = db.sql('''select f.strand, c.chromosome, a.accID, 
        startC = convert(int, f.startCoordinate), 
        endC = convert(int, f.endCoordinate) 
        from MAP_Coordinate cc, MAP_Coord_Feature f, MRK_Chromosome c, ACC_Accession a 
        where f._MGIType_key = 19 
	and f._Object_key = a._Object_key 
	and a._MGIType_key = 19 
        and f._Map_key = cc._Map_key 
	and cc._MGIType_key = 27 
	and cc._Object_key = c._Chromosome_key 
	and not exists (select 1 from SEQ_Marker_Cache smc where f._Object_key = smc._Sequence_key) 
	order by c.sequenceNum, f.startCoordinate, f.endCoordinate
	''', 'auto')

for r in results:

    fp.write(field1 % (r['chromosome']) + TAB)
    fp.write(field2 + TAB)
    fp.write(field3 + TAB)
    fp.write(str(r['startC']) + TAB)
    fp.write(str(r['endC']) + TAB)
    fp.write(field6 + TAB)
    fp.write(str(r['strand']) + TAB)
    fp.write(field8 + TAB)
    fp.write(idField + r['accID'] + CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
