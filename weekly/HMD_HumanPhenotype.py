#!/usr/local/bin/python

'''
#
# HMD_HumanPhenotype.py 
#
# Report:
#     A tab-delimited report of all mouse-human orthologous gene pairs, 
#     with all MP IDs associated with any allele of the mouse gene.
#
#     Output fields: 
#     
#     Human Gene Symbol 
#     Human EntrezGene ID 
#     Mouse Gene Symbol 
#     MGI Accession ID 
#     MP IDs associated with any allele of the mouse gene (space- or comma-separated) 
#
#     Sorted by:
#     Human Gene Symbol. 
#    
# Usage:
#       HMD_HumanPhenotype.py
#
# History:
#
# lec	04/23/2008
#	- TR 8963; only select primary MP Header Ids
#
# lec	01/04/2004
#	- TR 5939; EntrezGene->EntrezGene
#
# lec	08/27/2003
#	- TR 5082; make a public report
#
# dm    07/14/03
#       -Inital Creation
# 
'''
 
import sys
import os
import db
import reportlib
import string

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

def createDict(results, keyField, valueField):

        d = {}
        for r in results:
                key = r[keyField]
                value = r[valueField]
                if not d.has_key(key):
                        d[key] = ' '
                d[key] = string.join([d[key],value,])
        return d


#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# Get mouse to human orthologous marker pair's symbols and keys
#

db.sql('select distinct mouseKey = h1._Marker_key, mouseSym = m1.symbol, ' + \
		'humanKey = h2._Marker_key, humanSym = m2.symbol ' + \
		'into #homology ' + \
		'from MRK_Homology_Cache h1, MRK_Homology_Cache h2, ' + \
		'MRK_Marker m1, MRK_Marker m2 ' + \
		'where h1._Organism_key = 1 ' + \
		'and h1._Class_key = h2._Class_key ' + \
		'and h2._Organism_key = 2 ' + \
		'and h1._Marker_key = m1._Marker_key ' + \
		'and h2._Marker_key = m2._Marker_key ', None)

db.sql('create nonclustered index index_mouseKey on #homology(mouseKey)', None)
db.sql('create nonclustered index index_humanKey on #homology(humanKey)', None)
db.sql('create nonclustered index index_humanSym on #homology(humanSym)', None)

#
# Get the MGI IDs for the mouse markers
#

results = db.sql('select a._Object_key, a.accID ' + \
		'from #homology h, ACC_Accession a ' + \
		'where a._Object_key = h.mouseKey ' + \
		'and a._MGIType_key = 2 ' + \
		'and a._LogicalDB_key = 1 ' + \
		'and a.prefixPart = "MGI:" ' + \
		'and a.preferred = 1 ', 'auto')
mmgi = createDict(results, '_Object_key', 'accID')

#
# Get the EntrezGene for the human markers
#

results = db.sql('select a._Object_key, a.accID ' + \
                'from #homology h, ACC_Accession a ' + \
                'where a._Object_key = h.humanKey ' + \
                'and a._MGIType_key = 2 ' + \
                'and a._LogicalDB_key = 55', 'auto')
hlocus = createDict(results, '_Object_key', 'accID')

#
# Get the MP Header Terms for the mouse markers
# only select primary ids
#

results = db.sql('select distinct h.mouseKey, a.accID ' + \
                'from #homology h, GXD_AlleleGenotype g, VOC_AnnotHeader v, ACC_Accession a ' + \
                'where h.mouseKey = g._Marker_key ' + \
                'and g._Genotype_key = v._Object_key ' + \
                'and v._AnnotType_key = 1002 ' + \
                'and a._Object_key = v._Term_key ' + \
                'and a._MGIType_key = 13 ' + \
                'and a.preferred = 1 ' + \
                'and a._LogicalDB_key = 34 ', 'auto')
mpheno = createDict(results, 'mouseKey', 'accID')

results = db.sql('select * from #homology order by humanSym', 'auto')

for r in results:
    fp.write(r['humanSym'] + TAB)

    if hlocus.has_key(r['humanKey']):
        fp.write(hlocus[r['humanKey']] + TAB)
    else:
        fp.write(TAB)

    fp.write(r['mouseSym'] + TAB)

    if mmgi.has_key(r['mouseKey']):
        fp.write(mmgi[r['mouseKey']] + TAB)

    if mpheno.has_key(r['mouseKey']):
        fp.write(mpheno[r['mouseKey']] + TAB)
    else:
        fp.write(TAB)

    fp.write(CRT)

reportlib.finish_nonps(fp)