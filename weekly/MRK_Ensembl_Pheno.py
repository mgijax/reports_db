#!/usr/local/bin/python

'''
#
# TR 4533/TR 5485
#
# Report:
#       Produce a tab-delimited report of MP terms and associated
#       PubMed IDs for MGI markers that have associations to ENSEMBL records.
#
#       Display the following fields:
#
#           MGI ID
#           Symbol
#           Name
#           ENSEMBL ID
#           Classification ID
#           Classification Term
#           PubMed ID of the reference for the term/marker association
#
#       One row per MGI ID/Classification ID pair
#
#       Sort by MGI ID
#
# Usage:
#       MRK_Ensembl_Pheno.py
#
# History:
#
# lec	01/20/2004
#	- made a public report
#
# dbm   2/24/2003
#       - created
#
'''

import sys
import string
import os
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

#
# all markers that have an Ensembl ID
#

db.sql('select a1.accID, m._Marker_key, m.symbol, m.name ' + \
       'into #markers ' + \
       'from MRK_Marker m, ACC_Accession a1 ' + \
       'where m._Organism_key = 1 ' + \
       'and m._Marker_Type_key = 1 ' + \
       'and m._Marker_key = a1._Object_key ' + \
       'and a1._MGIType_key = 2 ' + \
       'and a1._LogicalDB_key = 1 ' + \
       'and a1.prefixPart = "MGI:" ' + \
       'and a1.preferred = 1 ' + \
       'and exists (select 1 from ACC_Accession a2 where m._Marker_key = a2._Object_key ' + \
       'and a2._MGIType_key = 2 ' + \
       'and a2._LogicalDB_key = 60 ' + \
       'and a2.preferred = 1)', None)

db.sql('create index idx1 on #markers(_Marker_key)', None)
db.sql('create index idx3 on #markers(symbol)', None)

#
# ensembl ids
#

results = db.sql('select m._Marker_key, a.accID from #markers m, ACC_Accession a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 60', 'auto')
ensembl = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if not ensembl.has_key(key):
	ensembl[key] = []
    ensembl[key].append(value)

#
# header terms for markers->genotype->MP annotations
#

results = db.sql('select distinct m._Marker_key, termID = a.accID, t.term ' + \
                 'from #markers m, GXD_AlleleGenotype g, VOC_AnnotHeader h, VOC_Term t, ACC_Accession a ' + \
                 'where m._Marker_key = g._Marker_key ' + \
		 'and g._Genotype_key = h._Object_key ' + \
		 'and h._AnnotType_key = 1002 ' + \
		 'and h._Term_key = t._Term_key ' + \
                 'and t._Term_key = a._Object_key ' + \
                 'and a._MGIType_key = 13 ' + \
                 'and a.preferred = 1 ', 'auto')
headerTerm = {}
for r in results:
    key = r['_Marker_key']
    value = r
    if not headerTerm.has_key(key):
	headerTerm[key] = []
    headerTerm[key].append(value)

results = db.sql('select * from #markers order by symbol', 'auto')

#
# print one line per marker/MP header term
#

for r in results:

    key = r['_Marker_key']

    if not headerTerm.has_key(key):
	continue

    for a in headerTerm[key]:
        fp.write(r['accID'] + TAB + \
	         r['symbol'] + TAB + \
	         r['name'] + TAB + \
                 string.join(ensembl[key], ',') + TAB + \
	         a['termID'] + TAB + \
	         a['term'] + CRT)

reportlib.finish_nonps(fp)

