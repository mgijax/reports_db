#!/usr/local/bin/python

'''
#
# TR 4533/TR 5485
#
# Report:
#       Produce a tab-delimited report of PhenoSlim terms and associated
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
#       MRK_Ensembl_PS.py
#
# Notes:
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

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []

cmds.append(
    'select a1.accID "MGI", m.symbol, m.name, ' + \
           'a2.accID "Ensembl", a._Annot_key, a._Term_key ' + \
    'into #mgi_ensembl ' + \
    'from ACC_Accession a1, ACC_Accession a2, MRK_Marker m, ' + \
         'GXD_AlleleGenotype g, VOC_Annot a ' + \
    'where a1._Object_key = a2._Object_key and ' + \
          'a1._Object_key = m._Marker_key and ' + \
          'm._Marker_key = g._Marker_key and ' + \
          'g._Genotype_key = a._Object_key and ' + \
          'a1._LogicalDB_key = 1 and ' + \
          'a1._MGIType_key = 2 and ' + \
          'a1.prefixPart = "MGI:" and ' + \
          'a1.preferred = 1 and ' + \
          'a2._LogicalDB_key = 33 and ' + \
          'a2._MGIType_key = 2 and ' + \
          'a2.preferred = 1 and ' + \
          'm._Marker_Type_key = 1 and ' + \
          'a._AnnotType_key = 1002')

cmds.append(
    'select distinct m.MGI, m.symbol, m.name, m.Ensembl, ' + \
           'a1.accID "classification", t.term, a2.accID "pubmed" ' + \
    'from #mgi_ensembl m, VOC_Term t, VOC_Evidence e, ' + \
         'ACC_Accession a1, ACC_Accession a2 ' + \
    'where m._Term_key = t._Term_key and ' + \
          't._Term_key = a1._Object_key and ' + \
          'a1._MGIType_key = 13 and ' + \
          'a1.preferred = 1 and ' + \
          'm._Annot_key = e._Annot_key and ' + \
          'e._Refs_key = a2._Object_key and ' + \
          'a2._LogicalDB_key = 29 and ' + \
          'a2.preferred = 1 ' + \
    'order by m.symbol')

results = db.sql(cmds, 'auto')

for r in results[1]:
    fp.write(r['MGI'] + TAB + r['symbol'] + TAB + r['name'] + TAB +
             r['Ensembl'] + TAB + r['classification'] + TAB +
             r['term'] + TAB + r['pubmed'] + CRT)

reportlib.finish_nonps(fp)

