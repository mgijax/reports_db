#!/usr/local/bin/python

'''
#
# MGI_GeneOMIM.py
#
# Report:
#
#	field 1: OMIM ID
#	field 2: MGI Marker ID
#
# Usage:
#       MGI_GeneOMIM.py
#
# History:
#
# lec	11/03/2016
#	TR12427/Disease Ontology (DO)
#
# sc    06/27/2014
#       - TR11560 Feature Relationships project
#       - exclude feature type 'mutation defined region'
#
# lec	03/01/2011
#	- TR10616/for Donna Maglott/NCBI
#
'''
 
import sys 
import os
import reportlib
import db

db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE(False)

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

#
# excluded markers
# feature type='mutation defined region', key=11928467
#
results = db.sql('''
        select _marker_key
        from MRK_MCV_Cache
        where qualifier = 'D'
        and _MCVTerm_key = 11928467
        ''', 'auto')
excludedMarkerList = []
for r in results:
        markerKey = r['_marker_key']
        excludedMarkerList.append(markerKey)

#
# Genotype/Markers that have OMIM annotations
#

results = db.sql('''select distinct a.accID omimID, ma.accID, ag._Marker_key
        from GXD_AlleleGenotype ag, VOC_Annot va, ACC_Accession a, ACC_Accession ma
        where ag._Genotype_key = va._Object_key
        and va._AnnotType_key = 1005
        and va._Term_key = a._Object_key
        and a._LogicalDB_key = 15
        and a.preferred = 1
        and ag._Marker_key = ma._Object_key 
        and ma._MGIType_key = 2
        and ma._LogicalDB_key = 1
        and ma.prefixPart = 'MGI:'
        and ma.preferred = 1
	order by ma.accID
	''', 'auto')

fp = reportlib.init('MGI_GeneOMIM', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

for r in results:
    if r['_Marker_key'] in excludedMarkerList:
	continue
    fp.write(r['accID'] + TAB)
    fp.write(r['omimID'] + CRT)

#
# only if someone asks us for it
# Genotype/Markers that have DO annotations
#

#results = db.sql('''select distinct a.accID doID, ma.accID, ag._Marker_key
#        from GXD_AlleleGenotype ag, VOC_Annot va, ACC_Accession a, ACC_Accession ma
#        where ag._Genotype_key = va._Object_key
#        and va._AnnotType_key = 1020
#        and va._Term_key = a._Object_key
#        and a._LogicalDB_key = 191
#        and a.preferred = 1
#        and ag._Marker_key = ma._Object_key 
#        and ma._MGIType_key = 2
#        and ma._LogicalDB_key = 1
#        and ma.prefixPart = 'MGI:'
#        and ma.preferred = 1
#	order by ma.accID
#	''', 'auto')
#
#fp2 = reportlib.init('MGI_GeneDO', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
#for r in results:
#    if r['_Marker_key'] in excludedMarkerList:
#	continue
#    fp2.write(r['accID'] + TAB)
#    fp2.write(r['doID'] + CRT)

reportlib.finish_nonps(fp)	# non-postscript file
#reportlib.finish_nonps(fp2)	# non-postscript file

