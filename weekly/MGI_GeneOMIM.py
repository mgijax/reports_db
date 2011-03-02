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
# lec	03/01/2011
#	- TR10616/for Donna Maglott/NCBI
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

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# Genotype/Markers that have OMIM annotations
#

results = db.sql('''select distinct omimID = a.accID, ma.accID
        from GXD_AlleleGenotype ag, VOC_Annot va, ACC_Accession a, ACC_Accession ma
        where ag._Genotype_key = va._Object_key
        and va._AnnotType_key = 1005
        and va._Term_key = a._Object_key
        and a._LogicalDB_key = 15
        and a.preferred = 1
        and ag._Marker_key = ma._Object_key 
        and ma._MGIType_key = 2
        and ma._LogicalDB_key = 1
        and ma.prefixPart = "MGI:"
        and ma.preferred = 1
	order by ma.accID
	''', 'auto')

for r in results:

    fp.write(r['accID'] + TAB)
    fp.write(r['omimID'] + CRT)

reportlib.finish_nonps(fp)	# non-postscript file

