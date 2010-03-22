#!/usr/local/bin/python

'''
#
# Report:
#       TR9239/BioType Conflict Report
#
# See TR directory, requirement #8 (below):
#
# A public report of BioType conflicts should be produced in 
# conjunction with each run of the BioType mismatch detection process.
#
# Report Details:
# For Markers that have a Biotype Conflict, 
# create a tab-delimmited report with the following columns:
# 
# Gene Symbol
# Source	 For gene models, use the provider; use "MGI" for the MGI Marker Type.
# Gene ID	 For gene models, use the accID of the gene model sequence, 
#	         for MGI Markers, use the MGI ID of the Marker. 
#                The IDs should link out to respective gene pages, as in the mockup.
# BioType	 For gene models, use the Raw BioType values, for MGI, use the MGI Marker Type.
# MGI_Rep_Gene_Model enter "Representative" if gene model sequence is representative genomic, 
#                    else, leave blank
# Order by Gene Symbol
# 
# Example:
# 
# Symbol   Source    Gene ID            BioType	                MGI_Rep_Gene_Model
# foo  	   MGI	     MGI:97511	        Gene	
# foo  	   VEGA	     OTTMUSG00000026053	Known protein coding	Representative
# foo  	   Ensembl   ENSMUSG00000021587	Known protein coding	
# foo  	   Ensembl   ENSMUSG00000027419	Known protein coding	
# foo  	   NCBI	     18548	        Pseudo	
#
# History:
#
# lec	02/24/2010
#	- created
#
'''
 
import sys 
import os
import db
import reportlib
import mgi_utils

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('#\n')
fp.write('# This report lists the Markers that have a Biotype Conflict that involves a pseudogene.\n')
fp.write('#\n')
fp.write('#  column 1: Symbol\n')
fp.write('#  column 2: Source:  for gene models, this is the provider\n')
fp.write('#                     for markers, this is "MGI"\n')
fp.write('#  column 3: Gene ID: for gene models, this is the sequence ID\n')
fp.write('#                     for markers, this is the MGI ID\n')
fp.write("#  column 4: Biotype: for gene models, this is the provider's biotype value\n")
fp.write('#                     for markers, this is the MGI marker type\n')
fp.write('#  column 5: MGI_Rep_Gene_Model: "Representative" if gene model sequence is the representative genomic\n')
fp.write('#\n\n')

results = db.sql('''
		 select s._Marker_key, s._Qualifier_key, s.accID, s.rawbiotype, 
			m.symbol, mgiID = a.accID, markerType = t.name,
			provider = v.term
		 from SEQ_Marker_Cache s, MRK_Marker m, ACC_Accession a, MRK_Types t, VOC_Term v
		 where s._BiotypeConflict_key = 5420767
		 and s._LogicalDB_key in (59, 60, 85)
		 and s._Marker_key = m._Marker_key
		 and s._Marker_key = a._Object_key
		 and a._MGIType_key = 2
		 and a.prefixPart = "MGI:"
		 and a._LogicalDB_key = 1
		 and a.preferred = 1
		 and s._Marker_Type_key = t._Marker_Type_key
		 and s._SequenceProvider_key = v._Term_key
		 order by m.symbol, m._Marker_key
		 ''', 'auto')

markerList = {}

for r in results:

    key = r['_Marker_key']
    value = r['symbol']

    # print one row of the marker record

    if not markerList.has_key(key):
        fp.write(r['symbol'] + TAB)
	fp.write('MGI' + TAB)
	fp.write(r['mgiID'] + TAB)
	fp.write(r['markerType'] + TAB)
	fp.write(CRT)
	markerList[key] = value

    # print row of the gene model sequence

    fp.write(r['symbol'] + TAB)
    fp.write(r['provider'] + TAB)
    fp.write(r['accID'] + TAB)
    fp.write(mgi_utils.prvalue(r['rawbiotype']) + TAB)

    if r['_Qualifier_key'] == 615419:
	fp.write('Representative')
    fp.write(CRT)

db.useOneConnection(0)
reportlib.finish_nonps(fp)	# non-postscript file

