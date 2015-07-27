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
# create a tab-delimitted report with the following columns:
# 
# MGI ID 
# Gene Symbol
# Source - For gene models, use the provider; use "MGI" for the MGI Marker Type.
# Gene ID -  For gene models, use the accID of the gene model sequence, 
#	     for MGI Markers, use the MGI ID of the Marker. 
#            The IDs should link out to respective gene pages, as in the mockup.
# BioType - For gene models, use the Raw BioType values, for MGI, 
#	    use the MGI Feature Type.
# MGI_Rep_Gene_Model - enter "Representative" if gene model sequence 
#	               is representative genomic, else, leave blank
# Order by MGI ID
# 
# Example:
# 
# MGI ID	Symbol   Source    Gene ID            BioType	                MGI_Rep_Gene_Model
# MGI:1913301	0610012G03Rik   MGI     MGI:1913301     pseudogene
# MGI:1913301	610012G03Rik   Ensembl Gene Model      ENSMUSG00000047112      pseudogene      Representative
# MGI:1913301	0610012G03Rik   NCBI Gene Model 106264  miscRNA
#
# History:
#
# sc	04/26/2011
#	- TR10336 - add new col 1 MGI ID, sort by MGI ID
#	  note for the gene row, MGI ID is repeated in column 4
#
# sc	11-15-2010
#	- TR10308 - update the biotype conflict algorithm
#
# lec	09/14/2010
#	- TR 10336
#	  include date of lthe last update
#	  include a count of the number of genes (NOT the number of rows)
#
# lec	02/24/2010
#	- created
#
'''
 
import sys 
import os
import mgi_utils
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

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('#\n')
fp.write('# This report lists the Markers that have a Biotype Conflict.\n')
fp.write('#\n')
fp.write('#  date report was generated:  %s\n#\n' % (mgi_utils.date()))
fp.write('#  column 1: MGI ID\n')
fp.write('#  column 2: Symbol\n')
fp.write('#  column 3: Source:  for gene models, this is the provider\n')
fp.write('#                     for markers, this is "MGI"\n')
fp.write('#  column 4: Gene ID: for gene models, this is the sequence ID\n')
fp.write('#                     for markers, this is the MGI ID\n')
fp.write("#  column 5: Biotype: for gene models, this is the provider's biotype value\n")
fp.write('#                     for markers, this is the MGI feature type\n')
fp.write('#  column 6: MGI_Rep_Gene_Model: "Representative" if gene model sequence is the representative genomic\n')
fp.write('#\n\n')

results = db.sql('''
		 select s._Marker_key, s._Qualifier_key, s.accID, s.rawbiotype, 
			m.symbol, a.accID as mgiID, mcv.term as featureType,
			v.term as provider
		 from SEQ_Marker_Cache s, MRK_Marker m, ACC_Accession a, 
			MRK_MCV_Cache mcv, VOC_Term v
		 where s._BiotypeConflict_key = 5420767
		 and s._LogicalDB_key in (59, 60, 85)
		 and s._Marker_key = m._Marker_key
		 and s._Marker_key = a._Object_key
		 and a._MGIType_key = 2
		 and a.prefixPart = 'MGI:'
		 and a._LogicalDB_key = 1
		 and a.preferred = 1
		 and s._Marker_key = mcv._Marker_key
		 and mcv.qualifier = 'D'
		 and s._SequenceProvider_key = v._Term_key
		 order by a.accid
		 ''', 'auto')

markerList = {}

for r in results:

    key = r['mgiID']
    value = r['symbol']

    # print one row of the marker record

    if not markerList.has_key(key):
	fp.write(r['mgiID'] + TAB)
        fp.write(r['symbol'] + TAB)
	fp.write('MGI' + TAB)
	fp.write(r['mgiID'] + TAB)
	fp.write(r['featureType'] + TAB)
	fp.write(CRT)
	markerList[key] = value

    # print row of the gene model sequence
    fp.write(r['mgiID'] + TAB)
    fp.write(r['symbol'] + TAB)
    fp.write(r['provider'] + TAB)
    fp.write(r['accID'] + TAB)
    fp.write(mgi_utils.prvalue(r['rawbiotype']) + TAB)

    if r['_Qualifier_key'] == 615419:
	fp.write('Representative')
    fp.write(CRT)

fp.write(CRT + '(%d genes affected)' % (len(markerList)) + CRT)

db.useOneConnection(0)
reportlib.finish_nonps(fp)	# non-postscript file

