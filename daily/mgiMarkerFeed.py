#!/usr/local/bin/python

'''
#
# Report:
#       Create several tab-delimited (bcp) files for Michael Rosenstein.
#	TR 2794
#
# bcp files
#	1. marker_type.bcp
#	2. species.bcp
#	3. marker.bcp
#	4. marker_label.bcp
#	5. accession_marker.bcp
#	6. allele_type.bcp
#	7. allele_cellline.bcp
#	8. allele_inheritance_mode.bcp
#	9. allele.bcp
#	10. allele_label.bcp
#	11. accession_allele.bcp
#
# Usage:
#       mgiMarkerFeed.py
#
# Notes:
#
# History:
#
# lec	08/14/2001
#	- created
#
'''
 
import sys
import os
import db
import reportlib
import mgi_utils

TAB = reportlib.TAB
CRT = reportlib.CRT

OUTPUTDIR = os.environ['REPORTOUTPUTDIR'] + '/mgimarkerfeed/'

def mgi_status(status):

	if status == "Not Applicable":
		return "NA"

	if status == "Not Specified":
		return "NS"

	if status == "approved":
		return "A"

	if status == "withdrawn":
		return "W"

	if status == "Approved":
		return "A"

	return ''

#
# marker_type
#

fp = open(OUTPUTDIR + 'marker_type.bcp', 'w')

results = db.sql('select _Marker_Type_key, name, ' + \
	'cdate = convert(char(10), creation_date, 101), ' + \
	'mdate = convert(char(10), modification_date, 101) ' + \
	'from MRK_Types', 'auto')
for r in results:
	fp.write(`r['_Marker_Type_key']` + TAB + \
		 r['name'] + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)
fp.close

#
# species
#

fp = open(OUTPUTDIR + 'species.bcp', 'w')

results = db.sql('select _Species_key, name, species, ' + \
	'cdate = convert(char(10), creation_date, 101), ' + \
	'mdate = convert(char(10), modification_date, 101) ' + \
	'from MRK_Species', 'auto')
for r in results:
	fp.write(`r['_Species_key']` + TAB + \
		 r['name'] + TAB + \
		 r['species'] + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)
fp.close

#
# allele_type
#

fp = open(OUTPUTDIR + 'allele_type.bcp', 'w')

results = db.sql('select _Allele_Type_key, alleleType, ' + \
	'cdate = convert(char(10), creation_date, 101), ' + \
	'mdate = convert(char(10), modification_date, 101) ' + \
	'from ALL_Type', 'auto')
for r in results:
	fp.write(`r['_Allele_Type_key']` + TAB + \
		 r['alleleType'] + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)
fp.close

#
# allele_cellline
#

fp = open(OUTPUTDIR + 'allele_cellline.bcp', 'w')

results = db.sql('select _CellLine_key, cellLine, cellLineStrain, ' + \
	'cdate = convert(char(10), creation_date, 101), ' + \
	'mdate = convert(char(10), modification_date, 101) ' + \
	'from ALL_CellLine_View', 'auto')
for r in results:
	fp.write(`r['_CellLine_key']` + TAB + \
		 r['cellLine'] + TAB + \
		 r['cellLineStrain'] + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)
fp.close

#
# allele_inheritance_mode
#

fp = open(OUTPUTDIR + 'allele_inheritance_mode.bcp', 'w')

results = db.sql('select _Mode_key, mode, ' + \
	'cdate = convert(char(10), creation_date, 101), ' + \
	'mdate = convert(char(10), modification_date, 101) ' + \
	'from ALL_Inheritance_Mode', 'auto')
for r in results:
	fp.write(`r['_Mode_key']` + TAB + \
		 r['mode'] + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)
fp.close

#
# marker
# marker_label
# accession_marker
#

fp1 = open(OUTPUTDIR + 'marker.bcp', 'w')
fp2 = open(OUTPUTDIR + 'marker_label.bcp', 'w')
fp3 = open(OUTPUTDIR + 'accession_marker.bcp', 'w')

cmds = []

cmds.append('select m._Marker_key into #markers ' + \
	'from MRK_Marker m ' + \
	'where m._Species_key = 1 ' + \
	'and exists (select 1 from ACC_Accession a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1) ')

cmds.append('select m._Marker_key, m._Species_key, m._Marker_Type_key, s.status, ' + \
	'm.symbol, m.name, m.chromosome, m.cytogeneticOffset, ' + \
	'cdate = convert(char(10), m.creation_date, 101), ' + \
	'mdate = convert(char(10), m.modification_date, 101) ' + \
	'from #markers k, MRK_Marker m, MRK_Status s ' + \
	'where k._Marker_key = m._Marker_key ' + \
	'and m._Marker_Status_key = s._Marker_Status_key ')

cmds.append('select m._Marker_key, m.label, m.labelType, s.status, ' + \
	'cdate = convert(char(10), m.creation_date, 101), ' + \
	'mdate = convert(char(10), m.modification_date, 101) ' + \
	'from #markers k, MRK_Label m, MRK_Status s ' + \
	'where k._Marker_key = m._Marker_key ' + \
	'and m._Marker_Status_key = s._Marker_Status_key')

cmds.append('select m.accID, m.LogicalDB, m._Object_key, m.preferred, ' + \
	'cdate = convert(char(10), m.creation_date, 101), ' + \
	'mdate = convert(char(10), m.modification_date, 101) ' + \
	'from #markers k, MRK_Acc_View m ' + \
	'where k._Marker_key = m._Object_key ' + \
	'and m.prefixPart = "MGI:"')

results = db.sql(cmds, 'auto')

for r in results[1]:
	fp1.write(`r['_Marker_key']` + TAB + \
	         `r['_Species_key']` + TAB + \
	         `r['_Marker_Type_key']` + TAB + \
		 mgi_status(r['status']) + TAB + \
		 r['symbol'] + TAB + \
		 r['name'] + TAB + \
		 r['chromosome'] + TAB + \
		 mgi_utils.prvalue(r['cytogeneticOffset']) + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[2]:
	fp2.write(`r['_Marker_key']` + TAB + \
		 r['label'] + TAB + \
		 r['labelType'] + TAB + \
		 mgi_status(r['status']) + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[3]:
	fp3.write(r['accID'] + TAB + \
		 r['LogicalDB'] + TAB + \
		 `r['_Object_key']` + TAB + \
		 `r['preferred']` + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

fp1.close
fp2.close
fp3.close

#
# allele
# allele_label
# accession_allele
#

fp1 = open(OUTPUTDIR + 'allele.bcp', 'w')
fp2 = open(OUTPUTDIR + 'allele_label.bcp', 'w')
fp3 = open(OUTPUTDIR + 'accession_allele.bcp', 'w')

cmds = []

cmds.append('select m._Allele_key into #alleles ' + \
	'from ALL_Allele m ' + \
	'where m._Allele_Status_key = 4 ')

cmds.append('select m._Allele_key, m._Marker_key, m._Mode_key, m._Allele_Type_key, m._CellLine_key, ' + \
	'm.status, m.strain, m.symbol, m.name, ' + \
	'cdate = convert(char(10), m.creation_date, 101), ' + \
	'mdate = convert(char(10), m.modification_date, 101) ' + \
	'from #alleles a, ALL_Allele_View m ' + \
	'where a._Allele_key = m._Allele_key ')

cmds.append('select m._Allele_key, m.name, labelType = "N", m.status, ' + \
	'cdate = convert(char(10), m.creation_date, 101), ' + \
	'mdate = convert(char(10), m.modification_date, 101) ' + \
	'from #alleles a, ALL_Allele_View m ' + \
	'where a._Allele_key = m._Allele_key ' + \
	'union ' + \
	'select m._Allele_key, m.symbol, labelType = "S", m.status, ' + \
	'cdate = convert(char(10), m.creation_date, 101), ' + \
	'mdate = convert(char(10), m.modification_date, 101) ' + \
	'from #alleles a, ALL_Allele_View m ' + \
	'where a._Allele_key = m._Allele_key ')

cmds.append('select m.accID, m.LogicalDB, m._Object_key, m.preferred, ' + \
	'cdate = convert(char(10), m.creation_date, 101), ' + \
	'mdate = convert(char(10), m.modification_date, 101) ' + \
	'from #alleles k, ALL_Acc_View m ' + \
	'where k._Allele_key = m._Object_key ' + \
	'and m.prefixPart = "MGI:"')

results = db.sql(cmds, 'auto')

for r in results[1]:
	fp1.write(`r['_Allele_key']` + TAB + \
	         `r['_Marker_key']` + TAB + \
	         `r['_Mode_key']` + TAB + \
	         `r['_Allele_Type_key']` + TAB + \
	         `r['_CellLine_key']` + TAB + \
		 mgi_status(r['status']) + TAB + \
		 r['strain'] + TAB + \
		 r['symbol'] + TAB + \
		 r['name'] + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[2]:
	fp2.write(`r['_Allele_key']` + TAB + \
		 r['name'] + TAB + \
		 r['labelType'] + TAB + \
		 mgi_status(r['status']) + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[3]:
	fp3.write(r['accID'] + TAB + \
		 r['LogicalDB'] + TAB + \
		 `r['_Object_key']` + TAB + \
		 `r['preferred']` + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

fp1.close
fp2.close
fp3.close

#
# nightly_reports.sh will compress and tar up the files
#
