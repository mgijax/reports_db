#!/usr/local/bin/python

'''
#
# Report:
#       Create several tab-delimited (bcp) files for Michael Rosenstein.
#	TR 2794
#	TR 3137 - add MRK_Offset
#	TR 3345 - add Strain info; incorporate marker offset into marker.bcp
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
#	12. strain_marker.bcp
#	13. strain_synonym.bcp
#	14. strain.bcp
#	15. strain_type.bcp
#	16. strain_strain_type.bcp
#	17. accession_strain.bcp
#	18. strain_species.bcp
#
# Usage:
#       mgiMarkerFeed.py
#
# Notes:
#
# History:
#
# lec	02/20/2003
#	- TR 1892; MRK_Label
#
# lec	04/09/2002
#	- removed needsReview and standard
#
# lec	08/14/2001
#	- created
#
'''
 
import sys
import os
import string
import db
import reportlib
import mgi_utils

TAB = reportlib.TAB
CRT = reportlib.CRT

OUTPUTDIR = os.environ['REPORTOUTPUTDIR'] + '/mgimarkerfeed/'

def strip_newline(s):

	if string.find (s, '\\n') > -1:
		s = string.join(string.split(s, '\\n'), '')
		return s
	elif string.find(s, '\012') > -1:
		s = string.join(string.split(s, '\012'), '')
		return s
	else:
		return s

def mgi_status(status):

	if status == "Approved":
		return "A"

	if status == "official":
		return "O"

	if status == "interim":
		return "I"

	if status == "withdrawn":
		return "W"

	return ''

#
# marker_type
#

fp = open(OUTPUTDIR + 'marker_type.bcp', 'w')

results = db.sql('select _Marker_Type_key, name, ' + \
	'cdate = convert(char(20), creation_date, 100), ' + \
	'mdate = convert(char(20), modification_date, 100) ' + \
	'from MRK_Types', 'auto', execute = 1)
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

results = db.sql('select _Organism_key, commonName, latinName, ' + \
	'cdate = convert(char(20), creation_date, 100), ' + \
	'mdate = convert(char(20), modification_date, 100) ' + \
	'from MGI_Organism_Marker_View', 'auto', execute = 1)
for r in results:
	fp.write(`r['_Organism_key']` + TAB + \
		 r['commonName'] + TAB + \
		 r['latinName'] + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)
fp.close

#
# allele_type
#

fp = open(OUTPUTDIR + 'allele_type.bcp', 'w')

results = db.sql('select _Allele_Type_key, alleleType, ' + \
	'cdate = convert(char(20), creation_date, 100), ' + \
	'mdate = convert(char(20), modification_date, 100) ' + \
	'from ALL_Type', 'auto', execute = 1)
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

results = db.sql('select _CellLine_key, cellLine, _Strain_key, ' + \
	'cdate = convert(char(20), creation_date, 100), ' + \
	'mdate = convert(char(20), modification_date, 100) ' + \
	'from ALL_CellLine', 'auto', execute = 1)
for r in results:
	fp.write(`r['_CellLine_key']` + TAB + \
		 r['cellLine'] + TAB + \
		 `r['_Strain_key']` + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)
fp.close

#
# allele_inheritance_mode
#

fp = open(OUTPUTDIR + 'allele_inheritance_mode.bcp', 'w')

results = db.sql('select _Mode_key, mode, ' + \
	'cdate = convert(char(20), creation_date, 100), ' + \
	'mdate = convert(char(20), modification_date, 100) ' + \
	'from ALL_Inheritance_Mode', 'auto', execute = 1)
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

#
# select all mouse markers which have a preferred MGI Accession ID.
# this will include any splits (since the MGI Acc ID stays with 
# the split symbol).
#

cmds.append('select m._Marker_key, m.symbol, m.name into #markers ' + \
	'from MRK_Marker m ' + \
	'where m._Organism_key = 1 ' + \
	'and exists (select 1 from ACC_Accession a ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.preferred = 1) ')

#
# select data fields for marker.bcp
#

cmds.append('select m._Marker_key, m._Organism_key, m._Marker_Type_key, s.status, ' + \
	'm.symbol, m.name, m.chromosome, m.cytogeneticOffset, o.offset, ' + \
	'cdate = convert(char(20), m.creation_date, 100), ' + \
	'mdate = convert(char(20), m.modification_date, 100) ' + \
	'from #markers k, MRK_Marker m, MRK_Status s, MRK_Offset o ' + \
	'where k._Marker_key = m._Marker_key ' + \
	'and m._Marker_Status_key = s._Marker_Status_key ' + \
	'and k._Marker_key = o._Marker_key ' + \
	'and o.source = 0')

#
# select data fields for marker_label.bcp
# status = 1 if the label matches the MRK_Marker.symbol or MRK_Marker.name
# field, else status = 0
#

cmds.append('select distinct m._Marker_key, m.label, m.labelType, status = 1, ' + \
	'cdate = convert(char(20), m.creation_date, 100), ' + \
	'mdate = convert(char(20), m.modification_date, 100) ' + \
	'from #markers k, MRK_Label m ' + \
	'where k._Marker_key = m._Marker_key ' + \
	'and m.labelType = "MS" ' + \
	'and k.symbol = m.label ')

cmds.append('select distinct m._Marker_key, m.label, m.labelType, status = 1, ' + \
	'cdate = convert(char(20), m.creation_date, 100), ' + \
	'mdate = convert(char(20), m.modification_date, 100) ' + \
	'from #markers k, MRK_Label m ' + \
	'where k._Marker_key = m._Marker_key ' + \
	'and m.labelType = "MN" ' + \
	'and k.name = m.label')

cmds.append('select distinct m._Marker_key, m.label, m.labelType, status = 0, ' + \
	'cdate = convert(char(20), m.creation_date, 100), ' + \
	'mdate = convert(char(20), m.modification_date, 100) ' + \
	'from #markers k, MRK_Label m ' + \
	'where k._Marker_key = m._Marker_key ' + \
	'and m.labelType in ("MS", "MN", "MY") ' + \
	'and k.symbol != m.label and k.name != m.label ')

#
# select data fields for accession_marker.bcp
# select MGI Accession IDs only
#

cmds.append('select m.accID, LogicalDB = l.name, m._Object_key, m.preferred, ' + \
	'cdate = convert(char(20), m.creation_date, 100), ' + \
	'mdate = convert(char(20), m.modification_date, 100) ' + \
	'from #markers k, ACC_Accession m, ACC_LogicalDB l ' + \
	'where k._Marker_key = m._Object_key ' + \
	'and m._MGIType_key = 2 ' + \
	'and m.prefixPart = "MGI:" ' + \
	'and m._LogicalDB_key = 1 ' + \
	'and m._LogicalDB_key = l._LogicalDB_key')

results = db.sql(cmds, 'auto', execute = 1)

for r in results[1]:
	fp1.write(`r['_Marker_key']` + TAB + \
	         `r['_Organism_key']` + TAB + \
	         `r['_Marker_Type_key']` + TAB + \
		 mgi_status(r['status']) + TAB + \
		 strip_newline(r['symbol']) + TAB + \
		 strip_newline(r['name']) + TAB + \
		 r['chromosome'] + TAB + \
		 mgi_utils.prvalue(r['cytogeneticOffset']) + TAB + \
		 `r['offset']` + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[2]:
	fp2.write(`r['_Marker_key']` + TAB + \
		 strip_newline(r['label']) + TAB + \
		 r['labelType'] + TAB + \
		 `r['status']` + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[3]:
	fp2.write(`r['_Marker_key']` + TAB + \
		 strip_newline(r['label']) + TAB + \
		 r['labelType'] + TAB + \
		 `r['status']` + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[4]:
	fp2.write(`r['_Marker_key']` + TAB + \
		 strip_newline(r['label']) + TAB + \
		 r['labelType'] + TAB + \
		 `r['status']` + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[5]:
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

#
# select all alleles with a status of "approved"
# all other statuses are private/confidential alleles
#

cmds.append('select m._Allele_key into #alleles ' + \
	'from ALL_Allele m ' + \
	'where m._Allele_Status_key = 4 ')

#
# select data fields for allele.bcp
#

cmds.append('select m._Allele_key, m._Marker_key, m._Mode_key, m._Allele_Type_key, m._CellLine_key, ' + \
	'm.status, m._Strain_key, m.symbol, m.name, ' + \
	'cdate = convert(char(20), m.creation_date, 100), ' + \
	'mdate = convert(char(20), m.modification_date, 100) ' + \
	'from #alleles a, ALL_Allele_View m ' + \
	'where a._Allele_key = m._Allele_key ')

#
# select data fields for allele_label.bcp
# this includes all allele symbols, names and synonyms
# status = 1 if the label matches the ALL_Allele.symbol or ALL_Allele.name
# field, else status = 0
#

cmds.append('select m._Allele_key, m.name, labelType = "AN", status = 1, ' + \
	'cdate = convert(char(20), m.creation_date, 100), ' + \
	'mdate = convert(char(20), m.modification_date, 100) ' + \
	'from #alleles a, ALL_Allele m ' + \
	'where a._Allele_key = m._Allele_key ' + \
	'union ' + \
	'select m._Allele_key, m.symbol, labelType = "AS", status = 1, ' + \
	'cdate = convert(char(20), m.creation_date, 100), ' + \
	'mdate = convert(char(20), m.modification_date, 100) ' + \
	'from #alleles a, ALL_Allele m ' + \
	'where a._Allele_key = m._Allele_key ' + \
	'union ' + \
	'select m._Allele_key, s.synonym, labelType = "AY", status = 0, ' + \
	'cdate = convert(char(20), m.creation_date, 100), ' + \
	'mdate = convert(char(20), m.modification_date, 100) ' + \
	'from #alleles a, ALL_Allele m, ALL_Synonym s ' + \
	'where a._Allele_key = m._Allele_key ' + \
	'and a._Allele_key = s._Allele_key ')

#
# select data fields for accession_allele.bcp
#

cmds.append('select m.accID, m.LogicalDB, m._Object_key, m.preferred, ' + \
	'cdate = convert(char(20), m.creation_date, 100), ' + \
	'mdate = convert(char(20), m.modification_date, 100) ' + \
	'from #alleles k, ACC_Accession m ' + \
	'where k._Allele_key = m._Object_key ' + \
	'and m._MGIType_key = 11 ' + \
	'and m.prefixPart = "MGI:" ' + \
	'and m._LogicalDB_key = 1')

results = db.sql(cmds, 'auto', execute = 1)

for r in results[1]:
	fp1.write(`r['_Allele_key']` + TAB + \
	         `r['_Marker_key']` + TAB + \
	         `r['_Mode_key']` + TAB + \
	         `r['_Allele_Type_key']` + TAB + \
	         `r['_CellLine_key']` + TAB + \
		 mgi_status(r['status']) + TAB + \
		 `r['_Strain_key']` + TAB + \
		 strip_newline(r['symbol']) + TAB + \
		 strip_newline(r['name']) + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[2]:
	fp2.write(`r['_Allele_key']` + TAB + \
		 strip_newline(r['name']) + TAB + \
		 r['labelType'] + TAB + \
		 `r['status']` + TAB + \
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
# strain.bcp
# strain_marker.bcp
# strain_synonym.bcp
# strain_type.bcp
# strain_strain_type.bcp
# accession_strain.bcp
# strain_species.bcp
#

fp1 = open(OUTPUTDIR + 'strain.bcp', 'w')
fp2 = open(OUTPUTDIR + 'strain_marker.bcp', 'w')
fp3 = open(OUTPUTDIR + 'strain_synonym.bcp', 'w')
fp4 = open(OUTPUTDIR + 'strain_type.bcp', 'w')
fp5 = open(OUTPUTDIR + 'strain_strain_type.bcp', 'w')
fp6 = open(OUTPUTDIR + 'accession_strain.bcp', 'w')
fp7 = open(OUTPUTDIR + 'strain_species.bcp', 'w')

#
# select all strains which have a Jax Registry ID or MMRRC ID
# plus all strains which are cross-referenced by Allele or Allele CellLine
#

cmds = []
cmds.append('select distinct s._Strain_key, m._Species_key, s.strain, s.private, ' + \
      'cdate = convert(char(20), s.creation_date, 100), ' + \
      'mdate = convert(char(20), s.modification_date, 100) ' + \
      'into #strains ' + \
      'from PRB_Strain s, MLP_Strain m, ACC_Accession a ' + \
      'where s._Strain_key = m._Strain_key ' + \
      'and s._Strain_key = a._Object_key ' + \
      'and a._MGIType_key = 10 ' + \
      'and a._LogicalDB_key in (22, 38) ')

#      'union ' + \
#      'select s._Strain_key, m._Species_key, s.strain, s.standard, s.needsReview, s.private, ' + \
#      'cdate = convert(char(20), s.creation_date, 100), ' + \
#      'mdate = convert(char(20), s.modification_date, 100) ' + \
#      'from PRB_Strain s, MLP_Strain m, All_Allele a ' + \
#      'where s._Strain_key = m._Strain_key ' + \
#      'and s._Strain_key = a._Strain_key ' + \
#      'union ' + \
#      'select s._Strain_key, m._Species_key, s.strain, s.standard, s.needsReview, s.private, ' + \
#      'cdate = convert(char(20), s.creation_date, 100), ' + \
#      'mdate = convert(char(20), s.modification_date, 100) ' + \
#      'from PRB_Strain s, MLP_Strain m, All_CellLine a ' + \
#      'where s._Strain_key = m._Strain_key ' + \
#      'and s._Strain_key = a._Strain_key')

cmds.append('create unique index index_strain on #strains(_Strain_key)')

cmds.append('select * from #strains')

cmds.append('select distinct m._Strain_key, m._Marker_key, m._Allele_key, s.private, ' + \
      'cdate = convert(char(20), m.creation_date, 100), ' + \
      'mdate = convert(char(20), m.modification_date, 100) ' + \
      'from #strains s, PRB_Strain_Marker m ' + \
      'where s._Strain_key = m._Strain_key')

cmds.append('select m._Synonym_key, m._Strain_key, m.synonym, s.private, ' + \
      'cdate = convert(char(20), m.creation_date, 100), ' + \
      'mdate = convert(char(20), m.modification_date, 100) ' + \
      'from #strains s, PRB_Strain_Synonym m ' + \
      'where s._Strain_key = m._Strain_key')

cmds.append('select _StrainType_key, strainType, ' + \
      'cdate = convert(char(20), creation_date, 100), ' + \
      'mdate = convert(char(20), modification_date, 100) ' + \
      'from MLP_StrainType')

cmds.append('select m._Strain_key, m._StrainType_key, s.private, ' + \
      'cdate = convert(char(20), m.creation_date, 100), ' + \
      'mdate = convert(char(20), m.modification_date, 100) ' + \
      'from #strains s, MLP_StrainTypes m ' + \
      'where s._Strain_key = m._Strain_key')

cmds.append('select distinct a.accID, a.LogicalDB, a._Object_key, a.preferred, s.private, ' + \
      'cdate = convert(varchar(20), a.creation_date, 100), ' + \
      'mdate = convert(varchar(20), a.modification_date, 100) ' + \
      'from #strains s, ACC_Accession a ' + \
      'where s._Strain_key = a._Object_key ' + \
      'and a._MGIType_key = 10')

cmds.append('select _Species_key, species, ' + \
      'cdate = convert(char(20), creation_date, 100), ' + \
      'mdate = convert(char(20), modification_date, 100) ' + \
      'from MLP_Species')

results = db.sql(cmds, 'auto', execute = 1)

for r in results[2]:
	fp1.write(`r['_Strain_key']` + TAB + \
	         `r['_Species_key']` + TAB + \
	         strip_newline(r['strain']) + TAB + \
	         `r['private']` + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[3]:
	fp2.write(`r['_Strain_key']` + TAB + \
	         `r['_Marker_key']` + TAB + \
	         mgi_utils.prvalue(r['_Allele_key']) + TAB + \
	         `r['private']` + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[4]:
	fp3.write(`r['_Synonym_key']` + TAB + \
	         `r['_Strain_key']` + TAB + \
	         `r['private']` + TAB + \
	         strip_newline(r['synonym']) + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[5]:
	fp4.write(`r['_StrainType_key']` + TAB + \
	         r['strainType'] + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[6]:
	fp5.write(`r['_Strain_key']` + TAB + \
	         `r['_StrainType_key']` + TAB + \
	         `r['private']` + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[7]:
	fp6.write(r['accID'] + TAB + \
		 r['LogicalDB'] + TAB + \
	         `r['_Object_key']` + TAB + \
	         `r['preferred']` + TAB + \
	         `r['private']` + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

for r in results[8]:
	fp7.write(`r['_Species_key']` + TAB + \
	         r['species'] + TAB + \
		 r['cdate'] + TAB + \
		 r['mdate'] + CRT)

fp1.close
fp2.close
fp3.close
fp4.close
fp5.close
fp6.close
fp7.close

