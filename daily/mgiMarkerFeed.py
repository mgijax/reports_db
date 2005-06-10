#!/usr/local/bin/python

'''
#
# Report:
#       Create several tab-delimited (bcp) files for BioDataMart feed
#	TR 2794
#	TR 3137 - add MRK_Offset
#	TR 3345 - add Strain info; incorporate marker offset into marker.bcp
#	TR 5565 - JaxStrain additions
#	MPR
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
#	9. allele_pairstate.bcp
#	10.allele.bcp
#	11. allele_label.bcp
#	12. allele_pair.bcp
#	13. allele_note.bcp
#	14. accession_allele.bcp
#	15. strain_marker.bcp
#	16. strain_synonym.bcp
#	17. strain.bcp
#	18. strain_type.bcp
#	19. strain_strain_type.bcp
#	20. strain_genotype.bcp
#	21. accession_strain.bcp
#	22. strain_species.bcp
#	23. reference.bcp
#	24. genotype_mpt_reference.bcp
#	25. allele_reference.bcp
#	26. marker_reference.bcp
#	27. strain_reference.bcp
#
# Usage:
#       mgiMarkerFeed.py
#
# Notes:
#
# History:
#
# lec	06/01/2005
#	- added all markers (not just mouse) per cjb
#
# lec	01/14/2004
#	- TR 5565; JaxStrain additions
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

COLDELIM = "&=&"
LINEDELIM = "\n#=#"

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

def vocabs():

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
    fp.close()

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
    fp.close()

    #
    # allele_type
    #

    fp = open(OUTPUTDIR + 'allele_type.bcp', 'w')

    results = db.sql('select _Term_key, term, ' + \
	    'cdate = convert(char(20), creation_date, 100), ' + \
	    'mdate = convert(char(20), modification_date, 100) ' + \
	    'from VOC_Term_ALLType_View', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # allele_cellline
    #

    fp = open(OUTPUTDIR + 'allele_cellline.bcp', 'w')
    
    results = db.sql('select _CellLine_key, cellLine, _Strain_key, provider, isMutant, ' + \
	    'cdate = convert(char(20), creation_date, 100), ' + \
	    'mdate = convert(char(20), modification_date, 100) ' + \
	    'from ALL_CellLine', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_CellLine_key']` + TAB + \
		     r['cellLine'] + TAB + \
		     `r['_Strain_key']` + TAB + \
		     mgi_utils.prvalue(r['provider']) + TAB + \
		     `r['isMutant']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # allele_inheritance_mode
    #
    
    fp = open(OUTPUTDIR + 'allele_inheritance_mode.bcp', 'w')

    results = db.sql('select _Term_key, term, ' + \
	    'cdate = convert(char(20), creation_date, 100), ' + \
	    'mdate = convert(char(20), modification_date, 100) ' + \
	    'from VOC_Term_ALLInheritMode_View', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # allele_pairstate
    #
    
    fp = open(OUTPUTDIR + 'allele_pairstate.bcp', 'w')

    results = db.sql('select _Term_key, term, ' + \
	    'cdate = convert(char(20), creation_date, 100), ' + \
	    'mdate = convert(char(20), modification_date, 100) ' + \
	    'from VOC_Term_ALLPairState_View', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # strain_type.bcp
    #
    
    fp = open(OUTPUTDIR + 'strain_type.bcp', 'w')
    
    results = db.sql('select _Term_key, term, ' + \
          'cdate = convert(char(20), creation_date, 100), ' + \
          'mdate = convert(char(20), modification_date, 100) ' + \
          'from VOC_Term_StrainType_View', 'auto')

    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
	             r['term'] + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # mp_term.bcp
    #
    
    fp = open(OUTPUTDIR + 'mp_term.bcp', 'w')
    
    results = db.sql('select x._Term_key, x.sequenceNum, note = rtrim(x.note) ' + \
          'from VOC_Term t, VOC_Text x ' + \
	  'where t._Vocab_key = 5 ' + \
	  'and t._Term_key = x._Term_key ' + \
	  'order by x._Term_key, x.sequenceNum', 'auto')

    notes = {}
    for r in results:
	key = r['_Term_key']
	value = r['note']
	if not notes.has_key(key):
	    notes[key] = []
	notes[key].append(value)

    results = db.sql('select a.accID, t.term, t._Term_key, ' +
          'cdate = convert(char(20), t.creation_date, 100), ' + \
          'mdate = convert(char(20), t.modification_date, 100) ' + \
          'from VOC_Term t, ACC_Accession a ' + \
	  'where t._Vocab_key = 5 ' + \
	  'and t._Term_key = a._Object_key ' + \
	  'and a._MGIType_key = 13 ' + \
	  'and a.preferred = 1', 'auto')

    for r in results:
	    key = r['_Term_key']

	    fp.write(`key` + TAB + \
		     r['accID'] + TAB + \
	             r['term'] + TAB)

            if notes.has_key(key):
		fp.write(string.join(notes[key], ''))
            fp.write(TAB)

	    fp.write(r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

def markers():

    #
    # marker
    # marker_label
    # accession_marker
    #

    #
    # select all mouse markers which have a preferred MGI Accession ID.
    # this will include any splits (since the MGI Acc ID stays with 
    # the split symbol).
    #

    db.sql('select m._Marker_key, m._Organism_key, m.symbol, m.name ' + \
	    'into #markers ' + \
	    'from MRK_Marker m ' + \
	    'where m._Organism_key = 1 ' + \
	    'and exists (select 1 from ACC_Accession a ' + \
	    'where m._Marker_key = a._Object_key ' + \
	    'and a._MGIType_key = 2 ' + \
	    'and a.prefixPart = "MGI:" ' + \
	    'and a._LogicalDB_key = 1 ' + \
	    'and a.preferred = 1) ' + \
	    'union ' + \
	    'select m._Marker_key, m._Organism_key, m.symbol, m.name ' + \
	    'from MRK_Marker m ' + \
	    'where m._Organism_key != 1', None)

    db.sql('create index idx1 on #markers(_Marker_key)', None)

    #
    # select data fields for marker.bcp
    #

    fp = open(OUTPUTDIR + 'marker.bcp', 'w')

    results = db.sql('select k._Marker_key, k._Organism_key, m._Marker_Type_key, s.status, ' + \
	    'm.symbol, m.name, m.chromosome, m.cytogeneticOffset, o.offset, ' + \
	    'cdate = convert(char(20), m.creation_date, 100), ' + \
	    'mdate = convert(char(20), m.modification_date, 100) ' + \
	    'from #markers k, MRK_Marker m, MRK_Status s, MRK_Offset o ' + \
	    'where k._Marker_key = m._Marker_key ' + \
	    'and m._Marker_Status_key = s._Marker_Status_key ' + \
	    'and k._Marker_key = o._Marker_key ' + \
	    'and o.source = 0 ' + \
	    'union ' + \
            'select m._Marker_key, m._Organism_key, m._Marker_Type_key, s.status, ' + \
	    'm.symbol, m.name, m.chromosome, m.cytogeneticOffset,  null, ' + \
	    'cdate = convert(char(20), m.creation_date, 100), ' + \
	    'mdate = convert(char(20), m.modification_date, 100) ' + \
	    'from #markers k, MRK_Marker m, MRK_Status s ' + \
	    'where k._Organism_key != 1 ' + \
	    'and k._Marker_key = m._Marker_key ' + \
	    'and m._Marker_Status_key = s._Marker_Status_key', 'auto')

    for r in results:
	    fp.write(`r['_Marker_key']` + TAB + \
	             `r['_Organism_key']` + TAB + \
	             `r['_Marker_Type_key']` + TAB + \
		     mgi_status(r['status']) + TAB + \
		     strip_newline(r['symbol']) + TAB + \
		     strip_newline(r['name']) + TAB + \
		     r['chromosome'] + TAB + \
		     mgi_utils.prvalue(r['cytogeneticOffset']) + TAB + \
		     mgi_utils.prvalue(r['offset']) + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)

    fp.close()

    #
    # select data fields for marker_label.bcp
    # status = 1 if the label matches the MRK_Marker.symbol or MRK_Marker.name
    # field, else status = 0
    #

    fp = open(OUTPUTDIR + 'marker_label.bcp', 'w')

    results = db.sql('select distinct m._Marker_key, m.label, m.labelType, status = 1, ' + \
	    'cdate = convert(char(20), m.creation_date, 100), ' + \
	    'mdate = convert(char(20), m.modification_date, 100) ' + \
	    'from #markers k, MRK_Label m ' + \
	    'where k._Marker_key = m._Marker_key ' + \
	    'and m.labelType = "MS" ' + \
	    'and k.symbol = m.label ', 'auto')

    for r in results:
	    fp.write(`r['_Marker_key']` + TAB + \
		     strip_newline(r['label']) + TAB + \
		     r['labelType'] + TAB + \
		     `r['status']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)

    results = db.sql('select distinct m._Marker_key, m.label, m.labelType, status = 1, ' + \
	    'cdate = convert(char(20), m.creation_date, 100), ' + \
	    'mdate = convert(char(20), m.modification_date, 100) ' + \
	    'from #markers k, MRK_Label m ' + \
	    'where k._Marker_key = m._Marker_key ' + \
	    'and m.labelType = "MN" ' + \
	    'and k.name = m.label', 'auto')

    for r in results:
	    fp.write(`r['_Marker_key']` + TAB + \
		     strip_newline(r['label']) + TAB + \
		     r['labelType'] + TAB + \
		     `r['status']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)

    results = db.sql('select distinct m._Marker_key, m.label, m.labelType, status = 0, ' + \
	    'cdate = convert(char(20), m.creation_date, 100), ' + \
	    'mdate = convert(char(20), m.modification_date, 100) ' + \
	    'from #markers k, MRK_Label m ' + \
	    'where k._Marker_key = m._Marker_key ' + \
	    'and m.labelType in ("MS", "MN", "MY") ' + \
	    'and k.symbol != m.label and k.name != m.label ', 'auto')

    for r in results:
	    fp.write(`r['_Marker_key']` + TAB + \
		     strip_newline(r['label']) + TAB + \
		     r['labelType'] + TAB + \
		     `r['status']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)

    fp.close()

    #
    # select data fields for accession_marker.bcp
    # select MGI Accession IDs only
    #

    fp = open(OUTPUTDIR + 'accession_marker.bcp', 'w')

    results = db.sql('select m.accID, LogicalDB = l.name, m._Object_key, m.preferred, ' + \
	    'cdate = convert(char(20), m.creation_date, 100), ' + \
	    'mdate = convert(char(20), m.modification_date, 100) ' + \
	    'from #markers k, ACC_Accession m, ACC_LogicalDB l ' + \
	    'where k._Marker_key = m._Object_key ' + \
	    'and m._MGIType_key = 2 ' + \
	    'and m.prefixPart = "MGI:" ' + \
	    'and m._LogicalDB_key = 1 ' + \
	    'and m._LogicalDB_key = l._LogicalDB_key ' + \
	    'union ' + \
            'select m.accID, LogicalDB = l.name, m._Object_key, m.preferred, ' + \
	    'cdate = convert(char(20), m.creation_date, 100), ' + \
	    'mdate = convert(char(20), m.modification_date, 100) ' + \
	    'from #markers k, ACC_Accession m, ACC_LogicalDB l ' + \
	    'where k._Organism_key != 1 ' + \
	    'and k._Marker_key = m._Object_key ' + \
	    'and m._MGIType_key = 2 ' + \
	    'and m._LogicalDB_key in (55, 64, 15) ' + \
	    'and m._LogicalDB_key = l._LogicalDB_key ', 'auto')

    for r in results:
	    fp.write(r['accID'] + TAB + \
		     r['LogicalDB'] + TAB + \
		     `r['_Object_key']` + TAB + \
		     `r['preferred']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

def alleles():

    #
    # allele
    # allele_label
    # allele_pair
    # allele_note
    # accession_allele
    #

    #
    # select all alleles with a status of "approved"
    # all other statuses are private/confidential alleles
    #

    db.sql('select m._Allele_key into #alleles ' + \
	'from ALL_Allele m, VOC_Term t ' + \
	'where m._Marker_key is not null and m._Allele_Status_key = t._Term_key and t.term = "Approved" ', None)
    db.sql('create index idx1 on #alleles(_Allele_key)', None)

    #
    # select data fields for allele.bcp
    #

    fp = open(OUTPUTDIR + 'allele.bcp', 'w')

    results = db.sql('select m._Allele_key, m._Marker_key, m._Mode_key, m._Allele_Type_key, ' + \
	'm._ESCellLine_key, m._MutantESCellLine_key, ' + \
	'm._Strain_key, m.isWildType, m.symbol, m.name, status = t.term, ' + \
	'cdate = convert(char(20), m.creation_date, 100), ' + \
	'mdate = convert(char(20), m.modification_date, 100) ' + \
	'from #alleles a, ALL_Allele m, VOC_Term t ' + \
	'where a._Allele_key = m._Allele_key ' + \
	'and m._Allele_Status_key = t._Term_key', 'auto')

    for r in results:
	    fp.write(`r['_Allele_key']` + TAB + \
	             `r['_Marker_key']` + TAB + \
	             `r['_Mode_key']` + TAB + \
	             `r['_Allele_Type_key']` + TAB + \
	             `r['_ESCellLine_key']` + TAB + \
	             `r['_MutantESCellLine_key']` + TAB + \
		     mgi_status(r['status']) + TAB + \
		     `r['_Strain_key']` + TAB + \
		     `r['isWildType']` + TAB + \
		     strip_newline(r['symbol']) + TAB + \
		     strip_newline(r['name']) + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # select data fields for allele_label.bcp
    # this includes all allele symbols, names and synonyms
    # status = 1 if the label matches the ALL_Allele.symbol or ALL_Allele.name
    # field, else status = 0
    #

    fp = open(OUTPUTDIR + 'allele_label.bcp', 'w')

    results = db.sql('select m._Allele_key, m.name, labelType = "AN", status = 1, ' + \
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
	    'from #alleles a, ALL_Allele m, MGI_Synonym s ' + \
	    'where a._Allele_key = m._Allele_key ' + \
	    'and a._Allele_key = s._Object_key ' + \
	    'and s._MGIType_key = 11', 'auto')

    for r in results:
	    fp.write(`r['_Allele_key']` + TAB + \
		     strip_newline(r['name']) + TAB + \
		     r['labelType'] + TAB + \
		     `r['status']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # select data fields for allele_pair.bcp
    #

    fp = open(OUTPUTDIR + 'allele_pair.bcp', 'w')

    results = db.sql('select p.*, ' + 
	    'cdate = convert(char(20), p.creation_date, 100), ' + \
	    'mdate = convert(char(20), p.modification_date, 100) ' + \
	    'from #alleles a, GXD_AllelePair p ' +  \
	    'where a._Allele_key = p._Allele_key_1 ' + \
	    'union ' + \
            'select p.*, ' + 
	    'cdate = convert(char(20), p.creation_date, 100), ' + \
	    'mdate = convert(char(20), p.modification_date, 100) ' + \
	    'from #alleles a, GXD_AllelePair p ' +  \
	    'where a._Allele_key = p._Allele_key_2 ' + \
	    'order by _Genotype_key, sequenceNum', 'auto')

    for r in results:
	    fp.write(`r['_AllelePair_key']` + TAB + \
		     `r['_Allele_key_1']` + TAB + \
		     mgi_utils.prvalue(r['_Allele_key_2']) + TAB + \
		     `r['_Marker_key']` + TAB + \
		     `r['_Genotype_key']` + TAB + \
		     `r['_PairState_key']` + TAB + \
		     `r['sequenceNum']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # select data fields for allele_note.bcp
    #

    fp = open(OUTPUTDIR + 'allele_note.bcp', 'w')

    results = db.sql('select a._Allele_key, note = rtrim(nc.note), nc.sequenceNum, nt.noteType, ' + \
	    'cdate = convert(char(20), n.creation_date, 100), ' + \
	    'mdate = convert(char(20), n.modification_date, 100) ' + \
	    'from #alleles a, MGI_Note n, MGI_NoteChunk nc, MGI_NoteType nt ' + \
	    'where a._Allele_key = n._Object_key ' + \
	    'and n._MGIType_key = 11 ' + \
	    'and n._Note_key = nc._Note_key ' + \
	    'and n._NoteType_key = nt._NoteType_key', 'auto')

    for r in results:
	    fp.write(`r['_Allele_key']` + COLDELIM + \
		     r['noteType'] + COLDELIM + \
		     `r['sequenceNum']` + COLDELIM + \
		     r['note'] + COLDELIM + \
		     r['cdate'] + COLDELIM + \
		     r['mdate'] + LINEDELIM)
    fp.close()

    #
    # select data fields for accession_allele.bcp
    #

    fp = open(OUTPUTDIR + 'accession_allele.bcp', 'w')

    results = db.sql('select m.accID, LogicalDB = l.name, m._Object_key, m.preferred, ' + \
	    'cdate = convert(char(20), m.creation_date, 100), ' + \
	    'mdate = convert(char(20), m.modification_date, 100) ' + \
	    'from #alleles k, ACC_Accession m, ACC_LogicalDB l ' + \
	    'where k._Allele_key = m._Object_key ' + \
	    'and m._MGIType_key = 11 ' + \
	    'and m.prefixPart = "MGI:" ' + \
	    'and m._LogicalDB_key = 1 ' + \
	    'and m._LogicalDB_key = l._LogicalDB_key', 'auto')

    for r in results:
	    fp.write(r['accID'] + TAB + \
		     r['LogicalDB'] + TAB + \
		     `r['_Object_key']` + TAB + \
		     `r['preferred']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

def strains():

    #
    # strain.bcp
    # strain_marker.bcp
    # strain_synonym.bcp
    # strain_strain_type.bcp
    # accession_strain.bcp
    # strain_species.bcp
    # strain_genotype.bcp
    #

    #
    # select all strains which have a Jax Registry ID or MMRRC ID
    # plus all strains which are cross-referenced by Allele or Allele CellLine
    # strain.bcp
    #

    db.sql('select distinct s._Strain_key, s._Species_key, s.strain, s.private, ' + \
          'cdate = convert(char(20), s.creation_date, 100), ' + \
          'mdate = convert(char(20), s.modification_date, 100) ' + \
          'into #strains ' + \
          'from PRB_Strain s, ACC_Accession a ' + \
          'where s._Strain_key = a._Object_key ' + \
          'and a._MGIType_key = 10 ' + \
          'and a._LogicalDB_key in (22, 38) ' + \
          'union ' + \
          'select distinct s._Strain_key, s._Species_key, s.strain, s.private, ' + \
          'cdate = convert(char(20), s.creation_date, 100), ' + \
          'mdate = convert(char(20), s.modification_date, 100) ' + \
          'from PRB_Strain s, ALL_Allele a ' + \
          'where s._Strain_key = a._Strain_key ' + \
          'union ' + \
          'select distinct s._Strain_key, s._Species_key, s.strain, s.private, ' + \
          'cdate = convert(char(20), s.creation_date, 100), ' + \
          'mdate = convert(char(20), s.modification_date, 100) ' + \
          'from PRB_Strain s, ALL_CellLine a ' + \
          'where s._Strain_key = a._Strain_key ', None)

    db.sql('create index idx1 on #strains(_Strain_key)', None)

    #
    # strain.bcp
    #

    fp = open(OUTPUTDIR + 'strain.bcp', 'w')

    results = db.sql('select * from #strains', 'auto')

    for r in results:
	    fp.write(`r['_Strain_key']` + TAB + \
	             `r['_Species_key']` + TAB + \
	             strip_newline(r['strain']) + TAB + \
	             `r['private']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # strain_marker.bcp
    #

    fp = open(OUTPUTDIR + 'strain_marker.bcp', 'w')

    results = db.sql('select distinct m._Strain_key, m._Marker_key, m._Allele_key, s.private, ' + \
          'cdate = convert(char(20), m.creation_date, 100), ' + \
          'mdate = convert(char(20), m.modification_date, 100) ' + \
          'from #strains s, PRB_Strain_Marker m ' + \
          'where s._Strain_key = m._Strain_key', 'auto')

    for r in results:
	    fp.write(`r['_Strain_key']` + TAB + \
	             `r['_Marker_key']` + TAB + \
	             mgi_utils.prvalue(r['_Allele_key']) + TAB + \
	             `r['private']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # strain_synonym.bcp
    #

    fp = open(OUTPUTDIR + 'strain_synonym.bcp', 'w')
    
    results = db.sql('select m._Synonym_key, m._Object_key, m.synonym, s.private, m.synonymType, ' + \
          'cdate = convert(char(20), m.creation_date, 100), ' + \
          'mdate = convert(char(20), m.modification_date, 100) ' + \
          'from #strains s, MGI_Synonym_Strain_View m ' + \
          'where s._Strain_key = m._Object_key', 'auto')

    for r in results:
	    fp.write(`r['_Synonym_key']` + TAB + \
	             `r['_Object_key']` + TAB + \
	             `r['private']` + TAB + \
	             strip_newline(r['synonym']) + TAB + \
	             strip_newline(r['synonymType']) + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # strain_strain_type.bcp
    #

    fp = open(OUTPUTDIR + 'strain_strain_type.bcp', 'w')

    results = db.sql('select m._Strain_key, m._StrainType_key, s.private, ' + \
          'cdate = convert(char(20), m.creation_date, 100), ' + \
          'mdate = convert(char(20), m.modification_date, 100) ' + \
          'from #strains s, PRB_Strain_Type m ' + \
          'where s._Strain_key = m._Strain_key', 'auto')

    for r in results:
	    fp.write(`r['_Strain_key']` + TAB + \
	             `r['_StrainType_key']` + TAB + \
	             `r['private']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # accession_strain.bcp
    #

    fp = open(OUTPUTDIR + 'accession_strain.bcp', 'w')
    
    results = db.sql('select distinct a.accID, LogicalDB = l.name, a._Object_key, a.preferred, s.private, ' + \
          'cdate = convert(varchar(20), a.creation_date, 100), ' + \
          'mdate = convert(varchar(20), a.modification_date, 100) ' + \
          'from #strains s, ACC_Accession a, ACC_LogicalDB l ' + \
          'where s._Strain_key = a._Object_key ' + \
          'and a._MGIType_key = 10 ' + \
          'and a._LogicalDB_key = l._LogicalDB_key', 'auto')

    for r in results:
	    fp.write(r['accID'] + TAB + \
		     r['LogicalDB'] + TAB + \
	             `r['_Object_key']` + TAB + \
	             `r['preferred']` + TAB + \
	             `r['private']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()
    
    #
    # strain_species.bcp
    #
    
    fp = open(OUTPUTDIR + 'strain_species.bcp', 'w')

    results = db.sql('select _Term_key, term, ' + \
          'cdate = convert(char(20), creation_date, 100), ' + \
          'mdate = convert(char(20), modification_date, 100) ' + \
          'from VOC_Term_StrainSpecies_View', 'auto')
    
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
	             r['term'] + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

def genotypes():

    db.sql('select distinct g._Genotype_key ' + \
	    'into #genotypes ' + \
	    'from #strains s, GXD_Genotype g, VOC_Annot a ' + \
	    'where s._Strain_key = g._Strain_key ' + \
	    'and g._Genotype_key = a._Object_key ' + \
	    'and a._AnnotType_key = 1002 ' + \
	    'union ' + \
	    'select distinct g._Genotype_key ' + \
	    'from #strains s, PRB_Strain_Genotype g ' + \
	    'where s._Strain_key = g._Strain_key', 'auto')

    db.sql('create index idex1 on #genotypes(_Genotype_key)', None)

    #
    # genotype.bcp
    #

    fp = open(OUTPUTDIR + 'genotype.bcp', 'w')

    results = db.sql('select g._Genotype_key, s.strain, p.isConditional, ' + \
          'cdate = convert(char(20), p.creation_date, 100), ' + \
          'mdate = convert(char(20), p.modification_date, 100) ' + \
	  'from #genotypes g, GXD_Genotype p, PRB_Strain s ' + \
	  'where g._Genotype_key = p._Genotype_key ' + \
	  'and p._Strain_key = s._Strain_key', 'auto')

    for r in results:
	fp.write(`r['_Genotype_key']` + TAB + \
		r['strain'] + TAB + \
		`r['isConditional']` + TAB + \
		r['cdate'] + TAB + \
		r['mdate'] + CRT)
    fp.close()

    #
    # genotype_mpt.bcp
    #

    fp = open(OUTPUTDIR + 'genotype_mpt.bcp', 'w')

    results = db.sql('select g._Genotype_key, a._Annot_key, a._Term_key, a.isNot, ' + \
          'cdate = convert(char(20), a.creation_date, 100), ' + \
          'mdate = convert(char(20), a.modification_date, 100) ' + \
	'from #genotypes g, VOC_Annot a ' + \
	'where g._Genotype_key = a._Object_key ' + \
	'and a._AnnotType_key = 1002 ', 'auto')

    for r in results:
	fp.write(`r['_Annot_key']` + TAB + \
	         `r['_Term_key']` + TAB + \
	         `r['_Genotype_key']` + TAB + \
		`r['isNot']` + TAB + \
		r['cdate'] + TAB + \
		r['mdate'] + CRT)

    fp.close()

    #
    # strain_genotype.bcp
    #
    
    fp = open(OUTPUTDIR + 'strain_genotype.bcp', 'w')

    results = db.sql('select p._Strain_key, p._Genotype_key, t.term, ' + \
          'cdate = convert(char(20), p.creation_date, 100), ' + \
          'mdate = convert(char(20), p.modification_date, 100) ' + \
          'from #strains s, #genotypes g, PRB_Strain_Genotype p, VOC_Term t ' + \
          'where s._Strain_key = p._Strain_key ' + \
	  'and g._Genotype_key = p._Genotype_key ' + \
          'and p._Qualifier_key = t._Term_key', 'auto')
    
    for r in results:
            fp.write(`r['_Strain_key']` + TAB + \
                     `r['_Genotype_key']` + TAB + \
                     r['term'] + TAB + \
                     r['cdate'] + TAB + \
                     r['mdate'] + CRT)
    fp.close()

def references():

    #
    # reference.bcp
    #

    fp = open(OUTPUTDIR + 'reference.bcp', 'w')

    #
    # references annotated to a Genotype
    #

    db.sql('select distinct e._Refs_key, g._Genotype_key, a._Annot_key, e._AnnotEvidence_key, ' + \
	    'cdate = convert(char(20), e.creation_date, 100), ' + \
	    'mdate = convert(char(20), e.modification_date, 100) ' + \
	    'into #genoreferences ' + \
	    'from #genotypes g, VOC_Annot a, VOC_Evidence e ' + \
	    'where g._Genotype_key = a._Object_key ' + \
	    'and a._AnnotType_key = 1002 ' + \
	    'and a._Annot_key = e._Annot_key', None)

    #
    # references annotated to an Allele via the Genotype
    #

    db.sql('select r._Refs_key, ag._Allele_key, rt.assocType, ' + \
	    'cdate = convert(char(20), r.creation_date, 100), ' + \
	    'mdate = convert(char(20), r.modification_date, 100) ' + \
	    'into #allreferences ' + \
	    'from #genotypes g, GXD_AlleleGenotype ag, MGI_Reference_Assoc r, MGI_RefAssocType rt ' + \
	    'where g._Genotype_key = ag._Genotype_key ' + \
	    'and ag._Allele_key = r._Object_key ' + \
	    'and r._MGIType_key = 11 ' + \
	    'and r._RefAssocType_key = rt._RefAssocType_key', None)

    #
    # references annotated to Strains
    #

    db.sql('select r._Refs_key, r._Object_key, rt.assocType, ' + \
            'cdate = convert(char(20), r.creation_date, 100), ' + \
            'mdate = convert(char(20), r.modification_date, 100) ' + \
            'into #strainreferences ' + \
            'from #strains s, MGI_Reference_Assoc r, MGI_RefAssocType rt ' + \
            'where s._Strain_key = r._Object_key ' + \
            'and r._MGIType_key = 10 ' + \
            'and r._RefAssocType_key = rt._RefAssocType_key', None)

    #
    # references annotated to a Marker via the Genotype
    #

    db.sql('select r._Refs_key, r._Marker_key, ' + \
	    'cdate = convert(char(20), r.creation_date, 100), ' + \
	    'mdate = convert(char(20), r.modification_date, 100) ' + \
	    'into #mrkreferences ' + \
	    'from #genotypes g, GXD_AlleleGenotype ag, MRK_Reference r ' + \
	    'where g._Genotype_key = ag._Genotype_key ' + \
	    'and ag._Marker_key = r._Marker_key', None)

    db.sql('create index idx1 on #genoreferences(_Refs_key)', None)
    db.sql('create index idx1 on #allreferences(_Refs_key)', None)
    db.sql('create index idx1 on #strainreferences(_Refs_key)', None)
    db.sql('create index idx1 on #mrkreferences(_Refs_key)', None)

    db.sql('select distinct _Refs_key into #references from #genoreferences ' + \
	    'union select distinct _Refs_key from #allreferences ' + \
	    'union select distinct _Refs_key from #strainreferences ' + \
	    'union select distinct _Refs_key from #mrkreferences', None)
    db.sql('create index idx1 on #references(_Refs_key)', None)

    results = db.sql('select r._Refs_key, b.refType, b.authors, b.authors2, ' + \
	    'b.title, b.title2, b.journal, b.vol, b.issue, b.pgs, b.year, ' + \
	    'b.isReviewArticle, ' + \
	    'k.book_au, k.book_title, k.publisher, k.place, k.series_ed, ' + \
	    'cdate = convert(char(20), b.creation_date, 100), ' + \
	    'mdate = convert(char(20), b.modification_date, 100) ' + \
	    'from #references r, BIB_Refs b, BIB_Books k ' + \
	    'where r._Refs_key = b._Refs_key ' + \
	    'and b._Refs_key *= k._Refs_key', 'auto')

    for r in results:
	    fp.write(`r['_Refs_key']` + TAB + \
	             r['refType'] + TAB + \
		     mgi_utils.prvalue(r['authors']) + mgi_utils.prvalue(r['authors2']) + TAB + \
		     mgi_utils.prvalue(r['title']) + mgi_utils.prvalue(r['title2']) + TAB + \
		     mgi_utils.prvalue(r['journal']) + TAB + \
		     mgi_utils.prvalue(r['vol']) + TAB + \
		     mgi_utils.prvalue(r['issue']) + TAB + \
		     mgi_utils.prvalue(r['pgs']) + TAB + \
	             `r['year']` + TAB + \
	             `r['isReviewArticle']` + TAB + \
		     mgi_utils.prvalue(r['book_au']) + TAB + \
		     mgi_utils.prvalue(r['book_title']) + TAB + \
		     mgi_utils.prvalue(r['publisher']) + TAB + \
		     mgi_utils.prvalue(r['place']) + TAB + \
		     mgi_utils.prvalue(r['series_ed']) + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # accession_reference.bcp
    #

    fp = open(OUTPUTDIR + 'accession_reference.bcp', 'w')
    results = db.sql('select a.accID, LogicalDB = l.name, a._Object_key, a.preferred, ' + \
	    'cdate = convert(char(20), a.creation_date, 100), ' + \
	    'mdate = convert(char(20), a.modification_date, 100) ' + \
	    'from #references r, ACC_Accession a, ACC_LogicalDB l ' + \
	    'where r._Refs_key = a._Object_key ' + \
	    'and a._MGIType_key = 1 ' + \
	    'and a._LogicalDB_key = l._LogicalDB_key', 'auto')

    for r in results:
	    fp.write(r['accID'] + TAB + \
		     r['LogicalDB'] + TAB + \
		     `r['_Object_key']` + TAB + \
		     `r['preferred']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # genotype_mpt_reference.bcp
    #

    fp = open(OUTPUTDIR + 'genotype_mpt_reference.bcp', 'w')

    results = db.sql('select g._Refs_key, g._Genotype_key, g._Annot_key, g.cdate, g.mdate ' + \
	'from #genoreferences g ' + \
	'where not exists (select 1 from MGI_Note_VocEvidence_View n ' + \
	'where g._AnnotEvidence_key = n._Object_key) ' + \
	'order by g._Genotype_key ', 'auto')
    for r in results:
	    fp.write(`r['_Annot_key']` + COLDELIM + \
		     `r['_Refs_key']` + COLDELIM + \
		     COLDELIM + \
		     COLDELIM + \
		     COLDELIM + \
		     COLDELIM + \
		     r['cdate'] + COLDELIM + \
		     r['mdate'] + LINEDELIM)

    results = db.sql('select g._Refs_key, g._Genotype_key, g._Annot_key, g._AnnotEvidence_key, g.cdate, g.mdate, ' + \
	'n._Note_key, n.noteType, n.sequenceNum, n.note ' + \
	'from #genoreferences g, MGI_Note_VocEvidence_View n ' + \
	'where g._AnnotEvidence_key = n._Object_key ' + \
	'order by g._Genotype_key, g._AnnotEvidence_key, n.sequenceNum', 'auto')
    for r in results:
	    fp.write(`r['_Annot_key']` + COLDELIM + \
		     `r['_Refs_key']` + COLDELIM + \
		     `r['_Note_key']` + COLDELIM + \
		     r['noteType'] + COLDELIM + \
		     `r['sequenceNum']` + COLDELIM + \
		     r['note'] + COLDELIM + \
		     r['cdate'] + COLDELIM + \
		     r['mdate'] + LINEDELIM)
    fp.close()

    #
    # allele_reference.bcp
    #

    fp = open(OUTPUTDIR + 'allele_reference.bcp', 'w')

    results = db.sql('select * from #allreferences', 'auto')
    for r in results:
	    fp.write(`r['_Allele_key']` + TAB + \
		     `r['_Refs_key']` + TAB + \
		     r['assocType'] + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

    #
    # strain_reference.bcp
    #

    fp = open(OUTPUTDIR + 'strain_reference.bcp', 'w')

    results = db.sql('select * from #strainreferences', 'auto')
    for r in results:
            fp.write(`r['_Object_key']` + TAB + \
                     `r['_Refs_key']` + TAB + \
                     r['assocType'] + TAB + \
                     r['cdate'] + TAB + \
                     r['mdate'] + CRT)
    fp.close()

    #
    # marker_reference.bcp
    #

    fp = open(OUTPUTDIR + 'marker_reference.bcp', 'w')

    results = db.sql('select * from #mrkreferences', 'auto')
    for r in results:
	    fp.write(`r['_Marker_key']` + TAB + \
		     `r['_Refs_key']` + TAB + \
		     r['cdate'] + TAB + \
		     r['mdate'] + CRT)
    fp.close()

#
# Main
#

db.useOneConnection(1)
db.set_sqlLogFunction(db.sqlLogAll)
vocabs()
alleles()
markers()
strains()
genotypes()
references()
db.useOneConnection(0)

