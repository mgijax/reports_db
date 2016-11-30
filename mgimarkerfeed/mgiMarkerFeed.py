#!/usr/local/bin/python

'''
#
# Report:
#       Create several tab-delimited (bcp) files for BioDataMart feed
#	TR12427/Disase Ontology (DO)
#	TR 2794
#	TR 3137 - add MRK_Offset
#	TR 3345 - add Strain info; incorporate marker offset into marker.bcp
#	TR 5565 - JaxStrain additions
#	TR 10460 - add mutant cell line accession ids
#	TR 11515 - added Allele/Collection
#	TR 12027 - added MGI_Relationship
#	MPR
#
# bcp files
#	1.  marker_type.bcp
#	2.  species.bcp
#	3.  allele_type.bcp
#	4.  allele_inheritance_mode.bcp
#	5.  allele_pairstate.bcp
#       6.  allele_transmission.bcp
#       7.  allele_collection.bcp
#       8.  allele_subtype.bcp
#	9.  allele_creator.bcp
#	10. allele_celllinetype.bcp
#	11. allele_vector.bcp
#	12. allele_vectortype.bcp
#       13. genotype_existsas.bcp
#	14. strain_type.bcp
#	15. mp_term.bcp
#	16. mp_synonym.bcp
#	17. mp_closure.bcp
#	18. marker.bcp
#	19. marker_label.bcp
#	20. accession_marker.bcp
#	21. allele_cellline.bcp
#	22. accession_allele_cellline.bcp
#	23. allele_derivation.bcp
#	24. allele.bcp
#	25. allele_allele_cellline.bcp
#	26. allele_label.bcp
#	27. allele_note.bcp
#	28. accession_allele.bcp
#	29. allele_subtype_assoc.bcp
#	30. strain.bcp
#	31. strain_marker.bcp
#	32. strain_synonym.bcp
#	33. strain_strain_type.bcp
#	34. accession_strain.bcp
#	35. strain_species.bcp
#       36. genotype.bcp
#       37. genotype_mpt.bcp
#       38. genotype_header.bcp
#       39. strain_genotype.bcp
#       40. allele_pair.bcp
#	41. reference.bcp
#	42. accession_reference.bcp
#	43. genotype_mpt_reference.bcp
#	44. allele_reference.bcp
#	45. strain_reference.bcp
#	46. marker_reference.bcp
#	47. marker_omim.bcp
#	48: mgi_relationship/mgi_relationship_category/mgi_relationship_property
#
# Usage:
#       mgiMarkerFeed.py
#
# History:
#
# lec	04/25/2016
#	- TR12027/mgi_relationship
#
# lec	01/05/2014
#	- TR11515/allele.bcp/allele_type.bcp/allele_collection.bcp, 
#		allele_subtype.bcp, allele_subtype_assoc.bcp
#
# lec	04/02/2013
#	- TR11343/genotype_mpt : update marker/omim cache query
#
# lec	12/28/2011
#	- changed non-ansi-standard query to left outer join
#
# lec	11/17/2010; add accession ids of cell lines
#	- TR 10460
#	accession_allele_cellline.bcp
#	
#	- TR 9887/remove note duplicate check
#
# lec	10/13/2009
#	- TR 9887; allele notes, check for duplicates
#	this fix can be removed once the alomrkload issue is resolved.
#
# lec	9/15/2009
#	- TR 9838;
#	missing strain_genotype, allele_pair:  
#		no changes; looks like the TR9831 fix solved this problem
#	missing references:  
#		added temp table "allrefs2" to add Allele/Marker references
#		and Allele/Derivation references
#	add "Autoload" to mgi_status()
#
# lec   09/10/2009
#	- TR 9831; add "Autoload" status to allele.bcp
#
# lec   05/17/2009
#       - TR 9405, gene trap less filling (TR7493)
#
#	new vocabularies:
#       allele_transmission.bcp
#	allele_creator.bcp
#       allele_celllinetype.bcp
#       allele_vector.bcp
#       allele_vectortype.bcp
#	allele_markerstatus.bcp
#	allele_markerqualifier.bcp
#       genotype_existsas.bcp
#
#	new joins:
#       allele_cellline.bcp
#	allele_derivation.bcp
#       allele_allele_cellline.bcp
#       allele_marker.bcp
#
#	modified:
#	allele.bcp:
#	genotype.bcp:  _ExistsAs_key
#
# lec	04/23/2008
#	- TR 8511; strain types changed to strain attribute
#
# lec	07/07/2006
#	- add OMIM annotations
#
# lec	03/21/2006	TR 7582
#	- exclude strain synonyms with 'nm####' nomenclature
#	- backed out per Carolyn Blake
#
# lec	10/04/2005
#	- TR 5188; GO Qualifier
#
# lec	08/2005
#	- added mp_term, mp_closure per csb
#	- added header id and header order to genotype_mpt
#
# lec	08/05/2005
#	- added mp_closure
#
# lec	06/01/2005
#	- added all markers (not just mouse) per csb
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
import reportlib
import mgi_utils
import db

db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE(False)

TAB = reportlib.TAB
CRT = reportlib.CRT

COLDELIM = "&=&"
LINEDELIM = "\n#=#"

OUTPUTDIR = os.environ['REPORTOUTPUTDIR'] + '/mgimarkerfeed/'

def strip_newline(s):

	if string.find (s, '\\n') > -1:
		s = string.join(string.split(s, '\\n'), '')

	if string.find(s, '\012') > -1:
		s = string.join(string.split(s, '\012'), '')

	if string.find (s, '\\t') > -1:
		s = string.join(string.split(s, '\\t'), '')

	return s

def mgi_status(status):

	if status == "Approved":
		return "A"

	if status == "Autoload":
		return "U"

	if status == "official":
		return "O"

	if status == "withdrawn":
		return "W"

	return ''

def vocabs():

    #
    # marker_type
    #

    fp = open(OUTPUTDIR + 'marker_type.bcp', 'w')

    results = db.sql('''
	    select _Marker_Type_key, name, 
		to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
		to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from MRK_Types
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Marker_Type_key']` + TAB + \
		     r['name'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # species
    #

    fp = open(OUTPUTDIR + 'species.bcp', 'w')

    results = db.sql('''
	    select _Organism_key, commonName, latinName, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from MGI_Organism_Marker_View
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Organism_key']` + TAB + \
		     r['commonName'] + TAB + \
		     r['latinName'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # allele_type
    #

    fp = open(OUTPUTDIR + 'allele_type.bcp', 'w')

    results = db.sql('''
	    select _Term_key, term, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from VOC_Term where _Vocab_key = 38
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # allele_inheritance_mode
    #
    
    fp = open(OUTPUTDIR + 'allele_inheritance_mode.bcp', 'w')

    results = db.sql('''
	    select _Term_key, term, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from VOC_Term where _Vocab_key = 35
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # allele_pairstate
    #
    
    fp = open(OUTPUTDIR + 'allele_pairstate.bcp', 'w')

    results = db.sql('''
	    select _Term_key, term, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from VOC_Term where _Vocab_key = 39
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # allele_transmission.bcp
    #

    fp = open(OUTPUTDIR + 'allele_transmission.bcp', 'w')

    results = db.sql('''
	    select _Term_key, term, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from VOC_Term where _Vocab_key = 61
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # allele_collection.bcp
    #

    fp = open(OUTPUTDIR + 'allele_collection.bcp', 'w')

    results = db.sql('''
	    select _Term_key, term, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from VOC_Term where _Vocab_key = 92
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # allele_subtype.bcp
    #

    fp = open(OUTPUTDIR + 'allele_subtype.bcp', 'w')

    results = db.sql('''
	    select _Term_key, term, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from VOC_Term where _Vocab_key = 93
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # allele_creator.bcp
    #

    fp = open(OUTPUTDIR + 'allele_creator.bcp', 'w')

    results = db.sql('''
	    select _Term_key, term, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from VOC_Term where _Vocab_key = 62
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # allele_celllinetype.bcp
    #

    fp = open(OUTPUTDIR + 'allele_celllinetype.bcp', 'w')

    results = db.sql('''
	    select _Term_key, term, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from VOC_Term where _Vocab_key = 63
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # allele_vector.bcp
    #

    fp = open(OUTPUTDIR + 'allele_vector.bcp', 'w')

    results = db.sql('''
	    select _Term_key, term, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from VOC_Term where _Vocab_key = 72
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # allele_vectortype.bcp
    #

    fp = open(OUTPUTDIR + 'allele_vectortype.bcp', 'w')

    results = db.sql('''
	    select _Term_key, term, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from VOC_Term where _Vocab_key = 64
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # genotype_existsas.bcp
    #

    fp = open(OUTPUTDIR + 'genotype_existsas.bcp', 'w')

    results = db.sql('''
	    select _Term_key, term, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from VOC_Term where _Vocab_key = 60
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
		     r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # strain_type.bcp
    #
    
    fp = open(OUTPUTDIR + 'strain_type.bcp', 'w')
    
    results = db.sql('''
	  select _Term_key, term, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
          from VOC_Term where _Vocab_key = 27
	  ''', 'auto')

    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
	             r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # mp_term.bcp
    #
    
    fp = open(OUTPUTDIR + 'mp_term.bcp', 'w')
    
    results = db.sql('''
	  select x._Term_key, x.sequenceNum, rtrim(x.note) as note
          from VOC_Term t, VOC_Text x 
	  where t._Vocab_key = 5 
	  and t._Term_key = x._Term_key 
	  order by x._Term_key, x.sequenceNum
	  ''', 'auto')

    notes = {}
    for r in results:
	key = r['_Term_key']
	value = r['note']
	if not notes.has_key(key):
	    notes[key] = []
	notes[key].append(value)

    #
    # vocabulary terms
    # MP (5) 
    # OMIM (44)
    # Disease Ontology (125)
    #

    results = db.sql('''
	  select a.accID, t.term, t._Term_key, v.name, 
                to_char(t.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(t.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
          from VOC_Vocab v, VOC_Term t, ACC_Accession a 
	  where t._Vocab_key in (5, 44, 125) 
	  and t._Vocab_key = v._Vocab_key 
	  and t._Term_key = a._Object_key 
	  and a._MGIType_key = 13 
	  and a.preferred = 1
	  ''', 'auto')

    for r in results:
	    key = r['_Term_key']

	    fp.write(`key` + TAB + \
		     r['accID'] + TAB + \
	             r['term'] + TAB + \
		     r['name'] + TAB)

            if notes.has_key(key):
		fp.write(string.join(notes[key], ''))
            fp.write(TAB)

	    fp.write(str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # mp_synonym.bcp
    # Synonyms for MP and OMIM terms
    #

    fp = open(OUTPUTDIR + 'mp_synonym.bcp', 'w')
    
    results = db.sql('''
	  select t._Term_key, s.synonym, 
                to_char(s.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(s.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
          from VOC_Term t, MGI_Synonym s 
	  where t._Vocab_key in (5, 44) 
	  and t._Term_key = s._Object_key 
	  and s._MGIType_key = 13 
	  ''', 'auto')

    for r in results:
	    key = r['_Term_key']

	    fp.write(`key` + TAB + \
		     r['synonym'] + TAB + \
	             str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # mp_closure.bcp
    #
    
    fp = open(OUTPUTDIR + 'mp_closure.bcp', 'w')

    results = db.sql('''
	select c._AncestorObject_key, c._DescendentObject_key, 
		c._AncestorLabel_key, c._DescendentLabel_key 
	from DAG_Closure c 
	where c._Dag_key = 4
	''', 'auto')

    for r in results:
	fp.write(`r['_AncestorObject_key']` + TAB + \
		 `r['_DescendentObject_key']` + TAB + \
		 `r['_AncestorLabel_key']` + TAB + \
		 `r['_DescendentLabel_key']` + CRT)

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

    db.sql('''select m._Marker_key, m._Organism_key, m.symbol, m.name 
	    into temporary table markers 
	    from MRK_Marker m 
	    where m._Organism_key = 1 
	    and exists (select 1 from ACC_Accession a 
	    where m._Marker_key = a._Object_key 
	    and a._MGIType_key = 2 
	    and a.prefixPart = 'MGI:' 
	    and a._LogicalDB_key = 1 
	    and a.preferred = 1) 
	    union 
	    select m._Marker_key, m._Organism_key, m.symbol, m.name 
	    from MRK_Marker m 
	    where m._Organism_key != 1
	    ''', None)

    db.sql('create index markers_idx1 on markers(_Marker_key)', None)

    #
    # select data fields for marker.bcp
    #

    fp = open(OUTPUTDIR + 'marker.bcp', 'w')

    results = db.sql('''
	    select k._Marker_key, k._Organism_key, m._Marker_Type_key, s.status, 
	    m.symbol, m.name, m.chromosome, m.cytogeneticOffset, o.cmoffset, 
            to_char(m.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
            to_char(m.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from markers k, MRK_Marker m, MRK_Status s, MRK_Offset o 
	    where k._Marker_key = m._Marker_key 
	    and m._Marker_Status_key = s._Marker_Status_key 
	    and k._Marker_key = o._Marker_key 
	    and o.source = 0 
	    union 
            select m._Marker_key, m._Organism_key, m._Marker_Type_key, s.status, 
	    m.symbol, m.name, m.chromosome, m.cytogeneticOffset,  null, 
            to_char(m.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
            to_char(m.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from markers k, MRK_Marker m, MRK_Status s 
	    where k._Organism_key != 1 
	    and k._Marker_key = m._Marker_key 
	    and m._Marker_Status_key = s._Marker_Status_key
	    ''', 'auto')

    for r in results:
	    fp.write(`r['_Marker_key']` + TAB + \
	             `r['_Organism_key']` + TAB + \
	             `r['_Marker_Type_key']` + TAB + \
		     mgi_status(r['status']) + TAB + \
		     strip_newline(r['symbol']) + TAB + \
		     strip_newline(r['name']) + TAB + \
		     r['chromosome'] + TAB + \
		     mgi_utils.prvalue(r['cytogeneticOffset']) + TAB + \
		     mgi_utils.prvalue(r['cmoffset']) + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)

    fp.close()

    #
    # select data fields for marker_label.bcp
    # status = 1 if the label matches the MRK_Marker.symbol or MRK_Marker.name
    # field, else status = 0
    #

    fp = open(OUTPUTDIR + 'marker_label.bcp', 'w')

    results = db.sql('''
	    select distinct m._Marker_key, m.label, m.labelType, 1 as status, 
            to_char(current_date, 'Mon DD YYYY HH:MIAM') as cdate
	    from markers k, MRK_Label m 
	    where k._Marker_key = m._Marker_key 
	    and m.labelType = 'MS'
	    and k.symbol = m.label 
	    ''', 'auto')

    for r in results:
	    fp.write(`r['_Marker_key']` + TAB + \
		     strip_newline(r['label']) + TAB + \
		     r['labelType'] + TAB + \
		     `r['status']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['cdate']) + CRT)

    results = db.sql('''
	    select distinct m._Marker_key, m.label, m.labelType, 1 as status, 
                to_char(current_date, 'Mon DD YYYY HH:MIAM') as cdate
	    from markers k, MRK_Label m 
	    where k._Marker_key = m._Marker_key 
	    and m.labelType = 'MN' 
	    and k.name = m.label
	    ''', 'auto')

    for r in results:
	    fp.write(`r['_Marker_key']` + TAB + \
		     strip_newline(r['label']) + TAB + \
		     r['labelType'] + TAB + \
		     `r['status']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['cdate']) + CRT)

    results = db.sql('''
	    select distinct m._Marker_key, m.label, m.labelType, 0 as status, 
                to_char(current_date, 'Mon DD YYYY HH:MIAM') as cdate
	    from markers k, MRK_Label m 
	    where k._Marker_key = m._Marker_key 
	    and m.labelType in ('MS', 'MN', 'MY') 
	    and k.symbol != m.label and k.name != m.label 
	    ''', 'auto')

    for r in results:
	    fp.write(`r['_Marker_key']` + TAB + \
		     strip_newline(r['label']) + TAB + \
		     r['labelType'] + TAB + \
		     `r['status']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['cdate']) + CRT)

    fp.close()

    #
    # select data fields for accession_marker.bcp
    # select MGI Accession IDs only
    #

    fp = open(OUTPUTDIR + 'accession_marker.bcp', 'w')

    results = db.sql('''
	    select m.accID, l.name as LogicalDB, m._Object_key, m.preferred, 
                to_char(m.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(m.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from markers k, ACC_Accession m, ACC_LogicalDB l 
	    where k._Marker_key = m._Object_key 
	    	and m._MGIType_key = 2 
	    	and m.prefixPart = 'MGI:' 
	    	and m._LogicalDB_key = 1 
	    	and m._LogicalDB_key = l._LogicalDB_key 
	    union 
            select m.accID, l.name as LogicalDB, m._Object_key, m.preferred, 
                to_char(m.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(m.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from markers k, ACC_Accession m, ACC_LogicalDB l 
	    where k._Organism_key != 1 
	    	and k._Marker_key = m._Object_key 
	    	and m._MGIType_key = 2 
	    	and m._LogicalDB_key in (55, 64, 15, 47) 
	    	and m._LogicalDB_key = l._LogicalDB_key 
	    ''', 'auto')

    for r in results:
	    fp.write(r['accID'] + TAB + \
		     r['LogicalDB'] + TAB + \
		     `r['_Object_key']` + TAB + \
		     `r['preferred']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

def alleles():

    #
    # allele_cellline
    # accession_allele_cellline
    # allele
    # allele_allele_cellline
    # allele_label
    # allele_pair
    # allele_note
    # accession_allele
    # allele_subtype_assoc
    #

    #
    # allele_cellline.bcp
    #

    fp = open(OUTPUTDIR + 'allele_cellline.bcp', 'w')
    
    results = db.sql('''
	    select _CellLine_key, cellLine, _CellLine_Type_key, _Strain_key, _Derivation_key, isMutant, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from ALL_CellLine
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_CellLine_key']` + TAB + \
		     r['cellLine'] + TAB + \
		     `r['_CellLine_Type_key']` + TAB + \
		     `r['_Strain_key']` + TAB + \
		     `r['_Derivation_key']` + TAB + \
		     `r['isMutant']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # accession_allele_cellline.bcp
    #

    fp = open(OUTPUTDIR + 'accession_allele_cellline.bcp', 'w')

    results = db.sql('''
	    select m.accID, l.name as LogicalDB, m._Object_key, m.preferred, 
                to_char(m.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(m.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from ALL_Cellline c, ACC_Accession m, ACC_LogicalDB l 
	    where c._CellLine_key = m._Object_key 
	    and m._MGIType_key = 28 
	    and m._LogicalDB_key = l._LogicalDB_key
	    ''', 'auto')

    for r in results:
	    fp.write(r['accID'] + TAB + \
		     r['LogicalDB'] + TAB + \
		     `r['_Object_key']` + TAB + \
		     `r['preferred']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # allele_derivation.bcp
    #

    fp = open(OUTPUTDIR + 'allele_derivation.bcp', 'w')
    
    results = db.sql('''
	select _Derivation_key, name, 
	    _Vector_key, _VectorType_key, 
	    _ParentCellLine_key, _DerivationType_key, 
	    _Creator_key, _Refs_key, 
            to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
            to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from ALL_CellLine_Derivation
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Derivation_key']` + TAB + \
		     r['name'] + TAB + \
		     `r['_Vector_key']` + TAB + \
		     `r['_VectorType_key']` + TAB + \
		     `r['_ParentCellLine_key']` + TAB + \
		     `r['_DerivationType_key']` + TAB + \
		     `r['_Creator_key']` + TAB + \
		     `r['_Refs_key']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # select all alleles with a status of 'approved' or 'autoload'
    # all other statuses are private/confidential alleles
    #

    db.sql('''
	select m._Allele_key into temporary table alleles 
	from ALL_Allele m, VOC_Term t 
	where m._Allele_Status_key = t._Term_key 
	and t.term in ('Approved', 'Autoload') 
	''', None)
    db.sql('create index alleles_idx1 on alleles(_Allele_key)', None)

    #
    # select data fields for allele.bcp
    #

    fp = open(OUTPUTDIR + 'allele.bcp', 'w')

    results = db.sql('''
	select m._Allele_key, m._Marker_key, m._Mode_key, m._Allele_Type_key, m._Transmission_key, 
	m._Strain_key, m._Collection_key, m.isWildType, m.isExtinct, m.isMixed, 
	m.symbol, m.name, t.term as status, 
        to_char(m.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
        to_char(m.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	from alleles a, ALL_Allele m, VOC_Term t 
	where a._Allele_key = m._Allele_key 
	and m._Allele_Status_key = t._Term_key
	''', 'auto')

    for r in results:
	    fp.write(`r['_Allele_key']` + TAB + \
	             `r['_Marker_key']` + TAB + \
	             `r['_Mode_key']` + TAB + \
	             `r['_Allele_Type_key']` + TAB + \
	             `r['_Transmission_key']` + TAB + \
		     mgi_status(r['status']) + TAB + \
		     `r['_Strain_key']` + TAB + \
	             `r['_Collection_key']` + TAB + \
		     `r['isWildType']` + TAB + \
		     `r['isExtinct']` + TAB + \
		     `r['isMixed']` + TAB + \
		     strip_newline(r['symbol']) + TAB + \
		     strip_newline(r['name']) + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # select data fields for allele_allele_cellline.bcp
    #

    fp = open(OUTPUTDIR + 'allele_allele_cellline.bcp', 'w')

    results = db.sql('''
	    select m._Assoc_key, m._Allele_key, m._MutantCellLine_key, 
                to_char(m.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(m.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from alleles a, ALL_Allele_Cellline m 
	    where a._Allele_key = m._Allele_key 
	    ''', 'auto')

    for r in results:
	    fp.write(`r['_Assoc_key']` + TAB + \
		     `r['_Allele_key']` + TAB + \
		     `r['_MutantCellLine_key']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # select data fields for allele_label.bcp
    # this includes all allele symbols, names and synonyms
    # status = 1 if the label matches the ALL_Allele.symbol or ALL_Allele.name
    # field, else status = 0
    #

    fp = open(OUTPUTDIR + 'allele_label.bcp', 'w')

    results = db.sql('''
	    select m._Allele_key, m.name, 'AN' as labelType, 1 as status, 
                to_char(m.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(m.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from alleles a, ALL_Allele m 
	    where a._Allele_key = m._Allele_key 
	    union 
	    select m._Allele_key, m.symbol, 'AS' as labelType, 1 as status, 
                to_char(m.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(m.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from alleles a, ALL_Allele m 
	    where a._Allele_key = m._Allele_key 
	    union 
	    select m._Allele_key, s.synonym, 'AY' as labelType, 0 as status, 
                to_char(m.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(m.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from alleles a, ALL_Allele m, MGI_Synonym s 
	    where a._Allele_key = m._Allele_key 
	    and a._Allele_key = s._Object_key 
	    and s._MGIType_key = 11
	    ''', 'auto')

    for r in results:
	    fp.write(`r['_Allele_key']` + TAB + \
		     strip_newline(r['name']) + TAB + \
		     r['labelType'] + TAB + \
		     `r['status']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # select data fields for allele_note.bcp
    #

    fp = open(OUTPUTDIR + 'allele_note.bcp', 'w')

    results = db.sql('''
	    select a._Allele_key, rtrim(nc.note) as note, nc.sequenceNum, nt.noteType, 
                to_char(n.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(n.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from alleles a, MGI_Note n, MGI_NoteChunk nc, MGI_NoteType nt 
	    where a._Allele_key = n._Object_key 
	    and n._MGIType_key = 11 
	    and n._Note_key = nc._Note_key 
	    and n._NoteType_key = nt._NoteType_key
	    and nt.private = 0
	    ''', 'auto')

    notes = {}
    for r in results:
	    fp.write(`r['_Allele_key']` + COLDELIM + \
		     r['noteType'] + COLDELIM + \
		     `r['sequenceNum']` + COLDELIM + \
		     mgi_utils.prvalue(r['note']) + COLDELIM + \
		     str(r['cdate']) + COLDELIM + \
		     str(r['mdate']) + LINEDELIM)
    fp.close()

    #
    # select data fields for accession_allele.bcp
    #

    fp = open(OUTPUTDIR + 'accession_allele.bcp', 'w')

    results = db.sql('''
	    select m.accID, l.name as LogicalDB, m._Object_key, m.preferred, 
                to_char(m.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(m.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from alleles k, ACC_Accession m, ACC_LogicalDB l 
	    where k._Allele_key = m._Object_key 
	    and m._MGIType_key = 11 
	    and m.prefixPart = 'MGI:' 
	    and m._LogicalDB_key = 1 
	    and m._LogicalDB_key = l._LogicalDB_key
	    ''', 'auto')

    for r in results:
	    fp.write(r['accID'] + TAB + \
		     r['LogicalDB'] + TAB + \
		     `r['_Object_key']` + TAB + \
		     `r['preferred']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # select data fields for allele_subtype_assoc.bcp
    #

    fp = open(OUTPUTDIR + 'allele_subtype_assoc.bcp', 'w')

    results = db.sql('''
	    select a._Allele_key, va._Term_key, 
                to_char(va.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(va.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from alleles a, VOC_Annot va
	    where a._Allele_key = va._Object_key
	    and va._AnnotType_key = 1014
	    ''', 'auto')

    for r in results:
	    fp.write(`r['_Allele_key']` + TAB + \
		     `r['_Term_key']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
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

    db.sql('''
	  select distinct s._Strain_key, s._Species_key, s._StrainType_key, s.strain, s.private, 
                to_char(s.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(s.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
          into temporary table strains 
          from PRB_Strain s, ACC_Accession a 
          where s._Strain_key = a._Object_key 
          and a._MGIType_key = 10 
          and a._LogicalDB_key in (22, 38) 
          union 
          select distinct s._Strain_key, s._Species_key, s._StrainType_key, s.strain, s.private, 
                to_char(s.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(s.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
          from PRB_Strain s, ALL_Allele a 
	  where s._Strain_key = a._Strain_key 
          union 
          select distinct s._Strain_key, s._Species_key, s._StrainType_key, s.strain, s.private, 
                to_char(s.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(s.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
          from PRB_Strain s, ALL_CellLine a 
	  where s._Strain_key = a._Strain_key
	  ''', None)

    db.sql('create index strains_idx1 on strains(_Strain_key)', None)

    #
    # strain.bcp
    #

    fp = open(OUTPUTDIR + 'strain.bcp', 'w')

    results = db.sql('select * from strains', 'auto')

    for r in results:
	    fp.write(`r['_Strain_key']` + TAB + \
	             `r['_Species_key']` + TAB + \
	             strip_newline(r['strain']) + TAB + \
	             `r['private']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # strain_marker.bcp
    #

    fp = open(OUTPUTDIR + 'strain_marker.bcp', 'w')

    results = db.sql('''
	  select distinct m._Strain_key, m._Marker_key, m._Allele_key, s.private, t.term as qualifier, 
                to_char(m.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(m.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
          from strains s, PRB_Strain_Marker m, VOC_Term t 
          where s._Strain_key = m._Strain_key 
	  and m._Qualifier_key = t._Term_key
	  ''', 'auto')

    for r in results:
	    fp.write(`r['_Strain_key']` + TAB + \
	             `r['_Marker_key']` + TAB + \
	             mgi_utils.prvalue(r['_Allele_key']) + TAB + \
	             `r['private']` + TAB + \
		     r['qualifier'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # strain_synonym.bcp
    #

    fp = open(OUTPUTDIR + 'strain_synonym.bcp', 'w')
    
    results = db.sql('''
	  select m._Synonym_key, m._Object_key, m.synonym, s.private, m.synonymType, 
                to_char(m.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(m.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
          from strains s, MGI_Synonym_Strain_View m 
          where s._Strain_key = m._Object_key
	  ''', 'auto')
#	  and m.synonym not like 'nm[0-9]%' ', 'auto')

    for r in results:
	    fp.write(`r['_Synonym_key']` + TAB + \
	             `r['_Object_key']` + TAB + \
	             `r['private']` + TAB + \
	             strip_newline(r['synonym']) + TAB + \
	             strip_newline(r['synonymType']) + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # strain_strain_type.bcp
    #

    fp = open(OUTPUTDIR + 'strain_strain_type.bcp', 'w')

    results = db.sql('''
	select m._Strain_key, m._Term_key, s.private, s.cdate, s.mdate 
        from strains s, PRB_Strain_Attribute_View m 
        where s._Strain_key = m._Strain_key
	''', 'auto')

    for r in results:
	    fp.write(`r['_Strain_key']` + TAB + \
	             `r['_Term_key']` + TAB + \
	             `r['private']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)

    fp.close()

    #
    # accession_strain.bcp
    #

    fp = open(OUTPUTDIR + 'accession_strain.bcp', 'w')
    
    results = db.sql('''
	select distinct a.accID, l.name as LogicalDB, a._Object_key, a.preferred, s.private, 
               to_char(a.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
               to_char(a.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
          from strains s, ACC_Accession a, ACC_LogicalDB l 
          where s._Strain_key = a._Object_key 
          and a._MGIType_key = 10 
          and a._LogicalDB_key = l._LogicalDB_key
	  ''', 'auto')

    for r in results:
	    fp.write(r['accID'] + TAB + \
		     r['LogicalDB'] + TAB + \
	             `r['_Object_key']` + TAB + \
	             `r['preferred']` + TAB + \
	             `r['private']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()
    
    #
    # strain_species.bcp
    #
    
    fp = open(OUTPUTDIR + 'strain_species.bcp', 'w')

    results = db.sql('''
	select _Term_key, term, 
                to_char(creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(modification_date, 'Mon DD YYYY HH:MIAM') as mdate
          from VOC_Term where _Vocab_key = 26
	  ''', 'auto')
    
    for r in results:
	    fp.write(`r['_Term_key']` + TAB + \
	             r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

def genotypes():

    #
    # genotype.bcp
    # genotype_mpt.bcp
    # genotype_header.bcp
    # strain_genotype.bcp
    # allele_pair.bcp
    #

    #
    # MP/Genotype annotations (1002)
    # OMIM/Genotype annotations (1005)
    # DO/Genotype annotations (1020)
    #

    db.sql('''
	   select distinct g._Genotype_key 
	   into temporary table genotypes 
	   from strains s, GXD_Genotype g, VOC_Annot a 
	   where s._Strain_key = g._Strain_key 
	   and g._Genotype_key = a._Object_key 
	   and a._AnnotType_key in (1002,1005,1020) 
	   union 
	   select distinct g._Genotype_key 
	   from strains s, PRB_Strain_Genotype g 
	   where s._Strain_key = g._Strain_key
	   ''', 'auto')

    db.sql('create index genotypes_idx1 on genotypes(_Genotype_key)', None)

    #
    # genotype.bcp
    #

    fp = open(OUTPUTDIR + 'genotype.bcp', 'w')

    results = db.sql('''
	  select g._Genotype_key, s.strain, p.isConditional, p._ExistsAs_key, c.note, 
                to_char(p.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(p.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	  from genotypes g, GXD_Genotype p, PRB_Strain s, MGI_Note n, MGI_NoteChunk c 
	  where g._Genotype_key = p._Genotype_key 
	  and p._Strain_key = s._Strain_key 
	  and g._Genotype_key = n._Object_key 
	  and n._NoteType_key = 1016 
	  and n._Note_key = c._Note_key
	  ''', 'auto')

    for r in results:

	note = string.replace(string.strip(r['note']), '\n', ' ')

	fp.write(`r['_Genotype_key']` + TAB + \
		r['strain'] + TAB + \
		`r['isConditional']` + TAB + \
	        `r['_ExistsAs_key']` + TAB + \
		note + TAB + \
		str(r['cdate']) + TAB + \
		str(r['mdate']) + CRT)
    fp.close()

    #
    # genotype_mpt.bcp
    #
    # select all genotypes
    #
    # for those that contain marker/DO annotations (mrk_do_cache)
    # include the omimCatergory3 ('None' if no DO annotation exists)
    #

    fp = open(OUTPUTDIR + 'genotype_mpt.bcp', 'w')

    results = db.sql('''
	(
	select distinct g._Genotype_key, a._AnnotType_key, a._Annot_key, a._Term_key, 
	q.term as qualifier, c.omimCategory3,
                to_char(a.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(a.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	from genotypes g
	INNER JOIN VOC_Annot a on (g._Genotype_key = a._Object_key and a._AnnotType_key in (1002,1005,1020))
	INNER JOIN VOC_Term q on (a._Qualifier_key = q._Term_key)
	LEFT OUTER JOIN MRK_DO_Cache c on (g._Genotype_key = c._Genotype_key and a._Term_key = c._Term_key)
	)
	order by _Genotype_key
	''', 'auto')

    for r in results:

	if r['omimCategory3'] == None:
	    category = '-1'
	else:
	    category = r['omimCategory3']

	fp.write(`r['_Annot_key']` + TAB + \
		 `r['_AnnotType_key']` + TAB + \
	         `r['_Term_key']` + TAB + \
	         `r['_Genotype_key']` + TAB + \
		string.strip(mgi_utils.prvalue(r['qualifier'])) + TAB + \
		str(category) + TAB + \
		str(r['cdate']) + TAB + \
		str(r['mdate']) + CRT)

    fp.close()

    #
    # genotype_header.bcp (MP only)
    #

    fp = open(OUTPUTDIR + 'genotype_header.bcp', 'w')

    results = db.sql('''
	select g._Genotype_key, h._Term_key as headerTerm, h.sequenceNum, 
                to_char(h.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(h.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	from genotypes g, VOC_AnnotHeader h 
	where g._Genotype_key = h._Object_key 
	and h._AnnotType_key = 1002 
	''', 'auto')

    for r in results:
	fp.write(`r['_Genotype_key']` + TAB + \
	         `r['headerTerm']` + TAB + \
	         `r['sequenceNum']` + TAB + \
		str(r['cdate']) + TAB + \
		str(r['mdate']) + CRT)

    fp.close()


    #
    # strain_genotype.bcp
    #
    
    fp = open(OUTPUTDIR + 'strain_genotype.bcp', 'w')

    results = db.sql('''
	  select p._Strain_key, p._Genotype_key, t.term, 
                to_char(p.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(p.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
          from strains s, genotypes g, PRB_Strain_Genotype p, VOC_Term t 
          where s._Strain_key = p._Strain_key 
	  and g._Genotype_key = p._Genotype_key 
          and p._Qualifier_key = t._Term_key
	  ''', 'auto')
    
    for r in results:
            fp.write(`r['_Strain_key']` + TAB + \
                     `r['_Genotype_key']` + TAB + \
                     r['term'] + TAB + \
                     str(r['cdate']) + TAB + \
                     str(r['mdate']) + CRT)
    fp.close()

    #
    # allele_pair.bcp
    #

    fp = open(OUTPUTDIR + 'allele_pair.bcp', 'w')

    results = db.sql('''
	    (
            select p._AllelePair_key, p._Allele_key_1, p._Allele_key_2, 
            p._Marker_key, p._Genotype_key, p._PairState_key, p.sequenceNum,
            to_char(p.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
            to_char(p.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from alleles a, genotypes g, GXD_AllelePair p 
	    where a._Allele_key = p._Allele_key_1 
	    and g._Genotype_key = p._Genotype_key 
	    union 
            select p._AllelePair_key, p._Allele_key_1, p._Allele_key_2, 
            p._Marker_key, p._Genotype_key, p._PairState_key, p.sequenceNum,
            to_char(p.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
            to_char(p.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from alleles a, genotypes g, GXD_AllelePair p 
	    where a._Allele_key = p._Allele_key_2 
	    and g._Genotype_key = p._Genotype_key 
	    )
	    order by _Genotype_key, sequenceNum
	    ''', 'auto')

    for r in results:
	    fp.write(`r['_AllelePair_key']` + TAB + \
		     `r['_Allele_key_1']` + TAB + \
		     mgi_utils.prvalue(r['_Allele_key_2']) + TAB + \
		     `r['_Marker_key']` + TAB + \
		     `r['_Genotype_key']` + TAB + \
		     `r['_PairState_key']` + TAB + \
		     `r['sequenceNum']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

def references():

    #
    # reference.bcp
    #

    fp = open(OUTPUTDIR + 'reference.bcp', 'w')

    #
    # references annotated to a Genotype
    #

    db.sql('''
	   select distinct e._Refs_key, g._Genotype_key, a._Annot_key, e._AnnotEvidence_key, 
           to_char(e.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
           to_char(e.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	   into temporary table genoreferences 
	   from genotypes g, VOC_Annot a, VOC_Evidence e 
	   where g._Genotype_key = a._Object_key 
	   and a._AnnotType_key in (1002,1005,1020) 
	   and a._Annot_key = e._Annot_key
	   ''', None)

    #
    # references annotated to an Allele via the Genotype
    # references annotated to an Allele via the Strain
    #

    db.sql('''
	   select distinct r._Refs_key, ag._Allele_key, rt.assocType, 
           to_char(r.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
           to_char(r.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	   into temporary table allrefs 
	   from genotypes g, GXD_AlleleGenotype ag, MGI_Reference_Assoc r, MGI_RefAssocType rt 
	   where g._Genotype_key = ag._Genotype_key 
	   and ag._Allele_key = r._Object_key 
	   and r._MGIType_key = 11 
	   and r._RefAssocType_key = rt._RefAssocType_key 
	   union 
           select distinct r._Refs_key, sm._Allele_key, rt.assocType, 
           to_char(r.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
           to_char(r.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	   from strains s, PRB_Strain_Marker sm, MGI_Reference_Assoc r, MGI_RefAssocType rt 
	   where s._Strain_key = sm._Strain_key 
	   and sm._Allele_key = r._Object_key 
	   and r._MGIType_key = 11 
	   and r._RefAssocType_key = rt._RefAssocType_key
	   ''', None)

    #
    # references annotated to an Allele via the Derivation
    #

    db.sql('''
           select distinct r._Refs_key, 
           to_char(r.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
           to_char(r.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	   into temporary table allrefs2
	   from ALL_CellLine_Derivation r 
	   where r._Refs_key is not null
	   ''', None)

    #
    # references annotated to Strains
    #

    db.sql('''
	   select r._Refs_key, r._Object_key, rt.assocType, 
           to_char(r.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
           to_char(r.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
           into temporary table strainreferences 
           from strains s, MGI_Reference_Assoc r, MGI_RefAssocType rt 
           where s._Strain_key = r._Object_key 
           and r._MGIType_key = 10 
           and r._RefAssocType_key = rt._RefAssocType_key
	   ''', None)

    #
    # references annotated to a Marker via the Genotype
    # references annotated to a Marker via the Strain
    #

    db.sql('''
	   select distinct r._Refs_key, r._Marker_key 
	   into temporary table mrkreferences 
	   from genotypes g, GXD_AlleleGenotype ag, MRK_Reference r 
	   where g._Genotype_key = ag._Genotype_key 
	   and ag._Marker_key = r._Marker_key 
	   union 
           select distinct r._Refs_key, sm._Marker_key 
	   from strains s, PRB_Strain_Marker sm, MGI_Reference_Assoc r 
	   where s._Strain_key = sm._Strain_key 
	   and sm._Marker_key = r._Object_key 
	   and r._MGIType_key = 2 
	   ''', None)

    #
    # references used in Human/OMIM annotations
    # references used in Human/DO annotations
    #

    db.sql('''
	   select distinct e._Refs_key, 
           to_char(e.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
           to_char(e.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	   into temporary table omimreferences 
	   from VOC_Annot a, VOC_Evidence e 
	   where a._AnnotType_key in (1006,1022)
	   and a._Annot_key = e._Annot_key
	   ''', None)

    db.sql('create index genoreferences_idx1 on genoreferences(_Refs_key)', None)
    db.sql('create index allrefs_idx1 on allrefs(_Refs_key)', None)
    db.sql('create index allrefs2_idx1 on allrefs2(_Refs_key)', None)
    db.sql('create index strainreferences_idx1 on strainreferences(_Refs_key)', None)
    db.sql('create index mrkreferences_idx1 on mrkreferences(_Refs_key)', None)
    db.sql('create index omimreferences_idx1 on omimreferences(_Refs_key)', None)

    db.sql('''
	   select distinct _Refs_key into temporary table refs from genoreferences 
	   union select distinct _Refs_key from allrefs 
	   union select distinct _Refs_key from allrefs2 
	   union select distinct _Refs_key from strainreferences 
	   union select distinct _Refs_key from mrkreferences 
	   union select distinct _Refs_key from omimreferences
	   ''', None)
    db.sql('create index references_idx1 on refs(_Refs_key)', None)

    results = db.sql('''
	    select r._Refs_key, b.refType, b.authors,
	    b.title, b.journal, b.vol, b.issue, b.pgs, b.year, 
	    b.isReviewArticle, 
	    k.book_au, k.book_title, k.publisher, k.place, k.series_ed, 
            to_char(b.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
            to_char(b.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from refs r, BIB_Refs b LEFT OUTER JOIN BIB_Books k on (b._Refs_key = k._Refs_key)
	    where r._Refs_key = b._Refs_key 
	    ''', 'auto')

    for r in results:
	    fp.write(`r['_Refs_key']` + TAB + \
	             r['refType'] + TAB + \
		     mgi_utils.prvalue(r['authors']) + TAB + \
		     mgi_utils.prvalue(r['title']) + TAB + \
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
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # accession_reference.bcp
    #

    fp = open(OUTPUTDIR + 'accession_reference.bcp', 'w')
    results = db.sql('''
	    select a.accID, l.name as LogicalDB, a._Object_key, a.preferred, 
            to_char(a.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
            to_char(a.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from refs r, ACC_Accession a, ACC_LogicalDB l 
	    where r._Refs_key = a._Object_key 
	    and a._MGIType_key = 1 
	    and a._LogicalDB_key = l._LogicalDB_key
	    ''', 'auto')

    for r in results:
	    fp.write(r['accID'] + TAB + \
		     r['LogicalDB'] + TAB + \
		     `r['_Object_key']` + TAB + \
		     `r['preferred']` + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # genotype_mpt_reference.bcp
    #

    fp = open(OUTPUTDIR + 'genotype_mpt_reference.bcp', 'w')

    results = db.sql('''
	select g._Refs_key, g._Genotype_key, g._Annot_key, g.cdate, g.mdate 
	from genoreferences g 
	where not exists (select 1 from MGI_Note_VocEvidence_View n 
	where g._AnnotEvidence_key = n._Object_key) 
	order by g._Genotype_key 
	''', 'auto')
    for r in results:
	    fp.write(`r['_Annot_key']` + COLDELIM + \
		     `r['_Refs_key']` + COLDELIM + \
		     COLDELIM + \
		     COLDELIM + \
		     COLDELIM + \
		     COLDELIM + \
		     str(r['cdate']) + COLDELIM + \
		     str(r['mdate']) + LINEDELIM)

    results = db.sql('''
	select g._Refs_key, g._Genotype_key, g._Annot_key, g._AnnotEvidence_key, g.cdate, g.mdate, 
	n._Note_key, n.noteType, n.sequenceNum, n.note 
	from genoreferences g, MGI_Note_VocEvidence_View n 
	where g._AnnotEvidence_key = n._Object_key 
	order by g._Genotype_key, g._AnnotEvidence_key, n.sequenceNum
	''', 'auto')
    for r in results:
	    fp.write(`r['_Annot_key']` + COLDELIM + \
		     `r['_Refs_key']` + COLDELIM + \
		     `r['_Note_key']` + COLDELIM + \
		     r['noteType'] + COLDELIM + \
		     `r['sequenceNum']` + COLDELIM + \
		     r['note'] + COLDELIM + \
		     str(r['cdate']) + COLDELIM + \
		     str(r['mdate']) + LINEDELIM)
    fp.close()

    #
    # allele_reference.bcp
    #

    fp = open(OUTPUTDIR + 'allele_reference.bcp', 'w')

    results = db.sql('select * from allrefs', 'auto')
    for r in results:
	    fp.write(`r['_Allele_key']` + TAB + \
		     `r['_Refs_key']` + TAB + \
		     r['assocType'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

    #
    # strain_reference.bcp
    #

    fp = open(OUTPUTDIR + 'strain_reference.bcp', 'w')

    results = db.sql('select * from strainreferences', 'auto')
    for r in results:
            fp.write(`r['_Object_key']` + TAB + \
                     `r['_Refs_key']` + TAB + \
                     r['assocType'] + TAB + \
                     str(r['cdate']) + TAB + \
                     str(r['mdate']) + CRT)
    fp.close()

    #
    # marker_reference.bcp
    #

    fp = open(OUTPUTDIR + 'marker_reference.bcp', 'w')

    results = db.sql('select * from mrkreferences', 'auto')
    for r in results:
	    fp.write(`r['_Marker_key']` + TAB + \
		     `r['_Refs_key']` + CRT)
    fp.close()

def omim():

    #
    # Human Marker/OMIM annotations (1006)
    # Human Marker/DO annotations (1022)
    #
    # marker_omim.bcp
    #

    fp = open(OUTPUTDIR + 'marker_omim.bcp', 'w')

    results = db.sql('''
	select a._Term_key, a._Object_key, e._Refs_key, 
               to_char(a.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
               to_char(a.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	from VOC_Annot a, VOC_Evidence e 
	where a._AnnotType_key in (1006,1022)
	and a._Annot_key = e._Annot_key
	''', 'auto')

    for r in results:
	fp.write(`r['_Term_key']` + TAB + \
	         `r['_Object_key']` + TAB + \
	         `r['_Refs_key']` + TAB + \
		 str(r['cdate']) + TAB + \
		 str(r['mdate']) + CRT)

    fp.close()

def mgi_relationship():

    #
    # mgi_relationship
    # mgi_relationship_category
    # mgi_relationship_property
    #
    # selecting info for 'mutation_involves' (1003)
    # and 'expresses_component' (1004)
    #

    fp = open(OUTPUTDIR + 'mgi_relationship.bcp', 'w')

    results = db.sql('''
    	select a.*,
	       to_char(a.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
	       to_char(a.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	from mgi_relationship a 
	where a._category_key in (1003,1004)
	''', 'auto')

    for r in results:
	fp.write(`r['_relationship_key']` + TAB + \
	         `r['_category_key']` + TAB + \
	         `r['_object_key_1']` + TAB + \
	         `r['_object_key_2']` + TAB + \
	         `r['_relationshipterm_key']` + TAB + \
	         `r['_qualifier_key']` + TAB + \
	         `r['_evidence_key']` + TAB + \
	         `r['_refs_key']` + TAB + \
		 str(r['cdate']) + TAB + \
		 str(r['mdate']) + CRT)

    fp.close()

    fp = open(OUTPUTDIR + 'mgi_relationship_category.bcp', 'w')

    results = db.sql('''
    	select a.*,
	       to_char(a.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
	       to_char(a.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	from mgi_relationship_category a
	where a._category_key in (1003,1004)
	''', 'auto')

    for r in results:
	fp.write(`r['_category_key']` + TAB + \
	         r['name'] + TAB + \
	         `r['_relationshipvocab_key']` + TAB + \
	         `r['_relationshipdag_key']` + TAB + \
	         `r['_mgitype_key_1']` + TAB + \
	         `r['_mgitype_key_2']` + TAB + \
	         `r['_qualifiervocab_key']` + TAB + \
	         `r['_evidencevocab_key']` + TAB + \
		 str(r['cdate']) + TAB + \
		 str(r['mdate']) + CRT)

    fp.close()

    fp = open(OUTPUTDIR + 'mgi_relationship_property.bcp', 'w')

    results = db.sql('''
    	select p.*,
	       to_char(p.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
	       to_char(p.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	from mgi_relationship_property p, mgi_relationship a
	where a._category_key in (1003,1004)
	and a._relationship_key = p._relationship_key
	''', 'auto')

    for r in results:
	fp.write(`r['_relationshipproperty_key']` + TAB + \
	         `r['_relationship_key']` + TAB + \
	         `r['_propertyname_key']` + TAB + \
	         r['value'] + TAB + \
	         `r['sequenceNum']` + TAB + \
		 str(r['cdate']) + TAB + \
		 str(r['mdate']) + CRT)

    fp.close()

    #
    # mgi_relationship_terms.bcp
    #

    fp = open(OUTPUTDIR + 'mgi_relationship_terms.bcp', 'w')

    results = db.sql('''
	    select v.name, t._Vocab_key, t._Term_key, t.term, 
                to_char(t.creation_date, 'Mon DD YYYY HH:MIAM') as cdate,
                to_char(t.modification_date, 'Mon DD YYYY HH:MIAM') as mdate
	    from VOC_Term t, VOC_Vocab v 
	    where t._Vocab_key in (94, 95, 96, 97)
	    and t._Vocab_key = v._Vocab_key
	    ''', 'auto', execute = 1)
    for r in results:
	    fp.write(`r['_Vocab_key']` + TAB + \
	             `r['_Term_key']` + TAB + \
	             r['name'] + TAB + \
		     r['term'] + TAB + \
		     str(r['cdate']) + TAB + \
		     str(r['mdate']) + CRT)
    fp.close()

#
# Main
#

db.useOneConnection(1)
vocabs()
alleles()
markers()
strains()
genotypes()
references()
omim()
mgi_relationship()
db.useOneConnection(0)

