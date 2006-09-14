#!/usr/local/bin/python

'''
#
# goa.py
#
# Report:
#	TR 7904
#
#	Takes the GOA file ${GOAINPUTFILE} and generates
#
#		1. goa.error
#			file of UniProtIDS that are not found in MGI
#	
#		2. goa.duplicates
#
#		   	file of GOA annotations that are duplicates to those in MGI
#
#		3. goa.forMGI
#
#			file of GOA annotations that are to be appended to the
#			GO_gene_association.py output file
#
# The GOA format has the following columns:
#
#   1.  Database designation (MGI)
#   2.  MGI Marker ID
#   3.  Symbol
#   4.  NOT
#   5.  GO id
#   6.  MGI ID of Reference (in format MGI:MGI:####) (double MGI: necessary)
#   7.  Evidence abbreviation
#   8.  Inferred From
#   9.  GO DAG Abbreviation (F, P, C)
#   10. Gene name
#   11. Gene synonym(s) - list of |-delimited synonyms
#   12. Marker Type (gene)
#   13. Species (taxon:10090)
#   14. Modification Date (YYYYMMDD)
#   15. Assigned By
#
# Usage:
#       goa.py
#
# History:
#
# lec   09/14/2005
#       - created
#
'''

import sys
import os
import string
import db
import mgi_utils
import reportlib

inFileName = os.environ['GOAINPUTFILE']

errorFile = reportlib.init('goa', outputdir = os.environ['GOADIR'], printHeading = 0, fileExt = '.error')
dupFile = reportlib.init('goa', outputdir = os.environ['GOADIR'], printHeading = 0, fileExt = '.duplicates')
mgiFile = reportlib.init('goa', outputdir = os.environ['GOADIR'], printHeading = 0, fileExt = '.mgi')

assoc = {}	# dictionary of ID:Marker key associations
marker = {}	# dictionary of Marker key:MGI ID
annot = []	# list of existing Marker key, GO ID, Evidence Code, Pub Med ID annotations
annotByGOID = []
annotByRef = []

mgiLine = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'

#
# SwissProt (13)
# TrEMBL (41)
# RefSeq (27)
# ENSEMBL (60)
# VEGA (85)
#

db.sql('select a._Object_key, a.accID into #assoc ' + \
	'from ACC_Accession a ' + \
	'where a._LogicalDB_key in (13, 41, 27, 60, 85) ' + \
	'and a._MGIType_key = 2', None)
db.sql('create index idx1 on #assoc(_Object_key)', None)

results = db.sql('select _Object_key, accID from #assoc', 'auto')
for r in results:
    key = r['accID']
    value = r['_Object_key']

    if not assoc.has_key(key):
	assoc[key] = []

    assoc[key].append(value)

results = db.sql('select distinct a._Object_key, a.accID, m.symbol, m.name, markerType = t.name ' + \
	'from #assoc u, ACC_Accession a, MRK_Marker m, MRK_Types t ' + \
	'where a._MGIType_key = 2 ' + \
	'and u._Object_key = a._Object_key ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1 ' + \
	'and a._Object_key = m._Marker_key ' + \
	'and m._Marker_Type_key = t._Marker_Type_key', 'auto')
for r in results:
    marker[r['_Object_key']] = r

#
# existing GO annotations that have pub med ids
# to detect duplicate annotations
#

results = db.sql('select goID = a.accID, t._Object_key, ec.term, refID = r.accID ' + \
	'from VOC_Annot t, ACC_Accession a, VOC_Evidence e, VOC_Term ec, ACC_Accession r ' + \
	'where t._AnnotType_key = 1000 ' + \
	'and t._Term_key = a._Object_key ' + \
	'and a._MGIType_key = 13 ' + \
	'and a.preferred = 1 ' + \
	'and t._Annot_key = e._Annot_key ' + \
	'and e._EvidenceTerm_key = ec._Term_key ' + \
	'and e._Refs_key = r._Object_key ' + \
	'and r._MGIType_key = 1 ' + \
	'and r._LogicalDB_key = 29', 'auto')
for r in results:
    annot.append((r['goID'], r['_Object_key'], r['term'], r['refID']))

    if r['goID'] not in annotByGOID:
        annotByGOID.append(r['goID'])

    if r['refID'] not in annotByRef:
        annotByRef.append(r['refID'])

inFile = open(inFileName, 'r')

for line in inFile.readlines():

    tokens = string.split(line[:-1], '\t')
    databaseID = tokens[0]
    goaID = tokens[1]	# translate
    goaSymbol = tokens[2]	# translate
    notValue = tokens[3]
    goID = tokens[4]
    refID = tokens[5]
    evidence = tokens[6]
    inferredFrom = tokens[7]
    dag = tokens[8]
    goaName = tokens[9]	# translate
    synonyms = tokens[10]
    markerType = tokens[11]	# translate
    taxID = tokens[12]
    modDate = tokens[13]
    assignedBy = tokens[14]

    # skip it if it's assigned by MGI; that means it came from us to begin with

    if assignedBy == 'MGI':
	continue

    # error if GOA id is not found

    if not assoc.has_key(goaID):
	errorFile.write(line + '\n')
	continue

    # error if GOA id maps to more than one MGI Marker

    if len(assoc[goaID]) > 1:
	errorFile.write(line + '\n')
	continue

    objectKey = assoc[goaID][0]
    m = marker[objectKey]

    # duplicate error if the annotation already exists in MGI

    if refID in annotByRef and goID in annotByGOID:
        goaAnnot = (goID, objectKey, evidence, refID)
        if goaAnnot in annot:
	    dupFile.write(line + '\n')
	    continue

    mgiFile.write(mgiLine % (databaseID, m['accID'], m['symbol'], notValue, goID, refID, evidence, inferredFrom,\
	dag, m['name'], synonyms, m['markerType'], taxID, modDate, assignedBy))

inFile.close()

reportlib.finish_nonps(errorFile)
reportlib.finish_nonps(dupFile)
reportlib.finish_nonps(mgiFile)
