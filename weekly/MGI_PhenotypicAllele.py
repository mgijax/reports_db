
'''
#
# MGI_PhenotypicAllele.py 05/04/2004
#
# Report:
#       Tab-delimited file
#       Allele Nomenclature
#
#	Fields:
#		MGI ID of Allele
#		Allele Symbol
#		Allele Name
#		Allele Type
#		Allele Attribute/Subtype
#		PubMed ID of Original Reference
#		MGI ID of Gene
#		Gene Symbol
#		RefSeq ID of Gene
#		Ensembl ID
#		MP IDs of MP annotations
#		Synonyms
#		Marker Name
#
# Usage:
#       MGI_PhenotypicAllele.py
#
# History:
#
# sc	08/10/2018
#	- TR12909 - add marker name as last column of report
#
# lec	07/06/2017
#	- TR12593/exclude obsolete MP terms
#
# lnh 07/14/2016
#     - TR12373 -  print gene ID or gene symbol for all allele types including Transgenes
#
# lec	02/03/2014
#	- TR11515/allele type change/add allele-attribute
#
# lec	08/31/2010
#	- TR10280 fix; _AlleleType_key was put in for testing but not removed for release
#
# lec	07/06/2010
#	- TR 10280/add synonyms
#
# lec	05/04/2004
#	- TR 5637
#
'''
 
import sys
import os
import string
import reportlib
import db

db.setTrace()

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('#\n#Allele symbols are usually of the form xxx<yy>, where the <> enclose the part of the symbol that is superscripted.\n#\n')
fp.write('#Transgene insertions, with symbols of the form Tg(aaa)##bbb, are included in this listing, but notably have no corresponding gene marker.\n#\n')
fp.write('#For details of nomenclature rules, see http://www.informatics.jax.org/mgihome/nomen/index.shtml\n#\n')

# Retrieve all Approved Alleles
# exclude wild types
# exclude QTLs

db.sql('''
       select a._Allele_key, a._Marker_key, a.symbol, a.name, t2.term as alleleType, m.symbol as marker, m.name as markerName, a._Allele_Type_key 
       into temporary table alleles 
       from ALL_Allele a, VOC_Term t1, VOC_Term t2, MRK_Marker m 
       where a._Allele_Status_key = t1._Term_key 
       and t1.term in ('Approved')
       and a.isWildType = 0 
       and a._Transmission_key != 3982953 
       and a._Allele_Type_key = t2._Term_key 
       and a._Marker_key = m._Marker_key 
       and m._Marker_Type_key != 6
       ''', None)
db.sql('create index idx1 on alleles(_Allele_key)', None)

# Retrieve MGI Accession number for Allele

results = db.sql('''
        select s._Allele_key, a.accID 
        from alleles s, ACC_Accession a 
        where s._Allele_key = a._Object_key 
        and a._MGIType_key = 11 
        and a._LogicalDB_key = 1 
        and a.prefixPart = 'MGI:'
        and a.preferred = 1
        ''', 'auto')
amgiIDs = {}
for r in results:
        amgiIDs[r['_Allele_key']] = r['accID']
        
# Retrieve MGI Accession number for Marker

results = db.sql('''
        select s._Marker_key, a.accID 
        from alleles s, ACC_Accession a 
        where s._Marker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a._LogicalDB_key = 1 
        and a.prefixPart = 'MGI:'
        and a.preferred = 1
        ''', 'auto')
mmgiIDs = {}
for r in results:
        mmgiIDs[r['_Marker_key']] = r['accID']

# Retrieve PubMed IDs for Original Reference

results = db.sql('''
        select s._Allele_key, b.accID 
        from alleles s, MGI_Reference_Assoc r, MGI_RefAssocType rt, ACC_Accession b 
        where s._Allele_key = r._Object_key 
        and r._MGIType_key = 11 
        and r._RefAssocType_key = rt._RefAssocType_key 
        and rt.assocType = 'Original' 
        and r._Refs_key = b._Object_key 
        and b._MGIType_key = 1 
        and b._LogicalDB_key = 29
        ''', 'auto')
pubIDs = {}
for r in results:
        pubIDs[r['_Allele_key']] = r['accID']

# Retrieve RefSeq ID for Gene

results = db.sql('''
        select s._Marker_key, a.accID 
        from alleles s, ACC_Accession a 
        where s._Marker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a._LogicalDB_key = 27 
        and a.prefixPart in ('NM_', 'XM_')
        ''', 'auto')
refIDs = {}
for r in results:
        refIDs[r['_Marker_key']] = r['accID']
        
# Retrieve Ensembl Gene Model ID for Gene

results = db.sql('''
        select s._Marker_key, a.accID 
        from alleles s, ACC_Accession a 
        where s._Marker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a._LogicalDB_key = 60 
        ''', 'auto')
ensemblIDs = {}
for r in results:
        ensemblIDs[r['_Marker_key']] = r['accID']
        
# Retrieve MP IDs for MP annotations

results = db.sql('''
        select distinct s._Allele_key, a.accID 
        from alleles s, GXD_AlleleGenotype ga, VOC_AnnotHeader na, ACC_Accession a, VOC_Term tt
        where s._Allele_key = ga._Allele_key 
        and ga._Genotype_key = na._Object_key 
        and na._AnnotType_key = 1002 
        and na._Term_key = a._Object_key 
        and a._MGIType_key = 13 
        and a.preferred = 1
        and na._Term_key = tt._Term_key
        and tt.isObsolete = 0
        ''', 'auto')

mpIDs = {}
for r in results:
        if r['_Allele_key'] not in mpIDs:
                mpIDs[r['_Allele_key']] = []
        mpIDs[r['_Allele_key']].append(r['accID'])
        
# Retrieve Synonyms
results = db.sql('''
        select s._Allele_key, ss.synonym 
        from alleles s, MGI_Synonym ss 
        where s._Allele_key = ss._Object_key 
        and ss._MGIType_key = 11 
        ''', 'auto')
synonym = {}
for r in results:
        key = r['_Allele_key']
        value = r['synonym']
        if key not in synonym:
                synonym[key] = []
        synonym[key].append(value)

# Retrieve Allele Attribute/Subtype
results = db.sql('''
        select s._Allele_key, t.term
        from alleles s, VOC_Annot a, VOC_Term t
        where s._Allele_key = a._Object_key 
        and a._AnnotType_key = 1014
        and a._Term_key = t._Term_key
        ''', 'auto')
subtype = {}
for r in results:
        key = r['_Allele_key']
        value = r['term']
        if key not in subtype:
                subtype[key] = []
        subtype[key].append(value)

#
# Main
#

results = db.sql('select * from alleles order by marker, symbol', 'auto')

for r in results:

        fp.write(amgiIDs[r['_Allele_key']] + reportlib.TAB + \
                 r['symbol'] + reportlib.TAB + \
                 r['name'] + reportlib.TAB + \
                 r['alleleType'] + reportlib.TAB)

        # allele attribute/sub-type/pipedelimited
        if r['_Allele_key'] in subtype:
                fp.write('|'.join(subtype[r['_Allele_key']]))
        fp.write(reportlib.TAB)

        if r['_Allele_key'] in pubIDs:
                fp.write(pubIDs[r['_Allele_key']])
        fp.write(reportlib.TAB)

        fp.write(mmgiIDs[r['_Marker_key']] + reportlib.TAB)
        fp.write(r['marker'] + reportlib.TAB)

        if r['_Marker_key'] in refIDs:
                fp.write(refIDs[r['_Marker_key']])
        fp.write(reportlib.TAB)

        if r['_Marker_key'] in ensemblIDs:
                fp.write(ensemblIDs[r['_Marker_key']])
        fp.write(reportlib.TAB)

        if r['_Allele_key'] in mpIDs:
                fp.write(','.join(mpIDs[r['_Allele_key']]))
        fp.write(reportlib.TAB)

        if r['_Allele_key'] in synonym:
                fp.write('|'.join(synonym[r['_Allele_key']]))
        fp.write(reportlib.TAB)
        fp.write(r['markerName'])
        fp.write(reportlib.CRT)

reportlib.finish_nonps(fp)
