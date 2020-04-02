
'''
#
# MGI_QTLAllele.py
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
#		PubMed ID of Original Reference
#		MGI ID of Gene
#		Gene Symbol
#		RefSeq ID of Gene
#		Ensembl ID
#		Marker Chromosome
#		Marker Start Coordinate
#		Marker End Coordinate
#		Genome Build
#		MP IDs of MP annotations
#
# History:
#
# lec	01/06/2010
#	- TR 9683
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
# include QTLs only

db.sql('''
        select a._Allele_key, a._Marker_key, a.symbol, a.name, t2.term as alleleType, m.symbol as marker
        into temporary table alleles 
        from ALL_Allele a, VOC_Term t1, VOC_Term t2, MRK_Marker m 
        where a._Allele_Status_key = t1._Term_key 
        and t1.term = 'Approved' 
        and a.isWildType = 0 
        and a._Allele_Type_key = t2._Term_key 
        and a._Marker_key = m._Marker_key 
        and m._Marker_Type_key = 6
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
        
# Retrieve marker coordinates

results = db.sql ('''select _Marker_key, chromosome, genomicChromosome,
                startCoordinate::int as startCoordinate,
                endCoordinate::int as endCoordinate,
                strand, version
        from MRK_Location_Cache
        where _Organism_key = 1
                and startCoordinate is not null''', 'auto')

coords = {}
for r in results:
        coords[r['_Marker_key']] = r

# Retrieve MP IDs for MP annotations

results = db.sql('''
        select distinct s._Allele_key, a.accID 
        from alleles s, GXD_AlleleGenotype ga, VOC_AnnotHeader na, ACC_Accession a 
        where s._Allele_key = ga._Allele_key 
        and ga._Genotype_key = na._Object_key 
        and na._AnnotType_key = 1002 
        and na._Term_key = a._Object_key 
        and a._MGIType_key = 13 
        and a.preferred = 1
        ''', 'auto')

mpIDs = {}
for r in results:
        if r['_Allele_key'] not in mpIDs:
                mpIDs[r['_Allele_key']] = []
        mpIDs[r['_Allele_key']].append(r['accID'])
        
results = db.sql('select * from alleles order by marker, symbol', 'auto')

for r in results:
        fp.write(amgiIDs[r['_Allele_key']] + reportlib.TAB + \
                 r['symbol'] + reportlib.TAB + \
                 r['name'] + reportlib.TAB + \
                 r['alleleType'] + reportlib.TAB)

        if r['_Allele_key'] in pubIDs:
                fp.write(pubIDs[r['_Allele_key']])
        fp.write(reportlib.TAB)

        # if Transgene, do not print gene ID or gene symbol
        if str.find(r['symbol'], 'Tg(') < 0:
                fp.write(mmgiIDs[r['_Marker_key']] + reportlib.TAB + \
                        r['marker'] + reportlib.TAB)
        else:
                fp.write(reportlib.TAB + reportlib.TAB)

        if r['_Marker_key'] in refIDs:
                fp.write(refIDs[r['_Marker_key']])
        fp.write(reportlib.TAB)

        if r['_Marker_key'] in ensemblIDs:
                fp.write(ensemblIDs[r['_Marker_key']])
        fp.write(reportlib.TAB)

        if r['_Marker_key'] in coords:
                coord = coords[r['_Marker_key']]

                # prefer a genomic chromosome, but fall back on a genetic one
                # if the marker has no coordinates
                chrom = coord['genomicChromosome']
                if not chrom:
                        chrom = coord['chromosome']

                start = coord['startCoordinate']
                end = coord['endCoordinate']
                strand = coord['strand']
                build = coord['version']	# not displayed for now

                for field in (chrom, start, end, build):
                        if field == None:
                                fp.write ('null' + reportlib.TAB)
                        else:
                                fp.write (str(field) + reportlib.TAB) 
        else:
                for i in range(0,4):
                        fp.write ('null' + reportlib.TAB)


        if r['_Allele_key'] in mpIDs:
                fp.write(','.join(mpIDs[r['_Allele_key']]))

        fp.write(reportlib.CRT)

reportlib.finish_nonps(fp)
