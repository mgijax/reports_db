
'''
# Name: MGI_DO.py
#
# Lists mouse markers and human markers with disease annotations.
# 
# Report columns:
#
# 1) DO Disease ID
# 2) DO Disease Name
# 3) OMIM IDs (piped-delimited)
# TR13349 removed: 4) HomoloGene ID (for the homology class of the marker, if one)
# 5) Common Organism Name
# 6) NCBI Taxon ID (for the organism)
# 7) Symbol (for the marker)
# 8) EntrezGene ID
# 9) Mouse MGI ID (mouse markers only)
#
# Sort by: 2, 3, 5 (above)
#
# History:
#
# sc    02/22/2021
#       TR13349 - B39 project. Update to use alliance direct homology
#               (was using Homologene)
#
# 06/05/2017	lec
#	- remove complicated sort/use simple sort
#
# 11/03/2016	lec
#	- TR12427/Disease Ontology (DO)
#
'''

import sys
import os
import db
import reportlib
import mgi_utils

#db.setTrace()

CRT = reportlib.CRT
TAB = reportlib.TAB

#
# get the disease associations for mouse markers (1023) and human (1022) markers.
#

cmd = '''(
        select distinct va._Object_key, va._Term_key, t.term, m.symbol, o.commonName, tx.accID as taxID, a.accID
        into temp table diseaseannot
        from VOC_Annot va, VOC_Term t, MRK_Marker m, MGI_Organism o, ACC_Accession tx, ACC_Accession a
        where va._AnnotType_key = 1023
        and va._Qualifier_key != 1614157
        and va._Term_key = t._Term_key
        and va._Object_key = m._Marker_key
        and m._Organism_key = o._Organism_key
        and o._Organism_key = tx._Object_key
        and tx._MGIType_key = 20
        and tx._LogicalDB_key = 32
        and va._Object_key = a._Object_key
        and a._LogicalDB_key = 1
        and a._MGIType_key = 2
        and a.preferred = 1
        union
        select distinct va._Object_key, va._Term_key, t.term, m.symbol, o.commonName, tx.accID as taxID, null
        from VOC_Annot va, VOC_Term t, MRK_Marker m, MGI_Organism o, ACC_Accession tx
        where va._AnnotType_key = 1022
        and va._Term_key = t._Term_key
        and va._Object_key = m._Marker_key
        and m._Organism_key = o._Organism_key
        and o._Organism_key = tx._Object_key
        and tx._MGIType_key = 20
        and tx._LogicalDB_key = 32
        )
        '''

db.sql(cmd, None)
db.sql('create index idx_obj on diseaseannot(_Object_key)', None)
db.sql('create index idx_term on diseaseannot(_Term_key)', None)

#
# cache DO disease terms and IDs, and OMIM IDs
#
results = db.sql('''select t._Term_key, t.term, a.accID, omim.accid as omimID
        from VOC_Term t, ACC_Accession a, ACC_Accession omim
        where t._Vocab_key = 125
        and t._Term_key = a._Object_key
        and a._MGIType_key = 13
        and a.preferred = 1
        and a._LogicalDB_key = 191
        and a.preferred = 1
        and a._Object_key = omim._Object_key
        and omim._LogicalDB_key  = 15
        union
        select t._Term_key, t.term, a.accID, null
        from VOC_Term t, ACC_Accession a
        where t._Vocab_key = 125
        and t._Term_key = a._Object_key
        and a._MGIType_key = 13
        and a.preferred = 1
        and a._LogicalDB_key = 191
        and a.preferred = 1
        and not exists (select 1 from ACC_Accession omim
                where a._Object_key = omim._Object_key
                and omim._LogicalDB_key  = 15
                )
        ''', 'auto')

diseaseLookup = {}
omimLookup = {}
for r in results:
    key = r['_Term_key']

    value = (r['term'], r['accID'])
    if key not in diseaseLookup:
        diseaseLookup[key] = value

    value = r['omimID']
    if value != None:
        if key not in omimLookup:
            omimLookup[key] = []
        omimLookup[key].append(value)

#
# cache the EntrezGene ID by Marker
#
entrezGeneLookup = {}
results = db.sql('''select m._Object_key, a.accID
        from diseaseannot m, ACC_Accession a
        where m._Object_key = a._Object_key
        and a._MGIType_key = 2
        and a.preferred = 1
        and a._LogicalDB_key = 55
        and a.private = 0''', 'auto')
for r in results:
    key = r['_Object_key']
    value = r['accID']
    entrezGeneLookup[key] = []
    entrezGeneLookup[key].append(value)

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)
fp.write('DO Disease ID' + TAB)
fp.write('DO Disease Name' + TAB)
fp.write('OMIM IDs' + TAB)
fp.write('Common Organism Name' + TAB)
fp.write('NCBI Taxon ID' + TAB)
fp.write('Symbol' + TAB)
fp.write('EntrezGene ID' + TAB)
fp.write('Mouse MGI ID' + CRT)

results = db.sql('select * from diseaseannot order by term, commonName', 'auto')
for r in results:

    fp.write(diseaseLookup[r['_Term_key']][1] + TAB)
    fp.write(diseaseLookup[r['_Term_key']][0] + TAB)

    if r['_Term_key'] in omimLookup:
        fp.write('|'.join(omimLookup[r['_Term_key']]))
    fp.write(TAB)

    fp.write(r['commonName'] + TAB)
    fp.write(r['taxID'] + TAB)
    fp.write(r['symbol'] + TAB)

    if r['_Object_key'] in entrezGeneLookup:
        fp.write('|'.join(entrezGeneLookup[r['_Object_key']]))
    fp.write(TAB)

    fp.write(mgi_utils.prvalue(r['accID']) + CRT)

reportlib.finish_nonps(fp)
