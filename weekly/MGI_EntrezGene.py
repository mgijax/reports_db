
'''
#
# Report:
#       Tab-delimited file of MGI Mouse Markers including Withdrawns
#	Prints current MGI accession ID of Marker record
#	and list of non-preferred MGI accession IDs of each Marker
#
# Usage:
#       MGI_EntrezGene.py
#
# Used by:
# 	TR 1283 - Donna Maglott at NCBI for EntrezGene
#
# Notes:
#
# Splits will retain their MGI Accession IDs.  The 3rd query is to find
# withdrawns which still have a preferred MGI Accession ID.
#
# History:
#
#
# sc    06/27/2014
#       - TR11560 Feature Relationships project
#       - exclude feature type 'mutation defined region'
#
# lec   04/23/2012
#       - TR11035/postgres cleanup; added feature type, coordinates & biotype
#
# lec	10/19/2011
#	- TR10885/add column 11/raw biotypes, column 12/feature type per Donna Maglott/NCBI
#
# lec	06/18/2002
#	- rewrote to use dictionaries for egID, other Acc IDs, other names
#
# lec	01/19/2000
#	- created
#
'''
 
import sys
import os
import mgi_utils
import reportlib
import db

db.setTrace()

CRT = reportlib.CRT
TAB = reportlib.TAB

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# 1. all Approved Marker records
# 2. all Withdrawn Marker records
#

db.sql('''
  select m._Marker_key, m.symbol, m.name, m._Marker_key as _Current_key, 
        m.cmoffset, m.chromosome, t.name as markerType, 1 as isPrimary, 
        upper(substring(s.status, 1, 1)) as markerStatus
  into temporary table markers 
  from MRK_Marker m, MRK_Types t, MRK_Status s 
  where m._Organism_key = 1 
  and m._Marker_Status_key = 1
  and m._Marker_Type_key = t._Marker_Type_key 
  and m._Marker_Status_key = s._Marker_Status_key 
  union 
  select m._Marker_key, m.symbol, m.name, c._Current_key, 
  m.cmoffset, m.chromosome, t.name as markerType, 0 as isPrimary, 
  upper(substring(s.status, 1, 1)) as markerStatus
  from MRK_Marker m, MRK_Current c, MRK_Types t, MRK_Status s 
  where m._Organism_key = 1 
  and m._Marker_Status_key = 2 
  and m._Marker_key = c._Marker_key 
  and m._Marker_Type_key = t._Marker_Type_key 
  and m._Marker_Status_key = s._Marker_Status_key 
  ''', None)
db.sql('create index idx_key on markers(_Marker_key)', None)

# MGI ids

results = db.sql('''
        select m._Current_key, a.accID 
        from markers m, ACC_Accession a 
        where m._Current_key = a._Object_key 
        and a._MGIType_key = 2 
        and a.prefixPart = 'MGI:' 
        and a._LogicalDB_key = 1 
        and a.preferred = 1 
        ''', 'auto')
mgiID = {}
for r in results:
    key = r['_Current_key']
    value = r['accID']
    mgiID[key] = value

# Get EntrezGene ID for Primary Markers
results = db.sql('''
        select m._Marker_key, a.accID 
        from markers m, ACC_Accession a 
        where m.isPrimary = 1 
        and m._Marker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a._LogicalDB_key = 55 
        ''', 'auto')
egID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    egID[key] = value

# Get Secondary MGI Ids for Primary Marker
results = db.sql('''
        select m._Marker_key, a.accID 
        from markers m, ACC_Accession a 
        where m.isPrimary = 1 
        and m._Marker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a.prefixPart = 'MGI:' 
        and a._LogicalDB_key = 1 
        and a.preferred = 0
        ''', 'auto')
otherAccId = {}
for r in results:
        if r['_Marker_key'] not in otherAccId:
                otherAccId[r['_Marker_key']] = []
        otherAccId[r['_Marker_key']].append(r['accID'])

# Get Synonyms for Primary Marker
results = db.sql('''
        select m._Marker_key, s.synonym 
        from markers m, MGI_Synonym s, MGI_SynonymType st 
        where m.isPrimary = 1 
        and m._Marker_key = s._Object_key 
        and s._MGIType_key = 2 
        and s._SynonymType_key = st._SynonymType_key 
        and st.synonymType = 'exact'
        ''', 'auto')
synonym = {}
for r in results:
        key = r['_Marker_key']
        value = r['synonym']
        if key not in synonym:
                synonym[key] = []
        synonym[key].append(value)

# Get BioType for Primary Marker
results = db.sql('''
        select distinct m._Marker_key, s.rawbiotype 
        from markers m, SEQ_Marker_Cache s 
        where m.isPrimary = 1 
        and m._Marker_key = s._Marker_key 
        and s.rawbiotype is not null
        ''', 'auto')
bioTypes = {}
for r in results:
        key = r['_Marker_key']
        value = r['rawbiotype']
        if key not in bioTypes:
                bioTypes[key] = []
        bioTypes[key].append(value)

# Get Feature Type for Primary Marker
results = db.sql('''
        select m._Marker_key, s.term, s._MCVTerm_key
        from markers m, MRK_MCV_Cache s 
        where m.isPrimary = 1 
        and m._Marker_key = s._Marker_key 
        and s.qualifier = 'D'
        ''', 'auto')

# {markerKey:[list of feature types]
featureTypes = {}
# {markerKey:[list of feature type keys]
featureTypeByKey = {}

for r in results:
        markerKey = r['_Marker_key']
        value = r['term']
        termKey = r['_MCVTerm_key']
        if markerKey not in featureTypes:
                featureTypes[markerKey] = []
        featureTypes[markerKey].append(value)
        if markerKey not in featureTypeByKey:
                featureTypeByKey[markerKey] = []
        featureTypeByKey[markerKey].append(termKey)


#
# coordinates
#
results = db.sql('''	
    select m._marker_key,
           c.strand, 
           c.startCoordinate::int as startC,
           c.endCoordinate::int as endC,
           c.genomicChromosome
    from markers m, MRK_Location_Cache c
    where m._marker_key = c._marker_key
        ''', 'auto')
coords = {}
for r in results:
    key = r['_marker_key']
    value = r
    if key not in coords:
        coords[key] = []
    coords[key].append(value)

#
# final query
#
results = db.sql('select * from markers order by _Current_key, isPrimary desc', 'auto')
for r in results:

        # column 1: mgi id
        # column 2: symbol
        # column 3: status
        # column 4: name
        # column 5: cmoffset
        # column 6: chromosome
        # column 7: marker type

        if r['_Current_key'] not in mgiID:
           continue

        key = r['_Marker_key']

        # if the marker's feature type is not
        # 'mutation defined region', key=11928467 write out to the report
        # default feature type
        fTypes = ''
        if key in featureTypes:
            mcvKeyList = featureTypeByKey[key]
            if 11928467 in mcvKeyList:
                continue
            else:
                fTypes = ('|'.join(featureTypes[key]))

        chromosome = None

        if r['_Marker_key'] in coords:
                # genomic chromosome, if marker has coordinates
                chromosome = coords[r['_Marker_key']][0]['genomicChromosome']

        if not chromosome:
                # genetic chromosome, if marker has no coordinates
                chromosome = r['chromosome']

        fp.write(mgiID[r['_Current_key']] + TAB + \
                r['symbol'] + TAB + \
                r['markerStatus'] + TAB + \
                r['name'] + TAB + \
                str(r['cmoffset']) + TAB + \
                chromosome + TAB + \
                r['markerType'] + TAB)

        if r['isPrimary']:

                # column 8: other accession ids
                if r['_Marker_key'] in otherAccId:	
                        fp.write('|'.join(otherAccId[r['_Marker_key']]))
                fp.write(TAB)

                # column 9: entrezgene ids
                if r['_Marker_key'] in egID:	
                        fp.write(egID[r['_Marker_key']])
                fp.write(TAB)

                # column 10: synonyms
                if r['_Marker_key'] in synonym:	
                        fp.write('|'.join(synonym[r['_Marker_key']]))
                fp.write(TAB)

                # column 11: feature types
                fp.write(fTypes + TAB)

                # column 12-13-14
                if r['_Marker_key'] in coords:
                    fp.write(mgi_utils.prvalue(coords[r['_Marker_key']][0]['startC']) + TAB)
                    fp.write(mgi_utils.prvalue(coords[r['_Marker_key']][0]['endC']) + TAB)
                    fp.write(mgi_utils.prvalue(coords[r['_Marker_key']][0]['strand']) + TAB)
                else:
                    fp.write(TAB + TAB + TAB)

                # column 15: biotypes
                if r['_Marker_key'] in bioTypes:	
                        fp.write('|'.join(bioTypes[r['_Marker_key']]))
                fp.write(CRT)

        else:
                # column 8-15 null
                fp.write(TAB)
                fp.write(TAB)
                fp.write(TAB)
                fp.write(TAB)
                fp.write(TAB)
                fp.write(TAB)
                fp.write(TAB)
                fp.write(CRT)

reportlib.finish_nonps(fp)
