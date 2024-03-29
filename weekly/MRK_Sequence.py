
'''
#
# MRK_Sequence.py 03/02/99
#
# Report:
#       Tab-delimited file
#
#	1:  MGI Marker Accession ID
#	2:  Marker Symbol
#	3:  Status
#	4:  Marker Type
#	5:  Marker Name
#	6:  cM position
#	7:  Chromosome
#	8:  Genome Coordinate Start
#	9:  Genome Coordinate End
#	10: Strand
#	11: GenBank ID
#	12: RefSeq transcript ID 
#	13: Ensembl transcript ID
#	14: UniProt ID
#	15: TrEMBL ID
#	16: Ensembl protein ID
#	17: RefSeq protein ID
#	18: UniGene ID
#	19: MCV Term
#
# Notes:
#
# History:
#
# sc	09/16/2015
#	-TR11957 - remove constraint of requiring Marker assoc to both
#		genbank and refseq sequences to be included in this 
#		report
#		- add feature type to column 21
#
# sc	06/07/2013
#	-TR11402 - add UniGene ID column (was removed in tag 5-0-0-18 
#		in May of 2012, but don't know why
#
# lec	05/14/2012
#	- TR11035/postgres cleanup/merge
#
# sc	03/12/2010
#	- TR9774 Add Ensembl transcripts
#
# lec	10/10/2006
#	- only include Markers that have at least one sequence.
#
# lec	09/19/2006
#	- make sure there are no deleted sequences in the report
#
# lec	01/27/2005
#	- TR 6529
#
# lec	12/18/2003
#	- TR 5440; added Marker Type
#
# lec	05/30/1999
#	- TR 1631; add UniGene Accession numbers
#
# lec	03/02/1999
#	- TR 130; now use direct Marker-Seq ID relationships
#
# lec	01/18/1999
#	- missing Segment records where relationship is null
#
# lec	01/13/98
#	- added comments section
#
'''
 
import sys
import os
import mgi_utils
import reportlib
import db

db.setTrace()

TAB = reportlib.TAB
CRT = reportlib.CRT

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# header
#
fp.write('MGI Marker Accession ID\t')
fp.write('Marker Symbol\t')
fp.write('Status\t')
fp.write('Marker Type\t')
fp.write('Marker Name\t')
fp.write('cM position\t')
fp.write('Chromosome\t')
fp.write('Genome Coordinate Start\t')
fp.write('Genome Coordinate End\t')
fp.write('Strand\t')
fp.write('GenBank IDs\t')
fp.write('RefSeq transcript IDs\t')
fp.write('Ensembl transcript IDs\t')
fp.write('UniProt IDs\t')
fp.write('TrEMBL IDs\t')
fp.write('Ensembl protein IDs\t')
fp.write('RefSeq protein IDs\t')
fp.write('UniGene IDs\t')
fp.write('Feature Type\n')

# deleted sequences

db.sql('select s._Sequence_key into temporary table deleted from SEQ_Sequence s where s._SequenceStatus_key = 316343', None)
db.sql('create index deleted_idx1 on deleted(_Sequence_key)', None)

db.sql('''
    select a.accID, a._LogicalDB_key 
    into temporary table deletedIDs from deleted d, ACC_Accession a 
    where d._Sequence_key = a._Object_key 
    and a._MGIType_key = 19
    ''', None)
db.sql('create index deletedIDS_idx1 on deletedIDs(accID)', None)
db.sql('create index deletedIDS_idx2 on deletedIDs(_LogicalDB_key)', None)

# all official mouse markers that have at least one Sequence ID
#
#	9	GenBank
#	27	RefSeq transcript (not "xp_" or "np_")
#	133	Ensembl transcript
#

db.sql('''
        select m._Marker_key, m.symbol, m.name, m.chromosome, 
        mlc.genomicChromosome,
        m.cmoffset, upper(substring(s.status, 1, 1)) as markerStatus, t.name as markerType        
        into temporary table markers 
        from MRK_Marker m, MRK_Status s, MRK_Types t,
                MRK_Location_Cache mlc
        where m._Organism_key = 1 
        and m._Marker_Status_key = 1 
        and m._Marker_key = mlc._Marker_key 
        and m._Marker_Status_key = s._Marker_Status_key 
        and m._Marker_Type_key = t._Marker_Type_key 
        and exists (select 1 from ACC_Accession a where m._Marker_key = a._Object_key 
        and a._MGIType_key = 2 and a._LogicalDB_key in (9, 27, 133) and a.prefixPart not in ('XP_', 'NP_'))
        ''', None)
db.sql('create index markers_idx1 on markers(_Marker_key)', None)
db.sql('create index markers_idx2 on markers(symbol)', None)

# MGI ids

results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 1 
      and a.prefixPart = 'MGI:' 
      and a.preferred = 1
      ''', 'auto')
mgiID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    mgiID[key] = value

# 
# Feature types
# 
results = db.sql('''
        select distinct _Marker_key, term
        from MRK_MCV_Cache
        where qualifier = 'D' ''', 'auto')
mcvTerms = {}
for r in results:
    key = r['_Marker_key']
    term = r['term']
    if key not in mcvTerms:
        mcvTerms[key] = []
    mcvTerms[key].append(term)
#
# coordinates
#
results = db.sql('''	
    select m._marker_key,
           c.strand, 
           c.startCoordinate::int as startC,
           c.endCoordinate::int as endC
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

# GenBank ids
results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 9 
      and not exists (select 1 from deletedIDs d where a.accID = d.accID and a._LogicalDB_key = d._LogicalDB_key)
      ''', 'auto')
gbID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if key not in gbID:
        gbID[key] = []
    gbID[key].append(value)
# UniGene ids

results = db.sql('''
      select distinct m._Marker_key, a.accID
      from markers m, ACC_Accession a
      where m._Marker_key = a._Object_key
      and a._MGIType_key = 2
      and a._LogicalDB_key = 23
      and not exists (select 1 from deletedIDs d where a.accID = d.accID and a. _LogicalDB_key = d._LogicalDB_key)
      ''', 'auto')
ugID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if key not in ugID:
        ugID[key] = []
    ugID[key].append(value)

# RefSeq transcript ids
results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 27 
      and a.prefixPart not in ('XP_', 'NP_') 
      and not exists (select 1 from deletedIDs d where a.accID = d.accID and a._LogicalDB_key = d._LogicalDB_key)
      ''', 'auto')
rstrans = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if key not in rstrans:
        rstrans[key] = []
    rstrans[key].append(value)

# RefSeq protein ids
results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 27 
      and a.prefixPart in ('XP_', 'NP_') 
      and not exists (select 1 from deletedIDs d where a.accID = d.accID and a._LogicalDB_key = d._LogicalDB_key)
      ''', 'auto')
rsprot = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if key not in rsprot:
        rsprot[key] = []
    rsprot[key].append(value)

# Ensembl transript IDs
results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 133 
      and not exists (select 1 from deletedIDs d where a.accID = d.accID and a._LogicalDB_key = d._LogicalDB_key)
      ''', 'auto')
enstrans = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if key not in enstrans:
        enstrans[key] = []
    enstrans[key].append(value)

# Ensembl protein IDs
results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 134 
      and not exists (select 1 from deletedIDs d where a.accID = d.accID and a._LogicalDB_key = d._LogicalDB_key)
      ''', 'auto')
ensprot = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if key not in ensprot:
        ensprot[key] = []
    ensprot[key].append(value)

# UniProt ids
results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 13
      and not exists (select 1 from deletedIDs d where a.accID = d.accID and a._LogicalDB_key = d._LogicalDB_key)
      ''', 'auto')
uniprotID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if key not in uniprotID:
        uniprotID[key] = []
    uniprotID[key].append(value)

# TrEMBL ids
results = db.sql('''
      select distinct m._Marker_key, a.accID 
      from markers m, ACC_Accession a 
      where m._Marker_key = a._Object_key 
      and a._MGIType_key = 2 
      and a._LogicalDB_key = 41
      and not exists (select 1 from deletedIDs d where a.accID = d.accID and a._LogicalDB_key = d._LogicalDB_key)
      ''', 'auto')
tremblID = {}
for r in results:
    key = r['_Marker_key']
    value = r['accID']
    if key not in tremblID:
        tremblID[key] = []
    tremblID[key].append(value)

# process

results = db.sql('select * from markers order by symbol', 'auto')

for r in results:

        key = r['_Marker_key']
        symbol = r['symbol']

        #
        # skipping markers that do not cotain a genbank and refseq transcript id
        # sc - 9/16/15, TR11957, remove this constraint
        #if not gbID.has_key(key) and not rstrans.has_key(key):
        #    #print 'not gb', symbol
        #    continue

        if r['cmoffset'] == -1.0:
                cmoffset = 'syntenic'
        elif r['cmoffset'] == -999.0:
                cmoffset = 'N/A'
        else:
                cmoffset = str(r['cmoffset'])

#	1:  MGI Marker Accession ID
#	2:  Marker Symbol
#	3:  Status
#	4:  Marker Type
#	5:  Marker Name
#	6:  cM position
#	7:  Chromosome

        if r['genomicChromosome']:
                chromosome = r['genomicChromosome']
        else:
                chromosome = r['chromosome']

        fp.write(mgiID[key] + TAB + \
                 r['symbol'] + TAB + \
                 r['markerStatus'] + TAB + \
                 r['markerType'] + TAB + \
                 r['name'] + TAB + \
                 cmoffset + TAB + \
                 chromosome + TAB)

        # genome coordinates: column 8-9-10
        if key in coords:
                fp.write(mgi_utils.prvalue(coords[r['_Marker_key']][0]['startC']) + TAB)
                fp.write(mgi_utils.prvalue(coords[r['_Marker_key']][0]['endC']) + TAB)
                fp.write(mgi_utils.prvalue(coords[r['_Marker_key']][0]['strand']) + TAB)
        else:
                fp.write(TAB + TAB + TAB)

#	11: GenBank ID
        if key in gbID:
                fp.write('|'.join(gbID[key]))
        fp.write(TAB)

#	12: RefSeq ID
        if key in rstrans:
                fp.write('|'.join(rstrans[key]))
        fp.write(TAB)

#	13: Ensembl transcript ID
        if key in enstrans:
                fp.write('|'.join(enstrans[key]))
        fp.write(TAB)

#	14: UniProt ID
        if key in uniprotID:
                fp.write('|'.join(uniprotID[key]))
        fp.write(TAB)

#	15: TrEMBL ID
        if key in tremblID:
                fp.write('|'.join(tremblID[key]))
        fp.write(TAB)

#	16: Ensembl protein ID
        if key in ensprot:
                fp.write('|'.join(ensprot[key]))
        fp.write(TAB)

#	17: RefSeq protein ID
        if key in rsprot:
                fp.write('|'.join(rsprot[key]))
        fp.write(TAB)
#	18: UniGene ID
        if key in ugID:
                fp.write('|'.join(ugID[key]))
        fp.write(TAB)

#	19: MCV Term
        if key in mcvTerms:
                fp.write('|'.join(mcvTerms[key]))
        fp.write(CRT)

reportlib.finish_nonps(fp)
