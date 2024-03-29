
'''
#
# MRK_Reference.py 01/23/2003
#
# Report:
#       TR 4454
#	Tab-delimited file of:
#		Mouse Marker MGI ID
#		Symbol
#		Name
#		Synonyms
#		PubMed IDs for References 
#	Sorted by Symbol
#
# Usage:
#       MRK_Reference.py
#
# History:
#
# lec	01/23/2003
#	- created
#
'''
 
import sys 
import os
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
# mouse markers
#
db.sql('''
        select a.accID, m._Marker_key, m.symbol, m.name 
        into temporary table markers 
        from MRK_Marker m, ACC_Accession a 
        where m._Organism_key = 1 
        and m._Marker_Status_key = 1 
        and m._Marker_key = a._Object_key 
        and a._MGIType_key = 2 
        and a.prefixPart = 'MGI:' 
        and a.preferred = 1
        and a._LogicalDB_key = 1
        ''', None)
db.sql('create index idx_marker on markers(_Marker_key)', None)

#
# refs
#
db.sql('''
        select distinct m._Marker_key, r._Refs_key 
        into temporary table refs 
        from markers m, MRK_Reference r 
        where m._Marker_key = r._Marker_key
        ''', None)
db.sql('create index idx_refs on refs(_Refs_key)', None)

#
# pub med ids
#
results = db.sql('''
        select r._Marker_key, a.accID 
        from refs r, ACC_Accession a 
        where a._MGIType_key = 1 
        and r._Refs_key = a._Object_key 
        and a._LogicalDB_key = 29
        ''', 'auto')
pubmed = {}
for r in results:
        key = r['_Marker_key']
        if key not in pubmed:
                pubmed[key] = []
        pubmed[key].append(r['accID'])

#
# synonyms
#
results = db.sql('''
        select m._Marker_key, s.synonym 
        from markers m, MGI_Synonym s, MGI_SynonymType st 
        where m._Marker_key = s._Object_key 
        and s._MGIType_key = 2 
        and s._SynonymType_key = st._SynonymType_key 
        and st.synonymType = 'exact'
        ''', 'auto')
syn = {}
for r in results:
        key = r['_Marker_key']
        if key not in syn:
                syn[key] = []
        syn[key].append(r['synonym'])

#
# final results
#
results = db.sql('select * from markers order by symbol', 'auto')
for r in results:

        # The list should include only publications with PubMed identifiers. (per TR)

        if r['_Marker_key'] in pubmed:

                fp.write(r['accID'] + TAB + \
                        r['symbol'] + TAB + \
                        r['name'] + TAB)

                if r['_Marker_key'] in syn:
                        fp.write('|'.join(syn[r['_Marker_key']]))

                fp.write(TAB)
                fp.write('|'.join(pubmed[r['_Marker_key']]) + CRT)

reportlib.finish_nonps(fp)	# non-postscript file
