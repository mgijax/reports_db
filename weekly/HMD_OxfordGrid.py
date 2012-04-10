#!/usr/local/bin/python

'''
#
# HMD_OxfordGrid.py 05/14/2003
#
# Report:
#	TR 4805
#	Requested by Donna Maglott (maglott@ncbi.nlm.nih.gov)
#
#	Produce tab-delimited report of:
#
#	Organism A
#	Organism B
#	Chr A
#	Chr B
#	Number of homologies
#
#	where Organism A is mouse, human, rat
#
# Usage:
#       HMD_OxfordGrid.py
#
# History:
#
# lec   05/14/2003
#       - created
#
'''

import sys
import os
import string
import mgi_utils
import reportlib

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db


CRT = reportlib.CRT
TAB = reportlib.TAB

def processOrganism(organismKey):

    # select all markers of given organism which have a homology

    db.sql('''
        select distinct m._Marker_key, m.chromosome, m._Organism_key, hm._Class_key, s.commonName, c.sequenceNum
        into #markers
        from MRK_Marker m, MRK_Homology_Cache hm, MGI_Organism s, MRK_Chromosome c
        where m._Organism_key = %s
        and m._Marker_key = hm._Marker_key
        and m._Organism_key = s._Organism_key
        and m._Organism_key = c._Organism_key
        and m.chromosome = c.chromosome
        ''' % (organismKey), None)

    db.sql('create index idx_key on #markers(_Class_key)', None)

    # select all distinct homology occurences for primary organism
    # verify numbers by checking Oxford Grid
    # sort by secondary organism, chromosome of primary organism, chromosome of secondary organism

    results = db.sql('''
        select distinct m._Class_key, 
        m.commonName as organismA, m.chromosome as chrA, m._Organism_key as organismAKey, 
        s2.commonName as organismB, m2.chromosome as chrB, m2._Organism_key as organismBKey,
        m.sequenceNum as msequenceNum,
        c2.sequenceNum as c2sequenceNum
        from #markers m, 
        MRK_Homology_Cache hm, MRK_Marker m2, MGI_Organism s2, MRK_Chromosome c2 
        where m._Class_key = hm._Class_key 
        and hm._Marker_key = m2._Marker_key 
        and m2._Organism_key != %s 
        and m2._Organism_key = s2._Organism_key 
        and m2._Organism_key = c2._Organism_key 
        and m2.chromosome = c2.chromosome 
        order by m2._Organism_key, msequenceNum, c2sequenceNum
        ''' % (organismKey), 'auto')

    # drop the temp table so the next call to this method can recreate it
    # for the next organism.

    db.sql('drop table #markers', None)

    count = 0
    prevKey = ''

    for r in results:

	key = str(r['organismAKey']) + ':' + r['chrA'] + ':' + r['chrB']

	if prevKey != key:

	    if prevKey != '':
		fp.write(str(count) + CRT)

            fp.write(r['organismA'] + TAB + \
	        r['organismB'] + TAB + \
	        r['chrA'] + TAB + \
	        r['chrB'] + TAB)

	    prevKey = key
	    count = 0
 
	count = count + 1

    fp.write(str(count) + CRT)

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

for organismKey in ['1', '2', '40']:
    processOrganism(organismKey)

reportlib.finish_nonps(fp)

