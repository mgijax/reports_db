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
# Notes:
#       - all reports use db default of public login
#       - all reports use server/database default of environment
#       - use lowercase for all SQL commands (i.e. select not SELECT)
#       - all public SQL reports require the header and footer
#       - all private SQL reports require the header
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
import db
import mgi_utils
import reportlib

CRT = reportlib.CRT
TAB = reportlib.TAB

def processOrganism(organismKey):

    # select all markers of given organism which have a homology

    db.sql('select distinct m._Marker_key, m.chromosome, m._Organism_key, hm._Class_key, s.commonName, c.sequenceNum ' + \
        'into #markers ' + \
	'from MRK_Marker m, MRK_Homology_Cache hm, MGI_Organism s, MRK_Chromosome c ' + \
	'where m._Organism_key = %s' % (organismKey) + \
	'and m._Marker_key = hm._Marker_key ' + \
	'and m._Organism_key = s._Organism_key ' + \
	'and m._Organism_key = c._Organism_key ' + \
	'and m.chromosome = c.chromosome', None)

    db.sql('create nonclustered index idx_key on #markers(_Class_key)', None)

    # select all distinct homology occurences for primary organism
    # verify numbers by checking Oxford Grid
    # sort by secondary organism, chromosome of primary organism, chromosome of secondary organism

    results = db.sql('select distinct m._Class_key, ' + \
	'organismA = m.commonName, chrA = m.chromosome, organismAKey = m._Organism_key, ' + \
	'organismB = s2.commonName, chrB = m2.chromosome, organismBKey = m2._Organism_key  ' + \
	'from #markers m, ' + \
	'MRK_Homology_Cache hm, MRK_Marker m2, MGI_Organism s2, MRK_Chromosome c2 ' + \
	'where m._Class_key = hm._Class_key ' + \
	'and hm._Marker_key = m2._Marker_key ' + \
	'and m2._Organism_key != %s ' % (organismKey) + \
	'and m2._Organism_key = s2._Organism_key ' + \
	'and m2._Organism_key = c2._Organism_key ' + \
	'and m2.chromosome = c2.chromosome ' + \
	'order by m2._Organism_key, m.sequenceNum, c2.sequenceNum', 'auto')

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

