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
#	Species A
#	Species B
#	Chr A
#	Chr B
#	Number of homologies
#
#	where Species A is mouse, human, rat
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

def processSpecies(speciesKey):

    cmds = []

    # select all markers of given species which have a homology

    cmds.append('select distinct m._Marker_key, m.chromosome, m._Species_key, h._Class_key, s.name, c.sequenceNum ' + \
        'into #markers ' + \
	'from MRK_Marker m, HMD_Homology_Marker hm, HMD_Homology h, MRK_Species s, MRK_Chromosome c ' + \
	'where m._Species_key = %s' % (speciesKey) + \
	'and m._Marker_key = hm._Marker_key ' + \
	'and hm._Homology_key = h._Homology_key ' + \
	'and m._Species_key = s._Species_key ' + \
	'and m._Species_key = c._Species_key ' + \
	'and m.chromosome = c.chromosome')

    cmds.append('create nonclustered index idx_key on #markers(_Class_key)')

    # select all distinct homology occurences for primary species
    # verify numbers by checking Oxford Grid
    # sort by secondary species, chromosome of primary species, chromosome of secondary species

    cmds.append('select distinct m._Class_key, ' + \
	'speciesA = m.name, chrA = m.chromosome, speciesAKey = m._Species_key, ' + \
	'speciesB = s2.name, chrB = m2.chromosome, speciesBKey = m2._Species_key  ' + \
	'from #markers m, ' + \
	'HMD_Homology h, HMD_Homology_Marker hm, MRK_Marker m2, MRK_Species s2, MRK_Chromosome c2 ' + \
	'where m._Class_key = h._Class_key ' + \
	'and h._Homology_key = hm._Homology_key ' + \
	'and hm._Marker_key = m2._Marker_key ' + \
	'and m2._Species_key != %s ' % (speciesKey) + \
	'and m2._Species_key = s2._Species_key ' + \
	'and m2._Species_key = c2._Species_key ' + \
	'and m2.chromosome = c2.chromosome ' + \
	'order by m2._Species_key, m.sequenceNum, c2.sequenceNum')

    results = db.sql(cmds, 'auto')

    count = 0
    prevKey = ''

    for r in results[-1]:

	key = str(r['speciesAKey']) + ':' + r['chrA'] + ':' + r['chrB']

	if prevKey != key:

	    if prevKey != '':
		fp.write(str(count) + CRT)

            fp.write(r['speciesA'] + TAB + \
	        r['speciesB'] + TAB + \
	        r['chrA'] + TAB + \
	        r['chrB'] + TAB)

	    prevKey = key
	    count = 0
 
	count = count + 1

    fp.write(str(count) + CRT)

    db.sql('drop table #markers', None)

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)
db.useOneConnection(1)

for speciesKey in ['1', '2', '40']:
    processSpecies(speciesKey)

db.useOneConnection(0)
reportlib.finish_nonps(fp)

