#!/usr/local/bin/python

'''
#
# HMD_Rat.py 02/04/2004
#
# Report:
#	TR 5539
#	RGD would like Rat/Human reports similar to Human/Human format.
#
# Usage:
#       HMD_RatHuman.py [1234]
#	1 = sort by Rat Chromosome
#	2 = sort by Human Chromosome
#	3 = sort by Rat Symbol
#	4 = sort by Human Symbol
#	none = generate ALL 4 reports
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
# lec	02/04/2004
#	- TR 5539; created; converted from HMD_Rat.py
#
'''

import sys
import os
import string
import regex
import db
import mgi_utils
import reportlib

CRT = reportlib.CRT
SPACE = reportlib.SPACE
REPORTNAME = 'HMD_RatHuman'

#
# NOTE: The following globals and getSortableOffset are variations on
# the original algorithm which can be found in the WI homology_report.cgi.py
# python module.
#

# special global variables used by getSortableOffset().  We declare them
# outside the function to yield a small speed benefit (we only create them
# once).

offset_cre = regex.compile ('\([pq]\)\([0-9\.]+\)')     # eg- q23.2
offset_special = {
        'pter'  : -10000.0,     # p terminus
        'p'     :     -1.0,     # p arm
        'cen'   :      0.0,     # centromere
        'q'     :      1.0,     # q arm
        'qter'  :  10000.0,     # q terminus
        'None'  :  12000.0,     # no defined offset - go to bottom (syntenic)
        'pter-qter'  :  12000.0 # equivalent to syntenic - go to bottom
        }

def getSortableOffset (cytogeneticOffset):
        # Purpose: use 'cytogeneticOffset' to generate and return a value
        #       which will sort as desired.  (using rat as a primary
        #       species requires special handling, as of TR 211)
        # Returns: a float
        # Assumes: nothing
        # Effects: nothing
        # Throws: nothing
        # Notes: 
        #       the markers on a
        #       chromosome sort by their cytogenetic offset as follows:
        #               1. pter         (p terminus)
        #               2. p<number>    (number in descending order)
        #               3. p            (specifies the whole p arm)
        #               4. cen          (centromere)
        #               5. q            (specifies the whole q arm)
        #               6. q<number>    (number in ascending order)
        #               7. qter         (q terminus)

	cyto_str = str(cytogeneticOffset)

	if offset_special.has_key (cyto_str):   # try to match whole offset
		return offset_special [cyto_str]

	cyto = string.split (cyto_str, '-' )[0]

	if offset_special.has_key (cyto):       # just match first band
                return offset_special [cyto]

        if offset_cre.match (cyto) != -1:
                pq, value = offset_cre.group (1,2)
                if pq == 'p':
                        factor = -1
                else:
                        factor = 1
                try:
                        float_value = string.atof (value)
                        return factor * (float_value + 5)
                except:
                        pass
        return 0.0

def processSort1():

	reportTitle = 'Homology - Rat vs. Human (Sorted by Rat Chromosome)'
	reportName = REPORTNAME + '1'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'])
	
	fp.write(string.ljust('Rat Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Rat LocusLink ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Rat Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Human Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Human LocusLink ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Human Symbol', 25))
	fp.write(SPACE)
	fp.write(CRT)

	fp.write(string.ljust('---------', 15))
	fp.write(SPACE)
	fp.write(string.ljust('------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('---------', 15))
	fp.write(SPACE)
	fp.write(string.ljust('------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(CRT)

	cmd = 'select distinct ratChr = m1.chromosome + m1.cytogeneticOffset, ' + \
		      'c.sequenceNum, ' + \
		      'm1.cytogeneticOffset, ' + \
                      'ratSymbol = m1.symbol, ' + \
		      'humanChr = m2.chromosome  + m2.cytogeneticOffset, ' + \
                      'humanSymbol = m2.symbol, ' + \
		      'ratLL = ha.accID, ' + \
		      'humanLL = ma.accID ' + \
		'from HMD_Homology r1, HMD_Homology_Marker h1, ' + \
		'HMD_Homology r2, HMD_Homology_Marker h2, ' + \
		'MRK_Marker m1, MRK_Marker m2, MRK_Chromosome c, ' + \
		'ACC_Accession ha, ACC_Accession ma ' + \
		'where m1._Species_key = 40 ' + \
		'and m1._Marker_key = h1._Marker_key ' + \
		'and h1._Homology_key = r1._Homology_key ' + \
		'and r1._Class_key = r2._Class_key ' + \
		'and r2._Homology_key = h2._Homology_key ' + \
		'and h2._Marker_key = m2._Marker_key ' + \
		'and m2._Species_key = 2 ' + \
		'and m1._Species_key = c._Species_key ' + \
		'and m1.chromosome = c.chromosome ' + \
		'and m1._Marker_key *= ha._Object_key ' + \
		'and ha._MGIType_key = 2 ' + \
		'and ha._LogicalDB_key = 24 ' + \
		'and m2._Marker_key *= ma._Object_key ' + \
		'and ma._MGIType_key = 2 ' + \
		'and ma._LogicalDB_key = 24 '

	results = db.sql(cmd, 'auto')
	count = 0

	#
	# initialize a list to sort the rat chromosome & offset values.
	# the first sort is by chromosome, second is by offset.
	# that is, sortKeys[0] holds the sequence number of the chromosome order
	# and sortKeys[1] holds the sort value of the cytogenetic offset
	#
	# store the sort key as a tuple in a dictionary (rows) so we can sort 
	# the dictionary keys
	#
	# the dictionary values will be set to the row tuple
	#

	sortKeys = [''] * 3	# initialize list to 3 'blanks'
	rows = {}

	for r in results:
		sortKeys[0] = r['sequenceNum']
		sortKeys[1] = getSortableOffset(r['cytogeneticOffset'])
		sortKeys[2] = r['ratSymbol']
		rows[tuple(sortKeys)] = r
		count = count + 1

	#
	# now sort the "rows" dictionary keys
	# and print out the dictionary values
	#

	keys = rows.keys()
	keys.sort()
	for key in keys:
		r = rows[key]
		fp.write(string.ljust(r['ratChr'], 15))
		fp.write(SPACE)
		fp.write(string.ljust(mgi_utils.prvalue(r['ratLL']), 30))
		fp.write(SPACE)
		fp.write(string.ljust(r['ratSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['humanChr'], 15))
		fp.write(SPACE)
		fp.write(string.ljust(mgi_utils.prvalue(r['humanLL']), 30))
		fp.write(SPACE)
		fp.write(string.ljust(r['humanSymbol'], 25))
		fp.write(SPACE)
		fp.write(CRT)

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.trailer(fp)
	reportlib.finish_nonps(fp)

def processSort2():

	reportTitle = 'Homology - Rat vs. Human (Sorted by Human Chromosome)'
	reportName = REPORTNAME + '2'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'])
	
	fp.write(string.ljust('Human Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Human LocusLink ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Human Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Rat Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Rat LocusLink ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Rat Symbol', 25))
	fp.write(SPACE)
	fp.write(CRT)

	fp.write(string.ljust('----------', 15))
	fp.write(SPACE)
	fp.write(string.ljust('------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 15))
	fp.write(SPACE)
	fp.write(string.ljust('------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(CRT)

	cmd = 'select distinct humanChr = m2.chromosome + m2.cytogeneticOffset, ' + \
	      'c.sequenceNum, ' + \
	      'm2.cytogeneticOffset, ' + \
              'humanSymbol = m2.symbol, ' + \
              'ratChr = m1.chromosome + m1.cytogeneticOffset, ' + \
              'ratSymbol = m1.symbol, ' + \
              'ratLL = ha.accID, ' + \
              'humanLL = ma.accID ' + \
              'from HMD_Homology r1, HMD_Homology_Marker h1, ' + \
              'HMD_Homology r2, HMD_Homology_Marker h2, ' + \
              'MRK_Marker m1, MRK_Marker m2, MRK_Chromosome c, ' + \
	      'ACC_Accession ha, ACC_Accession ma ' + \
              'where m1._Species_key = 40 ' + \
              'and m1._Marker_key = h1._Marker_key ' + \
              'and h1._Homology_key = r1._Homology_key ' + \
              'and r1._Class_key = r2._Class_key ' + \
              'and r2._Homology_key = h2._Homology_key ' + \
              'and h2._Marker_key = m2._Marker_key ' + \
              'and m2._Species_key = 2 ' + \
              'and m2._Species_key = c._Species_key ' + \
              'and m2.chromosome = c.chromosome ' + \
	      'and m1._Marker_key *= ha._Object_key ' + \
	      'and ha._MGIType_key = 2 ' + \
	      'and ha._LogicalDB_key = 24 ' + \
	      'and m2._Marker_key *= ma._Object_key ' + \
	      'and ma._MGIType_key = 2 ' + \
	      'and ma._LogicalDB_key = 24'

	results = db.sql(cmd, 'auto')
	count = 0

	#
	# initialize a list to sort the rat chromosome & offset values.
	# the first sort is by chromosome, second is by offset.
	# that is, sortKeys[0] holds the sequence number of the chromosome order
	# and sortKeys[1] holds the sort value of the cytogenetic offset
	#
	# store the sort key as a tuple in a dictionary (rows) so we can sort 
	# the dictionary keys
	#
	# the dictionary values will be set to the row tuple
	#

	sortKeys = [''] * 3	# initialize list to 3 'blanks'
	rows = {}

	for r in results:
		sortKeys[0] = r['sequenceNum']
		sortKeys[1] = getSortableOffset(r['cytogeneticOffset'])
		sortKeys[2] = r['humanSymbol']
		rows[tuple(sortKeys)] = r
		count = count + 1

	#
	# now sort the "rows" dictionary keys
	# and print out the dictionary values
	#

	keys = rows.keys()
	keys.sort()
	for key in keys:
		r = rows[key]
		fp.write(string.ljust(r['humanChr'], 15))
		fp.write(SPACE)
		fp.write(string.ljust(mgi_utils.prvalue(r['humanLL']), 30))
		fp.write(SPACE)
		fp.write(string.ljust(r['humanSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['ratChr'], 15))
		fp.write(SPACE)
		fp.write(string.ljust(mgi_utils.prvalue(r['ratLL']), 30))
		fp.write(SPACE)
		fp.write(string.ljust(r['ratSymbol'], 25))
		fp.write(CRT)

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.trailer(fp)
	reportlib.finish_nonps(fp)

def processSort3():

	reportTitle = 'Homology - Rat vs. Human (Sorted by Rat Symbol)'
	reportName = REPORTNAME + '3'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'])
	
	fp.write(string.ljust('Rat LocusLink ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Rat Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Rat Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Human LocusLink ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Human Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Human Chr', 15))
	fp.write(CRT)

	fp.write(string.ljust('------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 15))
	fp.write(SPACE)
	fp.write(string.ljust('------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 15))
	fp.write(CRT)

	cmd = 'select distinct humanChr = m2.chromosome + m2.cytogeneticOffset, ' + \
              'humanSymbol = m2.symbol, ' + \
              'ratChr = m1.chromosome + m1.cytogeneticOffset, ' + \
              'ratSymbol = m1.symbol, ' + \
              'ratLL = ha.accID, ' + \
              'humanLL = ma.accID ' + \
              'from HMD_Homology r1, HMD_Homology_Marker h1, ' + \
              'HMD_Homology r2, HMD_Homology_Marker h2, ' + \
              'MRK_Marker m1, MRK_Marker m2, ' + \
	      'ACC_Accession ha, ACC_Accession ma ' + \
              'where m1._Species_key = 40 ' + \
              'and m1._Marker_key = h1._Marker_key ' + \
              'and h1._Homology_key = r1._Homology_key ' + \
              'and r1._Class_key = r2._Class_key ' + \
              'and r2._Homology_key = h2._Homology_key ' + \
              'and h2._Marker_key = m2._Marker_key ' + \
              'and m2._Species_key = 2 ' + \
	      'and m1._Marker_key *= ha._Object_key ' + \
	      'and ha._MGIType_key = 2 ' + \
	      'and ha._LogicalDB_key = 24 ' + \
	      'and m2._Marker_key *= ma._Object_key ' + \
	      'and ma._MGIType_key = 2 ' + \
	      'and ma._LogicalDB_key = 24 ' + \
              'order by m1.symbol, m1.chromosome, m1.cytogeneticOffset'

	results = db.sql(cmd, 'auto')
	count = 0

	for r in results:
		fp.write(string.ljust(mgi_utils.prvalue(r['ratLL']), 30))
		fp.write(SPACE)
		fp.write(string.ljust(r['ratSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['ratChr'], 15))
		fp.write(SPACE)
		fp.write(string.ljust(mgi_utils.prvalue(r['humanLL']), 30))
		fp.write(SPACE)
		fp.write(string.ljust(r['humanSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['humanChr'], 15))
		fp.write(CRT)
		count = count + 1

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.trailer(fp)
	reportlib.finish_nonps(fp)

def processSort4():

	reportTitle = 'Homology - Rat vs. Human (Sorted by Human Symbol)'
	reportName = REPORTNAME + '4'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'])
	
	fp.write(string.ljust('Human LocusLink ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Human Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Human Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Rat LocusLink ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Rat Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Rat Chr', 15))
	fp.write(SPACE)
	fp.write(CRT)

	fp.write(string.ljust('------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 15))
	fp.write(SPACE)
	fp.write(string.ljust('------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 15))
	fp.write(SPACE)
	fp.write(CRT)

	cmd = 'select distinct humanChr = m2.chromosome + m2.cytogeneticOffset, ' + \
              'humanSymbol = m2.symbol, ' + \
              'ratChr = m1.chromosome + m1.cytogeneticOffset, ' + \
              'ratSymbol = m1.symbol, ' + \
              'ratLL = ha.accID, ' + \
              'humanLL = ma.accID ' + \
              'from HMD_Homology r1, HMD_Homology_Marker h1, ' + \
              'HMD_Homology r2, HMD_Homology_Marker h2, ' + \
              'MRK_Marker m1, MRK_Marker m2, ' + \
	      'ACC_Accession ha, ACC_Accession ma ' + \
              'where m1._Species_key = 40 ' + \
              'and m1._Marker_key = h1._Marker_key ' + \
              'and h1._Homology_key = r1._Homology_key ' + \
              'and r1._Class_key = r2._Class_key ' + \
              'and r2._Homology_key = h2._Homology_key ' + \
              'and h2._Marker_key = m2._Marker_key ' + \
              'and m2._Species_key = 2 ' + \
	      'and m1._Marker_key *= ha._Object_key ' + \
	      'and ha._MGIType_key = 2 ' + \
	      'and ha._LogicalDB_key = 24 ' + \
	      'and m2._Marker_key *= ma._Object_key ' + \
	      'and ma._MGIType_key = 2 ' + \
	      'and ma._LogicalDB_key = 24 ' + \
              'order by m2.symbol, m2.chromosome, m2.cytogeneticOffset'

	results = db.sql(cmd, 'auto')
	count = 0

	for r in results:
		fp.write(string.ljust(mgi_utils.prvalue(r['humanLL']), 30))
		fp.write(SPACE)
		fp.write(string.ljust(r['humanSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['humanChr'], 15))
		fp.write(SPACE)
		fp.write(string.ljust(mgi_utils.prvalue(r['ratLL']), 30))
		fp.write(SPACE)
		fp.write(string.ljust(r['ratSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['ratChr'], 15))
		fp.write(CRT)
		count = count + 1

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.trailer(fp)
	reportlib.finish_nonps(fp)

#
# Main
#

sortOption = None

if len(sys.argv) > 1:
	sortOption = sys.argv[1]

if sortOption == '1':
	processSort1()
elif sortOption == '2':
	processSort2()
elif sortOption == '3':
	processSort3()
elif sortOption == '4':
	processSort4()
else:
	processSort1()
	processSort2()
	processSort3()
	processSort4()
	
