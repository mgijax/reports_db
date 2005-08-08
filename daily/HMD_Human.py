#!/usr/local/bin/python

'''
#
# HMD_Human.py 04/20/2000
#
# Report:
#	TR 211
#	Originally 4 separate .sql reports (HMD_Human1.sql, etc.)
#	Re-implemeneted in Python due to sorting of human cytogenetic offsets.
#	The offset sorting algorithm is duplicated in the WI homology_report.cgi.py
#	module.
#
# Usage:
#       HMD_Human.py [1234]
#	1 = sort by Human Chromosome
#	2 = sort by Mouse Chromosome
#	3 = sort by Human Symbol
#	4 = sort by Mouse Symbol
#	none = generate AEG 4 reports
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
# lec	01/04/2004
#	- TR 5939; LocusLink->EntrezGene
#
# lec	07/23/2004
#	- TR 5611; added curated/calculated
#
# lec	07/01/2003
#	- TR 4945; added EntrezGene IDs for Mouse and Human
#
# lec   07/18/2000
#	- TR 1806; added MGI Accession ID for Mouse Symbols
#
# lec   04/20/2000
#       - created
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
TAB = reportlib.TAB
REPORTNAME = 'HMD_Human'

reportLegend = 'Data Attributes:  M - MGI curated, C - HomoloGene calculated, B - MGI curated and HomoloGene calculated'
mgiCurated = 'M'
hgCurated = 'C'
bothCurated = 'B'

curated = []
calculated = []
humanEG = {}
mouseEG = {}
mouseMGI = {}

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
        #       which will sort as desired.  (using human as a primary
        #       organism requires special handling, as of TR 211)
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

def runQueries():

	global curated, calculated, humanEG, mouseEG, mouseMGI

	cmd = 'select distinct h1._Marker_key ' + \
		'from HMD_Homology r1, HMD_Homology_Marker h1, ' + \
		'HMD_Homology r2, HMD_Homology_Marker h2, ' + \
		'MRK_Marker m1, MRK_Marker m2 ' + \
		'where r1._Refs_key = 91485 ' + \
		'and m1._Organism_key = 2 ' + \
		'and m1._Marker_key = h1._Marker_key ' + \
		'and h1._Homology_key = r1._Homology_key ' + \
		'and r1._Class_key = r2._Class_key ' + \
		'and r2._Homology_key = h2._Homology_key ' + \
		'and h2._Marker_key = m2._Marker_key ' + \
		'and m2._Organism_key = 1'

	results = db.sql(cmd, 'auto')
	for r in results:
		calculated.append(r['_Marker_key'])

	cmd = 'select distinct h1._Marker_key ' + \
		'from HMD_Homology r1, HMD_Homology_Marker h1, ' + \
		'HMD_Homology r2, HMD_Homology_Marker h2, ' + \
		'MRK_Marker m1, MRK_Marker m2 ' + \
		'where r1._Refs_key != 91485 ' + \
		'and m1._Organism_key = 2 ' + \
		'and m1._Marker_key = h1._Marker_key ' + \
		'and h1._Homology_key = r1._Homology_key ' + \
		'and r1._Class_key = r2._Class_key ' + \
		'and r2._Homology_key = h2._Homology_key ' + \
		'and h2._Marker_key = m2._Marker_key ' + \
		'and m2._Organism_key = 1'

	results = db.sql(cmd, 'auto')
	for r in results:
    		curated.append(r['_Marker_key'])

	##

	cmds = []
	cmds.append('select distinct humanMarkerKey = h1._Marker_key, mouseMarkerKey = h2._Marker_key ' + \
	'into #allhomologies ' + \
      	'from HMD_Homology r1, HMD_Homology_Marker h1, ' + \
      	'HMD_Homology r2, HMD_Homology_Marker h2, ' + \
      	'MRK_Marker m1, MRK_Marker m2 ' + \
      	'where m1._Organism_key = 2 ' + \
      	'and m1._Marker_key = h1._Marker_key ' + \
      	'and h1._Homology_key = r1._Homology_key ' + \
      	'and r1._Class_key = r2._Class_key ' + \
      	'and r2._Homology_key = h2._Homology_key ' + \
      	'and h2._Marker_key = m2._Marker_key ' + \
      	'and m2._Organism_key = 1')

	cmds.append('create nonclustered index idx_hkey on #allhomologies(humanMarkerKey)')
	cmds.append('create nonclustered index idx_mkey on #allhomologies(mouseMarkerKey)')

	db.sql(cmds, None)

	##

	cmds = []
	cmds.append('select h.humanMarkerKey, h.mouseMarkerKey, ' + \
		'humanOrganism = m1._Organism_key, ' + \
		'humanSymbol = m1.symbol, ' + \
		'humanChr = m1.chromosome + m1.cytogeneticOffset, ' + \
		'm1.chromosome, ' + \
		'm1.cytogeneticOffset, ' + \
		'mouseOrganism = m2._Organism_key, ' + \
		'mouseSymbol = m2.symbol, ' + \
		'mouseChr = m2.chromosome, ' + \
		'mouseBand = m2.cytogeneticOffset, ' + \
		'mouseName = substring(m2.name, 1, 75), ' + \
                'mouseCm = ' + \
        	'case ' + \
        	'when o.offset >= 0 then str(o.offset, 10, 2) ' + \
        	'when o.offset = -999.0 then "       N/A" ' + \
        	'when o.offset = -1.0 then "  syntenic" ' + \
        	'end,' + \
		'mouseOffset = o.offset ' + \
		'into #homologies ' + \
		'from #allhomologies h, MRK_Marker m1, MRK_Marker m2, MRK_Offset o ' + \
		'where h.humanMarkerKey = m1._Marker_key ' + \
		'and h.mouseMarkerKey = m2._Marker_key ' + \
		'and h.mouseMarkerKey = o._Marker_key ' + \
		'and o.source = 0 ')

	cmds.append('create nonclustered index idx_hkey1 on #homologies(humanOrganism)')
	cmds.append('create nonclustered index idx_mkey1 on #homologies(mouseOrganism)')
	cmds.append('create nonclustered index idx_hkey2 on #homologies(humanSymbol)')
	cmds.append('create nonclustered index idx_mkey2 on #homologies(mouseSymbol)')

	db.sql(cmds, None)

	##

	# human entrezgene ids

	results = db.sql('select h.humanMarkerKey, a.accID from #homologies h, ACC_Accession a ' + \
		'where h.humanMarkerKey = a._Object_key ' + \
		'and a._MGIType_key = 2 ' + \
		'and a._LogicalDB_key = 55 ', 'auto')
	for r in results:
		humanEG[r['humanMarkerKey']] = r['accID']

	# mouse entrezgene ids

	results = db.sql('select h.mouseMarkerKey, a.accID from #homologies h, ACC_Accession a ' + \
		'where h.mouseMarkerKey = a._Object_key ' + \
		'and a._MGIType_key = 2 ' + \
		'and a._LogicalDB_key = 55 ', 'auto')
	for r in results:
		mouseEG[r['mouseMarkerKey']] = r['accID']

	# mouse MGI 

	results = db.sql('select h.mouseMarkerKey, a.accID from #homologies h, ACC_Accession a ' + \
		'where h.mouseMarkerKey = a._Object_key ' + \
		'and a._MGIType_key = 2 ' + \
		'and a.prefixPart = "MGI:" ' + \
		'and a._LogicalDB_key = 1 ' + \
		'and a.preferred = 1 ', 'auto')
	for r in results:
		mouseMGI[r['mouseMarkerKey']] = r['accID']

	##

	# sorted by human chromosome

	results1 = db.sql('select h.*, c.sequenceNum ' + \
		'from #homologies h, MRK_Chromosome c ' + \
		'where h.humanOrganism = c._Organism_key ' + \
		'and h.chromosome = c.chromosome ', 'auto')

	# sorted by mouse chromosome

	results2 = db.sql('select h.*, c.sequenceNum ' + \
		'from #homologies h, MRK_Chromosome c ' + \
		'where h.mouseOrganism = c._Organism_key ' + \
		'and h.mouseChr = c.chromosome ' + \
		'order by c.sequenceNum, h.mouseOffset', 'auto')

	# sorted by human symbol

	results3 = db.sql('select * from #homologies order by humanSymbol', 'auto')

	# sorted by mouse symbol

	results4 = db.sql('select * from #homologies order by mouseSymbol', 'auto')
	
	return results1, results2, results3, results4

def printDataAttributes(fp, key):

	if key in curated and not key in calculated:
  		fp.write(mgiCurated)
	elif key not in curated and key in calculated:
  		fp.write(hgCurated)
	else:
  		fp.write(bothCurated)

def processSort1(results):

	reportTitle = 'Orthology - Human vs. Mouse (Sorted by Human Chromosome)'
	reportName = REPORTNAME + '1'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'])
	fp.write(reportLegend + CRT + CRT)

	fp.write(string.ljust('Human Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Human EntrezGene ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Human Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse MGI Acc ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Chr', 10))
	fp.write(SPACE)
	fp.write(string.ljust('  Mouse cM', 10))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse EntrezGene ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Name', 80))
	fp.write(SPACE)
	fp.write(string.ljust('Data Attributes', 15))
	fp.write(SPACE)
	fp.write(CRT)

	fp.write(string.ljust('---------', 15))
	fp.write(SPACE)
	fp.write(string.ljust('-------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('---------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('  --------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('-------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 80))
	fp.write(SPACE)
	fp.write(string.ljust('---------------', 15))
	fp.write(SPACE)
	fp.write(CRT)

	#
	# initialize a list to sort the human chromosome & offset values.
	# the first sort is by chromosome, second is by offset.
	# that is, sortKeys[0] holds the sequence number of the chromosome order
	# and sortKeys[1] holds the sort value of the cytogenetic offset
	#
	# store the sort key as a tuple in a dictionary (rows) so we can sort 
	# the dictionary keys
	#
	# the dictionary values will be set to the row tuple
	#

	count = 0
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

		if humanEG.has_key(r['humanMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(humanEG[r['humanMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['humanSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(mouseMGI[r['mouseMarkerKey']], 30))
		fp.write(SPACE)
		fp.write(string.ljust(r['mouseChr'], 10))
		fp.write(SPACE)
		fp.write(string.ljust(r['mouseCm'], 10))
		fp.write(SPACE)

		if mouseEG.has_key(r['mouseMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(mouseEG[r['mouseMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['mouseSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['mouseName'], 80))
		fp.write(SPACE)
		printDataAttributes(fp, r['humanMarkerKey'])
		fp.write(CRT)

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.trailer(fp)
	reportlib.finish_nonps(fp)

def processSort2(results):

	reportTitle = 'Orthology - Human vs. Mouse (Sorted by Mouse Chromosome)'
	reportName = REPORTNAME + '2'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'])
	fp.write(reportLegend + CRT + CRT)
	
	fp.write(string.ljust('Mouse MGI Acc ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Chr', 10))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse cM', 10))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse EntrezGene ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Human Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Human EntrezGene ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Human Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Name', 80))
	fp.write(SPACE)
	fp.write(string.ljust('Data Attributes', 15))
	fp.write(SPACE)
	fp.write(CRT)

	fp.write(string.ljust('----------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('-------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 15))
	fp.write(SPACE)
	fp.write(string.ljust('-------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 80))
	fp.write(SPACE)
	fp.write(string.ljust('---------------', 15))
	fp.write(SPACE)
	fp.write(CRT)

	count = 0

	for r in results:
		fp.write(string.ljust(mouseMGI[r['mouseMarkerKey']], 30))
		fp.write(SPACE)
		fp.write(string.ljust(r['mouseChr'], 10))
		fp.write(SPACE)
		fp.write(string.ljust(r['mouseCm'], 10))
		fp.write(SPACE)

		if mouseEG.has_key(r['mouseMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(mouseEG[r['mouseMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['mouseSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['humanChr'], 15))
		fp.write(SPACE)

		if humanEG.has_key(r['humanMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(humanEG[r['humanMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['humanSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['mouseName'], 80))
		fp.write(SPACE)
		printDataAttributes(fp, r['humanMarkerKey'])
		fp.write(CRT)
		count = count + 1

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.trailer(fp)
	reportlib.finish_nonps(fp)

def processSort3(results):

	reportTitle = 'Orthology - Human vs. Mouse (Sorted by Human Symbol)'
	reportName = REPORTNAME + '3'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'])
	fp.write(reportLegend + CRT + CRT)
	
	fp.write(string.ljust('Human EntrezGene ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Human Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Human Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse MGI Acc ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse EntrezGene ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Chr', 10))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse cM', 10))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Band', 10))
	fp.write(SPACE)
	fp.write(string.ljust('Data Attributes', 15))
	fp.write(SPACE)
	fp.write(CRT)

	fp.write(string.ljust('-------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 15))
	fp.write(SPACE)
	fp.write(string.ljust('----------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('-------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('---------------', 15))
	fp.write(SPACE)
	fp.write(CRT)

	count = 0

	for r in results:
		if humanEG.has_key(r['humanMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(humanEG[r['humanMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['humanSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['humanChr'], 15))
		fp.write(SPACE)
		fp.write(string.ljust(mouseMGI[r['mouseMarkerKey']], 30))
		fp.write(SPACE)

		if mouseEG.has_key(r['mouseMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(mouseEG[r['mouseMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['mouseSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['mouseChr'], 10))
		fp.write(SPACE)
		fp.write(string.ljust(r['mouseCm'], 10))
		fp.write(SPACE)
		fp.write(string.ljust(mgi_utils.prvalue(r['mouseBand']), 10))
		fp.write(SPACE)
		printDataAttributes(fp, r['humanMarkerKey'])
		fp.write(CRT)
		count = count + 1

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.trailer(fp)
	reportlib.finish_nonps(fp)

def processSort4(results):

	reportTitle = 'Orthology - Human vs. Mouse (Sorted by Mouse Symbol)'
	reportName = REPORTNAME + '4'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'])
	fp.write(reportLegend + CRT + CRT)
	
	fp.write(string.ljust('Mouse MGI Acc ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse EntrezGene ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Chr', 10))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse cM', 10))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Band', 10))
	fp.write(SPACE)
	fp.write(string.ljust('Human EntrezGene ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Human Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Human Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Data Attributes', 15))
	fp.write(SPACE)
	fp.write(CRT)

	fp.write(string.ljust('----------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('-------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('-------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 15))
	fp.write(SPACE)
	fp.write(string.ljust('---------------', 15))
	fp.write(SPACE)
	fp.write(CRT)

	count = 0

	for r in results:
		fp.write(string.ljust(mouseMGI[r['mouseMarkerKey']], 30))
		fp.write(SPACE)

		if mouseEG.has_key(r['mouseMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(mouseEG[r['mouseMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['mouseSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['mouseChr'], 10))
		fp.write(SPACE)
		fp.write(string.ljust(r['mouseCm'], 10))
		fp.write(SPACE)
		fp.write(string.ljust(mgi_utils.prvalue(r['mouseBand']), 10))
		fp.write(SPACE)

		if humanEG.has_key(r['humanMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(humanEG[r['humanMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['humanSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['humanChr'], 15))
		fp.write(SPACE)
		printDataAttributes(fp, r['humanMarkerKey'])
		fp.write(CRT)
		count = count + 1

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.trailer(fp)
	reportlib.finish_nonps(fp)

def processSort5(results):

	reportTitle = 'Orthology - Human vs. Mouse (Sorted by Human Chromosome)'
	reportName = REPORTNAME + '5'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'])
	fp.write(reportLegend + CRT + CRT)

	fp.write('Human Chr' + TAB)
	fp.write('Human EntrezGene ID' + TAB)
	fp.write('Human Symbol' + TAB)
	fp.write('Mouse MGI Acc ID' + TAB)
	fp.write('Mouse Chr' + TAB)
	fp.write('Mouse cM' + TAB)
	fp.write('Mouse EntrezGene ID' + TAB)
	fp.write('Mouse Symbol' + TAB)
	fp.write('Mouse Name' + TAB)
	fp.write('Data Attributes' + CRT * 2)

	#
	# initialize a list to sort the human chromosome & offset values.
	# the first sort is by chromosome, second is by offset.
	# that is, sortKeys[0] holds the sequence number of the chromosome order
	# and sortKeys[1] holds the sort value of the cytogenetic offset
	#
	# store the sort key as a tuple in a dictionary (rows) so we can sort 
	# the dictionary keys
	#
	# the dictionary values will be set to the row tuple
	#

	count = 0
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
		fp.write(r['humanChr'] + TAB)

		if humanEG.has_key(r['humanMarkerKey']):
			fp.write(mgi_utils.prvalue(humanEG[r['humanMarkerKey']]))
		fp.write(TAB)
		fp.write(r['humanSymbol'] + TAB)
		fp.write(mouseMGI[r['mouseMarkerKey']] + TAB)
		fp.write(r['mouseChr'] + TAB)
		fp.write(r['mouseCm'] + TAB)

		if mouseEG.has_key(r['mouseMarkerKey']):
			fp.write(mgi_utils.prvalue(mouseEG[r['mouseMarkerKey']]) + TAB)
		fp.write(TAB)

		fp.write(r['mouseSymbol'] + TAB)
		fp.write(r['mouseName'] + TAB)
		printDataAttributes(fp, r['humanMarkerKey'])
		fp.write(CRT)

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.trailer(fp)
	reportlib.finish_nonps(fp)

#
# Main
#

db.useOneConnection(1)
sortOption = None

if len(sys.argv) > 1:
	sortOption = sys.argv[1]

r1, r2, r3, r4 = runQueries()

if sortOption == '1':
	processSort1(r1)
elif sortOption == '2':
	processSort2(r2)
elif sortOption == '3':
	processSort3(r3)
elif sortOption == '4':
	processSort4(r4)
elif sortOption == '5':
	processSort5(r1)
else:
	processSort1(r1)
	processSort2(r2)
	processSort3(r3)
	processSort4(r4)
	processSort5(r1)

db.useOneConnection(0)
