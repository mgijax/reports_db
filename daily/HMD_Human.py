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
# lec	07/23/2004
#	- TR 5611; added curated/calculated
#
# lec	07/01/2003
#	- TR 4945; added LocusLink IDs for Mouse and Human
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
REPORTNAME = 'HMD_Human'

reportLegend = 'Data Attributes:  C = Curated, c = calculated'
isCurated = 'C'
isCalculated = 'c'

curated = []
calculated = []
humanLL = {}
mouseLL = {}
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

	global curated, calculated, humanLL, mouseLL, mouseMGI

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

	# human locus link ids

	cmds.append('select h.humanMarkerKey, a.accID from #homologies h, ACC_Accession a ' + \
		'where h.humanMarkerKey = a._Object_key ' + \
		'and a._MGIType_key = 2 ' + \
		'and a._LogicalDB_key = 24 ')

	# mouse locus link ids

	cmds.append('select h.mouseMarkerKey, a.accID from #homologies h, ACC_Accession a ' + \
		'where h.mouseMarkerKey = a._Object_key ' + \
		'and a._MGIType_key = 2 ' + \
		'and a._LogicalDB_key = 24 ')

	# mouse MGI 

	cmds.append('select h.mouseMarkerKey, a.accID from #homologies h, ACC_Accession a ' + \
		'where h.mouseMarkerKey = a._Object_key ' + \
		'and a._MGIType_key = 2 ' + \
		'and a.prefixPart = "MGI:" ' + \
		'and a._LogicalDB_key = 1 ' + \
		'and a.preferred = 1 ')

	# sorted by human chromosome

	cmds.append('select h.*, c.sequenceNum ' + \
		'from #homologies h, MRK_Chromosome c ' + \
		'where h.humanOrganism = c._Organism_key ' + \
		'and h.chromosome = c.chromosome ')

	# sorted by mouse chromosome

	cmds.append('select h.*, c.sequenceNum ' + \
		'from #homologies h, MRK_Chromosome c ' + \
		'where h.mouseOrganism = c._Organism_key ' + \
		'and h.mouseChr = c.chromosome ' + \
		'order by c.sequenceNum, h.mouseOffset')

	# sorted by human symbol

	cmds.append('select * from #homologies order by humanSymbol')

	# sorted by mouse symbol

	cmds.append('select * from #homologies order by mouseSymbol')
	
	results = db.sql(cmds, 'auto')

	for r in results[8]:
		humanLL[r['humanMarkerKey']] = r['accID']

	for r in results[9]:
		mouseLL[r['mouseMarkerKey']] = r['accID']

	for r in results[10]:
		mouseMGI[r['mouseMarkerKey']] = r['accID']

	return results

def processSort1(results):

	reportTitle = 'Orthology - Human vs. Mouse (Sorted by Human Chromosome)'
	reportName = REPORTNAME + '1'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'])
	fp.write(reportLegend + CRT + CRT)

	fp.write(string.ljust('Human Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Human LocusLink ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Human Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse MGI Acc ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Chr', 10))
	fp.write(SPACE)
	fp.write(string.ljust('  Mouse cM', 10))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse LocusLink ID', 30))
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
	fp.write(string.ljust('------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('---------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('  --------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('------------------', 30))
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

		if humanLL.has_key(r['humanMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(humanLL[r['humanMarkerKey']]), 30))
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

		if mouseLL.has_key(r['mouseMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(mouseLL[r['mouseMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['mouseSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['mouseName'], 80))
		fp.write(SPACE)

		if r['humanMarkerKey'] in curated:
		  fp.write(isCurated)
		if r['humanMarkerKey'] in calculated:
		  fp.write(isCalculated)

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
	fp.write(string.ljust('Mouse LocusLink ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Human Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Human LocusLink ID', 30))
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

		if mouseLL.has_key(r['mouseMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(mouseLL[r['mouseMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['mouseSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['humanChr'], 15))
		fp.write(SPACE)

		if humanLL.has_key(r['humanMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(humanLL[r['humanMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['humanSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['mouseName'], 80))
		fp.write(SPACE)

		if r['humanMarkerKey'] in curated:
		  fp.write(isCurated)
		if r['humanMarkerKey'] in calculated:
		  fp.write(isCalculated)

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
	
	fp.write(string.ljust('Human LocusLink ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Human Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Human Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse MGI Acc ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse LocusLink ID', 30))
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

	fp.write(string.ljust('------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 15))
	fp.write(SPACE)
	fp.write(string.ljust('----------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------------', 30))
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
		if humanLL.has_key(r['humanMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(humanLL[r['humanMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['humanSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['humanChr'], 15))
		fp.write(SPACE)
		fp.write(string.ljust(mouseMGI[r['mouseMarkerKey']], 30))
		fp.write(SPACE)

		if mouseLL.has_key(r['mouseMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(mouseLL[r['mouseMarkerKey']]), 30))
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

		if r['humanMarkerKey'] in curated:
		  fp.write(isCurated)
		if r['humanMarkerKey'] in calculated:
		  fp.write(isCalculated)

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
	fp.write(string.ljust('Mouse LocusLink ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Chr', 10))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse cM', 10))
	fp.write(SPACE)
	fp.write(string.ljust('Mouse Band', 10))
	fp.write(SPACE)
	fp.write(string.ljust('Human LocusLink ID', 30))
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
	fp.write(string.ljust('------------------', 30))
	fp.write(SPACE)
	fp.write(string.ljust('------------', 25))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('----------', 10))
	fp.write(SPACE)
	fp.write(string.ljust('------------------', 30))
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

		if mouseLL.has_key(r['mouseMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(mouseLL[r['mouseMarkerKey']]), 30))
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

		if humanLL.has_key(r['humanMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(humanLL[r['humanMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['humanSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['humanChr'], 15))
		fp.write(SPACE)

		if r['humanMarkerKey'] in curated:
		  fp.write(isCurated)
		if r['humanMarkerKey'] in calculated:
		  fp.write(isCalculated)

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

results = runQueries()

if sortOption == '1':
	processSort1(results[11])
elif sortOption == '2':
	processSort2(results[12])
elif sortOption == '3':
	processSort3(results[13])
elif sortOption == '4':
	processSort4(results[14])
else:
	processSort1(results[11])
	processSort2(results[12])
	processSort3(results[13])
	processSort4(results[14])
