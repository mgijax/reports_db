#!/usr/local/bin/python

'''
#
# HMD_Rat.py 04/20/2000
#
# Report:
#	TR 211
#	Originally 4 separate .sql reports (HMD_Rat1.sql, etc.)
#	Re-implemeneted in Python due to sorting of rat cytogenetic offsets.
#	The offset sorting algorithm is duplicated in the WI homology_report.cgi.py
#	module.
#
# Usage:
#       HMD_Rat.py [1234]
#	1 = sort by Rat Chromosome
#	2 = sort by Mouse Chromosome
#	3 = sort by Rat Symbol
#	4 = sort by Mouse Symbol
#	none = generate AEG 4 reports
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
#	- TR 4945; added LocusLink IDs for Mouse and Rat
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
import re
import mgi_utils
import reportlib

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
	db.setAutoTranslate(False)
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db


CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
REPORTNAME = 'HMD_Rat'

reportLegend = 'Data Attributes:  M - MGI curated, C - HomoloGene calculated, B - MGI curated and HomoloGene calculated'
mgiCurated = 'M'
hgCurated = 'C'
bothCurated = 'B'

curated = []
calculated = []
ratEG = {}
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

offset_cre = re.compile ('([pq])([0-9\.]+)')     # eg- q23.2
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

        offset_cre_result = offset_cre.match (cyto)
        if offset_cre_result is not None:
                pq, value = offset_cre_result.group (1,2)
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

	global curated, calculated, ratEG, mouseEG, mouseMGI

	cmd = '''
		select distinct h1._Marker_key 
		from MRK_Homology_Cache h1, MRK_Homology_Cache h2 
		where h1._Refs_key = 91485 
		and h1._Organism_key = 40 
		and h1._Class_key = h2._Class_key 
		and h2._Organism_key = 1
		'''

	results = db.sql(cmd, 'auto')
	for r in results:
		calculated.append(r['_Marker_key'])

	cmd = '''
		select distinct h1._Marker_key 
		from MRK_Homology_Cache h1, MRK_Homology_Cache h2 
		where h1._Refs_key != 91485 
		and h1._Organism_key = 40 
		and h1._Class_key = h2._Class_key 
		and h2._Organism_key = 1
		'''

	results = db.sql(cmd, 'auto')
	for r in results:
    		curated.append(r['_Marker_key'])

	db.sql('''
		select distinct h1._Marker_key as ratMarkerKey, h2._Marker_key as mouseMarkerKey
	        into #allhomologies 
	        from MRK_Homology_Cache h1, MRK_Homology_Cache h2 
	        where h1._Organism_key = 40 
	        and h1._Class_key = h2._Class_key 
	        and h2._Organism_key = 1
		''', None)

	db.sql('create index idx_hkey on #allhomologies(ratMarkerKey)', None)
	db.sql('create index idx_mkey on #allhomologies(mouseMarkerKey)', None)

	db.sql('''
		select h.ratMarkerKey, h.mouseMarkerKey, 
		m1._Organism_key as ratOrganism, 
		m1.symbol as ratSymbol, 
		m1.chromosome || m1.cytogeneticOffset as ratChr, 
		m1.chromosome, 
		m1.cytogeneticOffset, 
		m2._Organism_key as mouseOrganism, 
		m2.symbol as mouseSymbol, 
		m2.chromosome as mouseChr, 
		m2.cytogeneticOffset as mouseBand, 
		substring(m2.name, 1, 75) as mouseName, 
		o.offset as mouseOffset,
        	case 
        	when o.offset >= 0 then str(o.offset,10,2) 
        	when o.offset = -999.0 then "       N/A" 
        	when o.offset = -1.0 then "  syntenic" 
        	end as mouseCm
		into #homologies 
		from #allhomologies h, MRK_Marker m1, MRK_Marker m2, MRK_Offset o 
		where h.ratMarkerKey = m1._Marker_key 
		and h.mouseMarkerKey = m2._Marker_key 
		and h.mouseMarkerKey = o._Marker_key 
		and o.source = 0 
		''', None)

	db.sql('create index idx_hkey1 on #homologies(ratOrganism)', None)
	db.sql('create index idx_mkey1 on #homologies(mouseOrganism)', None)
	db.sql('create index idx_hkey2 on #homologies(ratSymbol)', None)
	db.sql('create index idx_mkey2 on #homologies(mouseSymbol)', None)

	# rat entrezgene ids

	results = db.sql('''
		select h.ratMarkerKey, a.accID
		from #homologies h, ACC_Accession a 
		where h.ratMarkerKey = a._Object_key 
		and a._MGIType_key = 2 
		and a._LogicalDB_key = 55 
		''', 'auto')

	for r in results:
		ratEG[r['ratMarkerKey']] = r['accID']

	# mouse entrezgene ids

	results = db.sql('''
		select h.mouseMarkerKey, a.accID
		from #homologies h, ACC_Accession a 
		where h.mouseMarkerKey = a._Object_key 
		and a._MGIType_key = 2 
		and a._LogicalDB_key = 55 
		''', 'auto')

	for r in results:
		mouseEG[r['mouseMarkerKey']] = r['accID']

	# mouse MGI 

	results = db.sql('''
		select h.mouseMarkerKey, a.accID
		from #homologies h, ACC_Accession a 
		where h.mouseMarkerKey = a._Object_key 
		and a._MGIType_key = 2 
		and a.prefixPart = 'MGI:' 
		and a._LogicalDB_key = 1 
		and a.preferred = 1 
		''', 'auto')

	for r in results:
		mouseMGI[r['mouseMarkerKey']] = r['accID']

	# sorted by rat chromosome

	results1 = db.sql('''
		select h.*, c.sequenceNum 
		from #homologies h, MRK_Chromosome c 
		where h.ratOrganism = c._Organism_key 
		and h.chromosome = c.chromosome 
		''', 'auto')

	# sorted by mouse chromosome

	results2 = db.sql('''
		select h.*, c.sequenceNum 
		from #homologies h, MRK_Chromosome c 
		where h.mouseOrganism = c._Organism_key 
		and h.mouseChr = c.chromosome 
		order by c.sequenceNum, h.mouseOffset
		''', 'auto')

	# sorted by rat symbol

	results3 = db.sql('select * from #homologies order by ratSymbol', 'auto')

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

	reportTitle = 'Orthology - Rat vs. Mouse (Sorted by Rat Chromosome)'
	reportName = REPORTNAME + '1'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'], printHeading = "MGI")
	fp.write(reportLegend + CRT + CRT)

	fp.write(string.ljust('Rat Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Rat EntrezGene ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Rat Symbol', 25))
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

	count = 0
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
		fp.write(string.ljust(mgi_utils.prvalue(r['ratChr']), 15))
		fp.write(SPACE)

		if ratEG.has_key(r['ratMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(ratEG[r['ratMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['ratSymbol'], 25))
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
		printDataAttributes(fp, r['ratMarkerKey'])
		fp.write(CRT)

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.finish_nonps(fp)

def processSort2(results):

	reportTitle = 'Orthology - Rat vs. Mouse (Sorted by Mouse Chromosome)'
	reportName = REPORTNAME + '2'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'], printHeading = "MGI")
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
	fp.write(string.ljust('Rat Chr', 15))
	fp.write(SPACE)
	fp.write(string.ljust('Rat EntrezGene ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Rat Symbol', 25))
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

		if mouseEG.has_key(r['mouseMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(mouseEG[r['mouseMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['mouseSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(mgi_utils.prvalue(r['ratChr']), 15))
		fp.write(SPACE)

		if ratEG.has_key(r['ratMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(ratEG[r['ratMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['ratSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['mouseName'], 80))
		fp.write(SPACE)
		printDataAttributes(fp, r['ratMarkerKey'])
		fp.write(CRT)
		count = count + 1

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.finish_nonps(fp)

def processSort3(results):

	reportTitle = 'Orthology - Rat vs. Mouse (Sorted by Rat Symbol)'
	reportName = REPORTNAME + '3'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'], printHeading = "MGI")
	fp.write(reportLegend + CRT + CRT)
	
	fp.write(string.ljust('Rat EntrezGene ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Rat Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Rat Chr', 15))
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
		if ratEG.has_key(r['ratMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(ratEG[r['ratMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['ratSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(mgi_utils.prvalue(r['ratChr']), 15))
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
		printDataAttributes(fp, r['ratMarkerKey'])
		fp.write(CRT)
		count = count + 1

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.finish_nonps(fp)

def processSort4(results):

	reportTitle = 'Orthology - Rat vs. Mouse (Sorted by Mouse Symbol)'
	reportName = REPORTNAME + '4'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'], printHeading = "MGI")
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
	fp.write(string.ljust('Rat EntrezGene ID', 30))
	fp.write(SPACE)
	fp.write(string.ljust('Rat Symbol', 25))
	fp.write(SPACE)
	fp.write(string.ljust('Rat Chr', 15))
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

		if ratEG.has_key(r['ratMarkerKey']):
			fp.write(string.ljust(mgi_utils.prvalue(ratEG[r['ratMarkerKey']]), 30))
		else:
			fp.write(string.ljust('', 30))

		fp.write(SPACE)
		fp.write(string.ljust(r['ratSymbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(mgi_utils.prvalue(r['ratChr']), 15))
		fp.write(SPACE)
		printDataAttributes(fp, r['ratMarkerKey'])
		fp.write(CRT)
		count = count + 1

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.finish_nonps(fp)

def processSort5(results):

	# same as processSort1, but tab-delimited

	reportTitle = 'Orthology - Rat vs. Mouse (Sorted by Rat Chromosome), tab-delimited'
	reportName = REPORTNAME + '5'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'], printHeading = "MGI")
	fp.write(reportLegend + CRT + CRT)

	fp.write('Rat Chr' + TAB)
	fp.write('Rat EntrezGene ID' + TAB)
	fp.write('Rat Symbol' + TAB)
	fp.write('Mouse MGI Acc ID' + TAB)
	fp.write('Mouse Chr' + TAB)
	fp.write('Mouse cM' + TAB)
	fp.write('Mouse EntrezGene ID' + TAB)
	fp.write('Mouse Symbol' + TAB)
	fp.write('Mouse Name' + TAB)
	fp.write('Data Attributes' + CRT * 2)

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

	#
	# now sort the "rows" dictionary keys
	# and print out the dictionary values
	#

	keys = rows.keys()
	keys.sort()
	for key in keys:
		r = rows[key]
		fp.write(mgi_utils.prvalue(r['ratChr']) + TAB)

		if ratEG.has_key(r['ratMarkerKey']):
			fp.write(mgi_utils.prvalue(ratEG[r['ratMarkerKey']]))
		fp.write(TAB)

		fp.write(r['ratSymbol'] + TAB)
		fp.write(mouseMGI[r['mouseMarkerKey']] + TAB)
		fp.write(r['mouseChr'] + TAB)
		fp.write(r['mouseCm'] + TAB)

		if mouseEG.has_key(r['mouseMarkerKey']):
			fp.write(mgi_utils.prvalue(mouseEG[r['mouseMarkerKey']]))
		fp.write(TAB)

		fp.write(r['mouseSymbol'] + TAB)
		fp.write(r['mouseName'] + TAB)
		printDataAttributes(fp, r['ratMarkerKey'])
		fp.write(CRT)

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
