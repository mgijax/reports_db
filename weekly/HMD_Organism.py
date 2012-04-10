#!/usr/local/bin/python

'''
#
# HMD_Organism.py
#
# History:
#
# 10/26/2006	lec
#	- TR 7492
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
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db


CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
REPORTPREFIX = 'HMD_'

reportLegend = 'Data Attributes:  M - MGI curated, C - HomoloGene calculated, B - MGI curated and HomoloGene calculated'
mgiCurated = 'M'
hgCurated = 'C'
bothCurated = 'B'

curated = []
calculated = []
otherEG = {}
mouseEG = {}
mouseMGI = {}

organismLookup = {10:'Chimpanzee', 13:'Dog'}

#
# NOTE: The following globals and getSortableOffset are variations on
# the original algorithm which can be found in the WI homology_report.cgi.py
# python module.
#

# special global variables used by getSortableOffset().  We declare them
# outside the function to yield a small speed benefit (we only create them
# once).

offset_cre = re.compile ('([pq])([0-9\.]+)')    # eg- q23.2

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
        #       which will sort as desired.  (using other as a primary
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

def runQueries(organismKey):

	global curated, calculated, otherEG, mouseEG, mouseMGI

	cmd = 'select distinct h1._Marker_key ' + \
		'from MRK_Homology_Cache h1, MRK_Homology_Cache h2 ' + \
		'where h1._Refs_key = 91485 ' + \
		'and h1._Organism_key = %s ' % (organismKey) + \
		'and h1._Class_key = h2._Class_key ' + \
		'and h2._Organism_key = 1'

	results = db.sql(cmd, 'auto')
	for r in results:
		calculated.append(r['_Marker_key'])

	cmd = 'select distinct h1._Marker_key ' + \
		'from MRK_Homology_Cache h1, MRK_Homology_Cache h2 ' + \
		'where h1._Refs_key != 91485 ' + \
		'and h1._Organism_key = %s ' % (organismKey) + \
		'and h1._Class_key = h2._Class_key ' + \
		'and h2._Organism_key = 1'

	results = db.sql(cmd, 'auto')
	for r in results:
    		curated.append(r['_Marker_key'])

	db.sql('select distinct otherMarkerKey = h1._Marker_key, mouseMarkerKey = h2._Marker_key ' + \
	       'into #allhomologies ' + \
	       'from MRK_Homology_Cache h1, MRK_Homology_Cache h2 ' + \
	       'where h1._Organism_key = %s ' % (organismKey) + \
	       'and h1._Class_key = h2._Class_key ' + \
	       'and h2._Organism_key = 1', None)

	db.sql('create index idx_hkey on #allhomologies(otherMarkerKey)', None)
	db.sql('create index idx_mkey on #allhomologies(mouseMarkerKey)', None)

	db.sql('''
		select h.otherMarkerKey, h.mouseMarkerKey, 
		m1._Organism_key as otherOrganism, 
		m1.symbol as otherSymbol, 
		m1.chromosome || m1.cytogeneticOffset as otherChr, 
		m1.chromosome, 
		m1.cytogeneticOffset, 
		m2._Organism_key as mouseOrganism, 
		m2.symbol as mouseSymbol, 
		m2.chromosome as mouseChr, 
		m2.cytogeneticOffset as mouseBand, 
		substring(m2.name, 1, 75) as mouseName, 
		o.offset as mouseOffset, 
        	case 
        	when o.offset >= 0 then str(o.offset, 10, 2) 
        	when o.offset = -999.0 then "       N/A" 
        	when o.offset = -1.0 then "  syntenic" 
        	end as mouseCm 
		into #homologies 
		from #allhomologies h, MRK_Marker m1, MRK_Marker m2, MRK_Offset o 
		where h.otherMarkerKey = m1._Marker_key 
		and h.mouseMarkerKey = m2._Marker_key 
		and h.mouseMarkerKey = o._Marker_key 
		and o.source = 0 
		''', None)

	db.sql('create index idx_hkey1 on #homologies(otherOrganism)', None)
	db.sql('create index idx_mkey1 on #homologies(mouseOrganism)', None)
	db.sql('create index idx_hkey2 on #homologies(otherSymbol)', None)
	db.sql('create index idx_mkey2 on #homologies(mouseSymbol)', None)

	# other entrezgene ids

	results = db.sql('select h.otherMarkerKey, a.accID from #homologies h, ACC_Accession a ' + \
		'where h.otherMarkerKey = a._Object_key ' + \
		'and a._MGIType_key = 2 ' + \
		'and a._LogicalDB_key = 55 ', 'auto')
	for r in results:
		otherEG[r['otherMarkerKey']] = r['accID']

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

	# sorted by other chromosome

	results = db.sql('select h.*, c.sequenceNum ' + \
		'from #homologies h, MRK_Chromosome c ' + \
		'where h.otherOrganism = c._Organism_key ' + \
		'and h.chromosome = c.chromosome ', 'auto')

	db.sql('drop table allhomologies', None)
	db.sql('drop table homologies', None)

	return results

def printDataAttributes(fp, key):

	if key in curated and not key in calculated:
  		fp.write(mgiCurated)
	elif key not in curated and key in calculated:
  		fp.write(hgCurated)
	else:
  		fp.write(bothCurated)

def process(organismKey, results):

	organism = organismLookup[organismKey]

	reportTitle = 'Orthology - %s vs. Mouse (Sorted by %s Chromosome)' % (organism, organism)
	reportName = REPORTPREFIX + organism
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'], printHeading = "MGI")
	fp.write(reportLegend + CRT + CRT)

	fp.write('%s Chr' % (organism) + TAB)
	fp.write('%s EntrezGene ID' % (organism) + TAB)
	fp.write('%s Symbol' % (organism) + TAB)
	fp.write('Mouse MGI Acc ID' + TAB)
	fp.write('Mouse Chr' + TAB)
	fp.write('Mouse cM' + TAB)
	fp.write('Mouse EntrezGene ID' + TAB)
	fp.write('Mouse Symbol' + TAB)
	fp.write('Mouse Name' + TAB)
	fp.write('Data Attributes' + CRT * 2)

	#
	# initialize a list to sort the other chromosome & offset values.
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
		sortKeys[2] = r['otherSymbol']
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
		fp.write(str(r['otherChr']) + TAB)

		if otherEG.has_key(r['otherMarkerKey']):
			fp.write(mgi_utils.prvalue(otherEG[r['otherMarkerKey']]))
		fp.write(TAB)
		fp.write(r['otherSymbol'] + TAB)
		fp.write(mouseMGI[r['mouseMarkerKey']] + TAB)
		fp.write(r['mouseChr'] + TAB)
		fp.write(r['mouseCm'] + TAB)

		if mouseEG.has_key(r['mouseMarkerKey']):
			fp.write(mgi_utils.prvalue(mouseEG[r['mouseMarkerKey']]))
		fp.write(TAB)

		fp.write(r['mouseSymbol'] + TAB)
		fp.write(r['mouseName'] + TAB)
		printDataAttributes(fp, r['otherMarkerKey'])
		fp.write(CRT)

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.finish_nonps(fp)

#
# Main
#

for organismKey in organismLookup.keys():
    db.useOneConnection(1)
    results = runQueries(organismKey)
    process(organismKey, results)

db.useOneConnection(0)

