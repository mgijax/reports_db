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
#       HMD_Human.py
#
# History:
#
# lec	01/04/2004
#	- TR 5939; LocusLink->EntrezGene
#
# lec	07/23/2004
#	- TR 5611; added humanCurated/humanCalculated
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
REPORTNAME = 'HMD_Human'

reportLegend = 'Data Attributes:  M - MGI humanCurated, C - HomoloGene humanCalculated, B - MGI humanCurated and HomoloGene humanCalculated'
mgiCurated = 'M'
hgCurated = 'C'
bothCurated = 'B'

humanCurated = []
humanCalculated = []
humanEG = {}
humanHGNC = {}
mouseEG = {}
mouseMGI = {}
mouseCoords = {}

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

	global humanCurated, humanCalculated, humanEG, humanHGNC, mouseEG, mouseMGI
	global mouseCoords

	#
	# human orthologs
	#
	cmd = '''
		select distinct h1._Marker_key 
		from MRK_Homology_Cache h1, MRK_Homology_Cache h2 
		where h1._Refs_key = 91485 
		and h1._Organism_key = 2 
		and h1._Class_key = h2._Class_key 
		and h2._Organism_key = 1
		'''
	results = db.sql(cmd, 'auto')
	for r in results:
		humanCalculated.append(r['_Marker_key'])

	cmd = '''
		select distinct h1._Marker_key 
		from MRK_Homology_Cache h1, MRK_Homology_Cache h2 
		where h1._Refs_key != 91485 
		and h1._Organism_key = 2 
		and h1._Class_key = h2._Class_key 
		and h2._Organism_key = 1
		'''
	results = db.sql(cmd, 'auto')
	for r in results:
    		humanCurated.append(r['_Marker_key'])

	#
	# human/mouse
	#
	db.sql('''
	       select distinct h1._Marker_key as humanMarkerKey, h2._Marker_key as mouseMarkerKey
	       into #humanHomologies 
	       from MRK_Homology_Cache h1, MRK_Homology_Cache h2 
	       where h1._Organism_key = 2 
	       and h1._Class_key = h2._Class_key 
	       and h2._Organism_key = 1
	       ''', None)

	db.sql('create index idx_hkey on #humanHomologies(humanMarkerKey)', None)
	db.sql('create index idx_mkey on #humanHomologies(mouseMarkerKey)', None)

	db.sql('''
		select h.humanMarkerKey, h.mouseMarkerKey, 
		       m1._Organism_key as humanOrganism, 
		       m1.symbol as humanSymbol, 
		       m1.chromosome || m1.cytogeneticOffset as humanChr, 
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
		from #humanHomologies h, MRK_Marker m1, MRK_Marker m2, MRK_Offset o 
		where h.humanMarkerKey = m1._Marker_key 
		      and h.mouseMarkerKey = m2._Marker_key 
		      and h.mouseMarkerKey = o._Marker_key 
		      and o.source = 0 
		''', None)

	db.sql('create index idx_hkey1 on #homologies(humanOrganism)', None)
	db.sql('create index idx_mkey1 on #homologies(mouseOrganism)', None)
	db.sql('create index idx_hkey2 on #homologies(humanSymbol)', None)
	db.sql('create index idx_mkey2 on #homologies(mouseSymbol)', None)

	##

	# human entrezgene ids

	results = db.sql('''
		select h.humanMarkerKey, a.accID 
		from #homologies h, ACC_Accession a 
		where h.humanMarkerKey = a._Object_key 
		and a._MGIType_key = 2 
		and a._LogicalDB_key = 55 
		''', 'auto')
	for r in results:
		humanEG[r['humanMarkerKey']] = r['accID']

	# human hgnc ids

	results = db.sql('''
		select h.humanMarkerKey, a.accID 
		from #homologies h, ACC_Accession a 
		where h.humanMarkerKey = a._Object_key 
		and a._MGIType_key = 2 
		and a._LogicalDB_key = 64 
		''', 'auto')
	for r in results:
		humanHGNC[r['humanMarkerKey']] = r['accID']

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
		and a.prefixPart = "MGI:" 
		and a._LogicalDB_key = 1 
		and a.preferred = 1 
		''', 'auto')
	for r in results:
		mouseMGI[r['mouseMarkerKey']] = r['accID']

	#
	# coordinates
	#
	results = db.sql('''	
    		select h.mousemarkerkey,
           	       c.strand, 
	   	       convert(int, c.startCoordinate) as startC,
	   	       convert(int, c.endCoordinate) as endC
    		from #homologies h, MRK_Location_Cache c
    		where h.mousemarkerkey = c._marker_key
		''', 'auto')
	mouseCoords = {}
	for r in results:
    		key = r['mousemarkerkey']
    		value = r
    		if not mouseCoords.has_key(key):
			mouseCoords[key] = []
    		mouseCoords[key].append(value)

	##

	# sorted by human chromosome

	results = db.sql('''
		select h.*, c.sequenceNum
		from #homologies h, MRK_Chromosome c 
		where h.humanOrganism = c._Organism_key 
		and h.chromosome = c.chromosome 
		''', 'auto')

	return results

def printDataAttributes(fp, key):

	if key in humanCurated and not key in humanCalculated:
  		fp.write(mgiCurated)
	elif key not in humanCurated and key in humanCalculated:
  		fp.write(hgCurated)
	else:
  		fp.write(bothCurated)

def processSort(results):

	reportTitle = 'Orthology - Human vs. Mouse (Sorted by Human Chromosome)'
	reportName = REPORTNAME + '5'
	
	fp = reportlib.init(reportName, reportTitle, os.environ['REPORTOUTPUTDIR'], printHeading = "MGI")
	fp.write(reportLegend + CRT + CRT)

	fp.write('Human EntrezGene ID' + TAB)
	fp.write('Human Symbol' + TAB)
	fp.write('Human Chr' + TAB)
	fp.write('Human HGNC ID' + TAB)
	fp.write('Mouse MGI Acc ID' + TAB)
	fp.write('Mouse EntrezGene ID' + TAB)
	fp.write('Mouse Symbol' + TAB)
	fp.write('Mouse Name' + TAB)
	fp.write('Mouse Chr' + TAB)
	fp.write('Mouse cM' + TAB)
	fp.write('Mouse Genome Coordinate Start' + TAB)
	fp.write('Mouse Genome Coordinate End' + TAB)
	fp.write('Mouse Strand' + TAB)
	fp.write('Mouse Cytogenetic Band' + TAB)
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

		if humanEG.has_key(r['humanMarkerKey']):
			fp.write(mgi_utils.prvalue(humanEG[r['humanMarkerKey']]))
		fp.write(TAB)
		fp.write(r['humanSymbol'] + TAB)
		fp.write(mgi_utils.prvalue(r['humanChr']) + TAB)

		if humanHGNC.has_key(r['humanMarkerKey']):
			fp.write(mgi_utils.prvalue(humanHGNC[r['humanMarkerKey']]))
		fp.write(TAB)

		# add rat EG
		# add rat synbol
		# add rat chr

		fp.write(mouseMGI[r['mouseMarkerKey']] + TAB)
		if mouseEG.has_key(r['mouseMarkerKey']):
			fp.write(mgi_utils.prvalue(mouseEG[r['mouseMarkerKey']]))
		fp.write(TAB)
		fp.write(r['mouseSymbol'] + TAB)
		fp.write(r['mouseName'] + TAB)
		fp.write(r['mouseChr'] + TAB)
		fp.write(r['mouseCm'] + TAB)

    		if mouseCoords.has_key(r['mouseMarkerKey']):
        		fp.write(mgi_utils.prvalue(mouseCoords[r['mouseMarkerKey']][0]['startC']) + TAB)
        		fp.write(mgi_utils.prvalue(mouseCoords[r['mouseMarkerKey']][0]['endC']) + TAB)
        		fp.write(mgi_utils.prvalue(mouseCoords[r['mouseMarkerKey']][0]['strand']) + TAB)
    		else:
        		fp.write(TAB + TAB + TAB)

		fp.write(mgi_utils.prvalue(r['mouseBand']) + TAB)
		
		printDataAttributes(fp, r['humanMarkerKey'])
		fp.write(CRT)

	fp.write(CRT + '(%d rows affected)' % (count) + CRT)
	reportlib.finish_nonps(fp)

#
# Main
#

db.useOneConnection(1)
processSort(runQueries())
db.useOneConnection(0)

