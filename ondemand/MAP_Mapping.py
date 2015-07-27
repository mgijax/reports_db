#!/usr/local/bin/python

'''
#
# MAP_Mapping.py 04/20/2010
#
# Report:
#       Tab-delimited reports for all of the mapping data in the database.
#
# Usage:
#       MAP_Mapping.py
#
# Used by:
#
# Output format:
#
#
# History:
#
# 04/25/2012
#	- TR11035/postgres cleanup/PI cleanup
#	- move to 'ondemand' directory
#
# 12/28/2011	lec
#	- changed non-ansi-standard query to left outer join
#
# mhall 04/20/2010
# Initial version
#
'''

import sys
import os
import string
import mgi_utils
import reportlib
import db

db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE(False)

TAB = reportlib.TAB
CRT = reportlib.CRT

db.useOneConnection(1)

allPanels = ['1839', '8', '7', '2944', '3241', '4867', '4868', '4377', '4378', '2941', '3383', '4347', '6', '4869']
allChromosomes = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X','Y','XY']

typings = {}
typings2 = {}
cross = {}
strain = {}
lastMod = {}
abbrevHT = {}
abbrevHO = {}
strainHT = {}
strainHO = {}
header = {}
header2 = {}
refs = {}
maxLength = 0
preheader = {}
preheader2 = {}
mergedID = {}

results = db.sql('''
	select _Cross_key, rowNumber, colNumber, data 
	from CRS_Typings
	order by _Cross_key, rowNumber, colNumber
	''', 'auto')
for item in results:
	key = 'cross_key' + str(item['_Cross_key']) + 'rowNumber' + str(item['rowNumber'])
	itemList = item['data'].rstrip().split(' ')
	
	for item in itemList:
		if key in typings:
			if item != '':
				typings[key] = typings[key] + TAB + item 
				typings2[key] = typings2[key] + str(item).ljust(6)
		else:
			typings[key] = item 
			typings2[key] = str(item).ljust(6)

for key in typings:
	temp = key.split('rowNumber')
	temp2 = temp[0].split('cross_key')
	crsKey = temp2[1]
	temp3 = typings[key].split('	')
	
	preheader[crsKey] = 'CHR	MGI ID	SYMBOL	ANIMAL#' + TAB*(len(temp3)) + 'J:NUMBER'
	
	# Increadibly finicky formatting.  There has to be a better way to do this, but I've come
	# too far now
	
	preheader2[crsKey] = \
		'CHR'.ljust(8) + \
		'MGI ID'.ljust(12) + \
		'SYMBOL'.ljust(20) + \
		'ANIMAL#'.ljust(11) + \
		' '*6*(len(temp3)-2) + \
		' J:NUMBER'
	preheader2[crsKey] = \
		preheader2[crsKey] + '\n' + \
		'-'*7 + ' ' + \
		'-'*11 + ' ' + \
		'-'*19 + ' '*6*len(temp3) + \
		' ' + '-'*16
	
results = db.sql('''
	select _Cross_key, whoseCross, femaleStrain, maleStrain,
		to_char(g.modification_date, 'YYYYMMDD') as mDate,
		abbrevHT, abbrevHO, strainHT, strainHO 
	from crs_cross_view
	where _Cross_key in (1839, 8, 7, 2944, 3241, 4867, 4868, 4377, 4378, 2941, 3383, 4347, 6, 4869)
	''', 'auto')
for item in results:
	key = str(item['_Cross_key'])
	cross[key] = item['whoseCross']
	
	if item['_Cross_key'] != 4347:
		strain[key] = item['femaleStrain'] + ' x ' + item['maleStrain']
	else:
		strain[key] = item['femaleStrain']
		
	lastMod[key] = str(item['modification_date'])
	abbrevHT[key] = item['abbrevHT']
	abbrevHO[key] = item['abbrevHO']
	strainHT[key] = item['strainHT']
	strainHO[key] = item['strainHO']

results = db.sql('''
	select _Cross_key, sequenceNum, name 
	from CRS_Progeny
	order by _Cross_key, sequenceNum
	''', 'auto')
for item in results:
	key = 'cross_key' + str(item['_Cross_key'])
	if key in header:
		header[key] = header[key] + item['name'] + TAB
		header2[key] = header2[key] + str(item['name']).ljust(6)
	else:
		header[key] = TAB + TAB + TAB + item['name'] + TAB
		header2[key] = ' '*40 + str(item['name']) + ' '*(6-len(item['name']))	

results = db.sql('''
	select bv.jnumID, cr._Cross_key, cr._Marker_key 
	from CRS_References cr, BIB_View bv
	where cr._Refs_key != null
	and cr._Refs_key = bv._Refs_key
	order by cr._Cross_key, cr._Marker_key, bv.jnumID
	''', 'auto')
for item in results:
	key = 'cross_key' + str(item['_Cross_key']) + 'marker_key' + str(item['_Marker_key'])
	if key in refs: 
		refs[key] = refs[key] + ', ' + item['jnumID']
	else:
		refs[key] = item['jnumID']


db.sql('''
	select distinct cm._Marker_key, a.accID 
	into #tmp_accID
	from CRS_Matrix cm 
		LEFT OUTER JOIN ACC_Accession a on (
			cm._Marker_key = a._Object_key 
			and a._MGIType_key = 2 
			and a.prefixPart = 'MGI:'
			and a.private = 0 and a.preferred = 1)
	''', None)

results = db.sql('''
	select mc._Marker_key, a.accID
	from MRK_Current mc, ACC_Accession a
	where mc._Marker_key in (select _Marker_key from #tmp_accID where accID = null)
	and mc._Current_key = a._Object_key 
	and a.prefixPart = 'MGI:'
	and a.private = 0 
	and a.preferred = 1 
	and a._MGIType_key = 2
	''', 'auto')
for item in results:
	mergedID[item['_Marker_key']] = item['accID']

cmd = '''
	select distinct cm.chromosome, cm._Marker_key, cm._Cross_key, cm.rowNumber, 
		mm.symbol, cm.otherSymbol, a.accID, mm._Marker_Status_key
	from CRS_Matrix cm
	     LEFT OUTER JOIN MRK_Marker mm on (cm._Marker_key = mm._Marker_key)
	     LEFT OUTER JOIN ACC_Accession a on (cm._Marker_key = a._Object_key
	         and a.prefixPart = 'MGI:' 
		 and a.private = 0 
		 and a.preferred = 1 and a._MGIType_key = 2)
	     LEFT OUTER JOIN CRS_References cr on (cm._Cross_key = cr._Cross_key
	         and cm._Marker_key = cr._Marker_key)
	where cm._Cross_key = %s
	and cm.chromosome = '%s'
	order by rowNumber
	'''

for panel in allPanels:

	fp = reportlib.init('MGI_' + \
		cross[panel].replace(' ', '_'	).replace('(', '').replace(')', '') + \
		'_Panel', fileExt = '.rpt', outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

	fp.write('Mapping Panel: ' + cross[panel] + '\n')

	fp.write('Cross designation: ' + strain[panel] + '\n')

	fp.write('Last Modified: ' + lastMod[panel] + '\n')

	fp.write('WARNING:  The MGI IDs and/or SYMBOL may be outdated\n')

	if panel != '3383':
		fp.write('Legend:\n')
		fp.write('    ' + abbrevHT[panel] + ' indicates allele from ' + strainHT[panel] + '\n')
		fp.write('    ' + abbrevHO[panel] + ' indicates allele from ' + strainHO[panel] + '\n')
		fp.write('    . indicates animal was not typed.\n\n')	
	else:
		fp.write('\n\n')

	fp.write(preheader[panel] + '\n')
	fp.write(header['cross_key' +panel] + '\n')

	for chr in allChromosomes:
		#print 'Panel: ' + panel + ' Chr: ' + chr
		results = None
		results = db.sql(cmd % (panel, chr), 'auto')
		
		for item in results:
			key = 'cross_key' + str(item['_Cross_key']) + 'rowNumber' + str(item['rowNumber'])
			#print item
			line = item['chromosome']
			line2 = item['chromosome'].ljust(8)

			# Calculate the accID String

			accID = item['accID']
			marker_status = item['_Marker_Status_key']
			if accID != None:
				if marker_status == 2:
					line = line + TAB + 'withdrawn'
					line2 = line2 + 'withdrawn'.ljust(12)														
				else:
					line = line + TAB + accID
					line2 = line2 + accID.ljust(12)
			else:
				if marker_status == 2 and mergedID.has_key(item['_Marker_key']):
					line = line + TAB + mergedID[item['_Marker_key']]
					line2 = line2 + mergedID[item['_Marker_key']].ljust(12)					
				else:							
					line = line + TAB + 'n/a'
					line2 = line2 + 'n/a'.ljust(12)
			
			# Calculate the symbol, either the mgi one, or the data supplied one
			
			symbol = item['symbol']
			if symbol != None:
				line = line + TAB + symbol
				line2 = line2 + symbol.ljust(20)
			else:
				line = line + TAB + item['otherSymbol']
				line2 = line2 + item['otherSymbol'].ljust(20)
			
			# Add in the typings
			
			line = line + TAB + typings[key]
			line2 = line2 + typings2[key]
			
			# Add in the optional JNUM Column
			
			refKey = 'cross_key' + str(item['_Cross_key']) + 'marker_key' + str(item['_Marker_key'])
			if refKey in refs:
				line = line + TAB + refs[refKey]
				line2 = line2 + refs[refKey].ljust(16)
			else:
				line = line + TAB + ''
				line2 = line2 + ''.ljust(16)			
			
			# Finally print the line
			
			fp.write(line + '\n')
			
	fp.close()
		
#for item in typings:
#	print 'item: ' + item
#	print 'data: ' + typings[item]
		
reportlib.finish_nonps(fp)
db.useOneConnection(0)
