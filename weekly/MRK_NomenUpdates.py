#!/usr/local/bin/python

'''
#
# TR 830 - Recent Nomenclature Events Report in HTML format
#
# Usage:
#       MRK_NomenUpdates.py
#
# Notes:
#	- all reports use db default of public login
#	- all reports use server/database default of environment
#	- use lowercase for all SQL commands (i.e. select not SELECT)
#	- all public SQL reports require the header and footer
#	- all private SQL reports require the header
#
# History:
#
# lec	10/20/2000
#	- TR 2023; display J: as "J:#####"
#
# lec	06/08/2000
#	- TR 1650; added Sequence ID
#
# lec	01/04/99
#	- created
#
'''
 
import sys 
import os
import string
import db
import mgi_utils
import reportlib

#
# Main
#

currentReport = 'Nomenclature-current.html'

# remove current report link if it exists
if os.path.isfile('%s/%s' % (os.environ['FTPREPORTDIR'], currentReport)):
	os.remove('%s/%s' % (os.environ['FTPREPORTDIR'], currentReport))

# move existing Nomen reports to the archive
os.system('mv %s/Nomenclature-*.html %s/archive/nomen' % (os.environ['FTPREPORTDIR'], os.environ['FTPREPORTDIR']))

reportName = 'Nomenclature-' + mgi_utils.date('%Y-%m-%d')

currentDate = mgi_utils.date('%m/%d/%Y')

results = db.sql('select convert(varchar(25), dateadd(day, -7, "%s"))' % (currentDate), 'auto')
bdate = results[0]['']

results = db.sql('select convert(varchar(25), dateadd(day, 0, "%s"))' % (currentDate), 'auto')
edate = results[0]['']

title = 'Updates to Mouse Nomenclature from %s to %s' % (bdate, edate)
fp = reportlib.init(reportName, title, os.environ['FTPREPORTDIR'], isHTML = 1)

fp.write('J:23000 generally indicates gene family nomenclature revision event.\n\n')

cmd = []

cmd.append('select m._Marker_key, m.mgiID, c.sequenceNum, h._Refs_key ' + \
'into #m1 ' + \
'from MRK_Mouse_View m, MRK_History h, MRK_Chromosome c ' + \
'where m.creation_date between dateadd(day, -7, "%s") ' % (currentDate) + \
'and dateadd(day, 1, "%s") ' % (currentDate) + \
'and m._Marker_key = h._Marker_key ' + \
'and m._Marker_key = h._History_key ' + \
'and h._Marker_Event_key = 1 ' + \
'and m.chromosome = c.chromosome ' + \
'and m._Species_key = c._Species_key ' + \
'union ' + \
'select h._History_key, m.mgiID, sequenceNum = 100, h._Refs_key ' + \
'from MRK_History h, MRK_Mouse_View m ' + \
'where h.event_date between dateadd(day, -7, "%s") ' % (currentDate) + \
'and dateadd(day, 1, "%s") ' % (currentDate) + \
'and h._Marker_key = m._Marker_key ' + \
'and h._Marker_Event_key in (2,3,4,5)'
)

cmd.append('select m._Marker_key, ' + \
'chr = substring(r.chromosome,1,2), ' + \
'm.mgiID, r.symbol,  ' + \
'name = substring(r.name,1,25),  ' + \
'jnumID, ' + \
'author = substring(b._primary, 1, 16) ' + \
'into #m2 ' + \
'from #m1 m, MRK_Marker r, BIB_All_View b ' + \
'where m._Marker_key = r._Marker_key ' + \
'and m._Refs_key = b._Refs_key ' + \
'order by m.sequenceNum, r.symbol'
)

cmd.append('select m.*, a.accID ' + \
'from #m2 m, MRK_Acc_View a ' + \
'where m._Marker_key = a._Object_key ' + \
'and a._LogicalDB_key = 9 '
)

results = db.sql(cmd, 'auto')

fp.write('%-2s %-25s %-25s %-10s %-20s %-25s\n' % ('Ch', 'Symbol', 'Gene Name', 'J#', 'First Author    ', 'Sequence ID'))
fp.write('%-2s %-25s %-25s %-10s %-20s %-25s\n' % ('--', '------', '---------', '--', '----------------', '-----------'))

rows = 0
prevMarker = ''
sequence = []

for r in results[2]:

	if prevMarker != r['_Marker_key']:

		if len(sequence) > 0:
			fp.write(string.join(sequence, ','))
		sequence = []

		if prevMarker != '':
			fp.write(reportlib.CRT)

		fp.write('%-2s ' % (r['chr']))

		if r['mgiID'] != "None":
			fp.write('%s%-25s%s ' % (reportlib.create_accession_anchor(r['mgiID']), r['symbol'], reportlib.close_accession_anchor()))
		else:
			fp.write('%-25s ' % (r['symbol']))
		
		fp.write('%-25s ' % (r['name']))
		fp.write('%s%-10s%s ' % (reportlib.create_accession_anchor(r['jnumID']), r['jnumID'], reportlib.close_accession_anchor()))
		fp.write('%-20s ' % (r['author']))

		prevMarker = r['_Marker_key']
		rows = rows + 1

	sequence.append(r['accID'])

fp.write(string.join(sequence, ',') + reportlib.CRT)
fp.write(reportlib.CRT + '(%d rows affected)' % (rows) + reportlib.CRT)

reportlib.trailer(fp)
reportlib.finish_nonps(fp, isHTML = 1)	# non-postscript file

# re-create a symbolic link between the new file and the current file

os.chdir(os.environ['FTPREPORTDIR'])
os.symlink(reportName + '.html', currentReport)

