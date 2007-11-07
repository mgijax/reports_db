#!/usr/local/bin/python

'''
#
# TR 830 - Recent Nomenclature Events Report in HTML format
#
# Usage:
#       MRK_NomenUpdates.py
#
# History:
#
# lec	06/25/2002
#	- TR 3827; add Human Ortholog 
#
# lec	06/14/2000
#	- TR 2613; display current MGI IDs and Other MGI IDs
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
import mgi_html
import mgi_utils
import reportlib

#
# Main
#

reportDir = os.environ['FTPREPORTDIR']
reportName = 'Nomenclature-' + mgi_utils.date('%Y-%m-%d')
currentReport = 'Nomenclature-current'
currentDate = mgi_utils.date('%m/%d/%Y')

# remove current report links if they exists
if os.path.isfile('%s/%s' % (reportDir, currentReport + '.html')):
	os.remove('%s/%s' % (reportDir, currentReport + '.html'))
if os.path.isfile('%s/%s' % (reportDir, currentReport + '.rpt')):
	os.remove('%s/%s' % (reportDir, currentReport + '.rpt'))

# move existing Nomen reports to the archive
os.system('mv %s/Nomenclature-*.html %s/archive/nomen' % (reportDir, reportDir))
os.system('mv %s/Nomenclature-*.rpt %s/archive/nomen' % (reportDir, reportDir))

results = db.sql('select convert(varchar(25), dateadd(day, -7, "%s"))' % (currentDate), 'auto')
bdate = results[0]['']

results = db.sql('select convert(varchar(25), dateadd(day, 0, "%s"))' % (currentDate), 'auto')
edate = results[0]['']

title = 'Updates to Mouse Nomenclature from %s to %s' % (bdate, edate)
fpHTML = reportlib.init(reportName, title, os.environ['FTPREPORTDIR'], isHTML = 1, printHeading = "MGI")
fpRpt = reportlib.init(reportName, outputdir = os.environ['FTPREPORTDIR'], printHeading = None)

fpHTML.write('J:23000 generally indicates gene family nomenclature revision event.\n\n')

cmd = []

cmd.append('select h._Marker_key, m.symbol, name = substring(m.name,1,35), ' + \
'chromosome = substring(m.chromosome,1,2), ' + \
'c.sequenceNum, ' + \
'b.jnumID, author = substring(b._primary,1,16) ' + \
'into #markers ' + \
'from MRK_Marker m, MRK_History h, MRK_Chromosome c, BIB_All_View b ' + \
'where m._Organism_key = 1 ' + \
'and m._Marker_key = h._History_key ' + \
'and h.event_date between dateadd(day, -7, "%s") ' % (currentDate) + \
'and dateadd(day, 1, "%s") ' % (currentDate) + \
'and m.chromosome = c.chromosome ' + \
'and m._Organism_key = c._Organism_key ' + \
'and h._Refs_key = b._Refs_key')

cmd.append('create nonclustered index idx_key on #markers(_Marker_key)')

# get primary MGI ids (2)
cmd.append('select distinct m._Marker_key, a.accID ' + \
'from #markers m, ACC_Accession a ' + \
'where m._Marker_key = a._Object_key ' + \
'and a._MGIType_key = 2 ' + \
'and a.prefixPart = "MGI:" ' + \
'and a._LogicalDB_key = 1 ' + \
'and a.preferred = 1')

# other MGI ids (3)
cmd.append('select distinct m._Marker_key, a.accID ' + \
'from #markers m, ACC_Accession a ' + \
'where m._Marker_key = a._Object_key ' + \
'and a._MGIType_key = 2 ' + \
'and a.prefixPart = "MGI:" ' + \
'and a._LogicalDB_key = 1 ' + \
'and a.preferred = 0')

# get sequence ids (4)
cmd.append('select distinct m._Marker_key, a.accID ' + \
'from #markers m, ACC_Accession a ' + \
'where m._Marker_key = a._Object_key ' + \
'and a._MGIType_key = 2 ' + \
'and a._LogicalDB_key = 9')

# get human ortholog (5)
cmd.append('select distinct m._Marker_key, humanSymbol = h.symbol ' + \
'from #markers m, MRK_Marker h, MRK_Homology_Cache hm1, MRK_Homology_Cache hm2 ' + \
'where m._Marker_key = hm1._Marker_key ' + \
'and hm1._Homology_key = hm2._Homology_key ' + \
'and hm2._Marker_key = h._Marker_key ' + \
'and hm2._Organism_key = 1')

# retrieve markers, sort (6)
cmd.append('select distinct * from #markers order by sequenceNum, symbol')

results = db.sql(cmd, 'auto')

fpHTML.write('%-2s %-25s %-35s %-10s %-20s %-25s %-75s %-15s %-25s\n' % ('Ch', 'Symbol', 'Gene Name', 'J#', 'First Author    ', 'MGI ID', 'Sequence ID', 'Human Ortholog', 'Other MGI IDs'))
fpHTML.write('%-2s %-25s %-35s %-10s %-20s %-25s %-75s %-15s %-25s\n' % ('--', '------', '---------', '--', '----------------', '------', '-----------', '--------------', '-------------'))

primaryID = {}
for r in results[2]:
	primaryID[r['_Marker_key']] = r['accID']

otherIDs = {}
for r in results[3]:
	if not otherIDs.has_key(r['_Marker_key']):
		otherIDs[r['_Marker_key']] = []
	otherIDs[r['_Marker_key']].append(r['accID'])

seqIDs = {}
for r in results[4]:
	if not seqIDs.has_key(r['_Marker_key']):
		seqIDs[r['_Marker_key']] = []
	seqIDs[r['_Marker_key']].append(r['accID'])

human = {}
for r in results[5]:
	human[r['_Marker_key']] = r['humanSymbol']

rows = 0
for r in results[6]:

	key = r['_Marker_key']
	symbol = mgi_html.escape(r['symbol'])

	fpHTML.write('%-2s ' % (r['chromosome']))

	fpHTML.write('%s%-25s%s ' % (reportlib.create_accession_anchor(primaryID[key]), symbol, reportlib.close_accession_anchor()))

	fpHTML.write('%-35s ' % (r['name']))

	fpHTML.write('%s%-10s%s ' % (reportlib.create_accession_anchor(r['jnumID']), r['jnumID'], reportlib.close_accession_anchor()))

	fpHTML.write('%-20s ' % (r['author']))

	fpHTML.write('%s%-25s%s ' % (reportlib.create_accession_anchor(primaryID[key]), primaryID[key], reportlib.close_accession_anchor()))

	fpRpt.write('%s\t%s\t%s\t%s\t%s\t%s\t' % (r['chromosome'], symbol, r['name'], r['jnumID'], r['author'], primaryID[key]))

	if seqIDs.has_key(key):
		fpHTML.write('%-75s ' % (string.join(seqIDs[key], ',')))
		fpRpt.write('%s\t' % (string.join(seqIDs[key], ',')))
	else:
		fpHTML.write('%-75s ' % (''))
		fpRpt.write('\t')

	if human.has_key(key):
		fpHTML.write('%-15s ' % (human[key]))
		fpRpt.write('%s\t' % (human[key]))
	else:
		fpHTML.write('%-15s ' % (''))
		fpRpt.write('\t')

	if otherIDs.has_key(key):
		fpHTML.write(string.join(otherIDs[key], ','))
		fpRpt.write('%s' % (string.join(otherIDs[key], ',')))

	fpHTML.write(reportlib.CRT)
	fpRpt.write(reportlib.CRT)

	rows = rows + 1

fpHTML.write(reportlib.CRT + '(%d rows affected)' % (rows) + reportlib.CRT)
reportlib.finish_nonps(fpHTML, isHTML = 1)	# non-postscript file
reportlib.finish_nonps(fpRpt)			# non-postscript file

# re-create a symbolic link between the new file and the current file
os.chdir(reportDir)
os.symlink(reportName + '.html', currentReport + '.html')
os.symlink(reportName + '.rpt', currentReport + '.rpt')
