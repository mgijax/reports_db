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

currentReport = 'Nomenclature-current.html'

# remove current report link if it exists
if os.path.isfile('%s/%s' % (os.environ['FTPREPORTDIR2'], currentReport)):
	os.remove('%s/%s' % (os.environ['FTPREPORTDIR2'], currentReport))

# move existing Nomen reports to the archive
os.system('mv %s/Nomenclature-*.html %s/archive/nomen' % (os.environ['FTPREPORTDIR2'], os.environ['FTPREPORTDIR2']))

if len(sys.argv) > 1:
	reportName = 'Nomenclature-' + sys.argv[1]
	currentDate = sys.argv[2]
else:
	reportName = 'Nomenclature-' + mgi_utils.date('%Y-%m-%d')
	currentDate = mgi_utils.date('%m/%d/%Y')

results = db.sql('select convert(varchar(25), dateadd(day, -7, "%s"))' % (currentDate), 'auto')
bdate = results[0]['']

results = db.sql('select convert(varchar(25), dateadd(day, 0, "%s"))' % (currentDate), 'auto')
edate = results[0]['']

title = 'Updates to Mouse Nomenclature from %s to %s' % (bdate, edate)
fp = reportlib.init(reportName, title, os.environ['FTPREPORTDIR2'], isHTML = 1)

fp.write('J:23000 generally indicates gene family nomenclature revision event.\n\n')

cmd = []

cmd.append('select h._Marker_key, m.symbol, name = substring(m.name,1,35), ' + \
'chromosome = substring(m.chromosome,1,2), ' + \
'c.sequenceNum, ' + \
'b.jnumID, author = substring(b._primary,1,16) ' + \
'into #markers ' + \
'from MRK_Marker m, MRK_History h, MRK_Chromosome c, BIB_All_View b ' + \
'where m._Species_key = 1 ' + \
'and m._Marker_key = h._History_key ' + \
'and h.event_date between dateadd(day, -7, "%s") ' % (currentDate) + \
'and dateadd(day, 1, "%s") ' % (currentDate) + \
'and m.chromosome = c.chromosome ' + \
'and m._Species_key = c._Species_key ' + \
'and h._Refs_key = b._Refs_key')

cmd.append('create nonclustered index idx_key on #markers(_Marker_key)')

# get primary MGI ids (2)
cmd.append('select distinct m._Marker_key, a.accID ' + \
'from #markers m, MRK_Acc_View a ' + \
'where m._Marker_key = a._Object_key ' + \
'and a.prefixPart = "MGI:" ' + \
'and a._LogicalDB_key = 1 ' + \
'and a.preferred = 1')

# other MGI ids (3)
cmd.append('select distinct m._Marker_key, a.accID ' + \
'from #markers m, MRK_Acc_View a ' + \
'where m._Marker_key = a._Object_key ' + \
'and a.prefixPart = "MGI:" ' + \
'and a._LogicalDB_key = 1 ' + \
'and a.preferred = 0')

# get sequence ids (4)
cmd.append('select distinct m._Marker_key, a.accID ' + \
'from #markers m, MRK_Acc_View a ' + \
'where m._Marker_key = a._Object_key ' + \
'and a._LogicalDB_key = 9')

# get human ortholog (5)
cmd.append('select distinct m._Marker_key, n.humanSymbol ' + \
'from #markers m, MRK_Acc_View ma, NOM_Acc_View na, NOM_Marker n ' + \
'where m._Marker_key = ma._Object_key ' + \
'and ma.prefixPart = "MGI:" ' + \
'and ma._LogicalDB_key = 1 ' + \
'and ma.preferred = 1 ' + \
'and ma.accID = na.accID ' + \
'and na._Object_key = n._Nomen_key ' + \
'and n.humanSymbol is not null')

# retrieve markers, sort (6)
cmd.append('select distinct * from #markers order by sequenceNum, symbol')

results = db.sql(cmd, 'auto')

fp.write('%-2s %-25s %-35s %-10s %-20s %-25s %-75s %-15s %-25s\n' % ('Ch', 'Symbol', 'Gene Name', 'J#', 'First Author    ', 'MGI ID', 'Sequence ID', 'Human Ortholog', 'Other MGI IDs'))
fp.write('%-2s %-25s %-35s %-10s %-20s %-25s %-75s %-15s %-25s\n' % ('--', '------', '---------', '--', '----------------', '------', '-----------', '--------------', '-------------'))

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

	fp.write('%-2s ' % (r['chromosome']))

	fp.write('%s%-25s%s ' % (reportlib.create_accession_anchor(primaryID[key]), symbol, reportlib.close_accession_anchor()))
	fp.write('%-35s ' % (r['name']))
	fp.write('%s%-10s%s ' % (reportlib.create_accession_anchor(r['jnumID']), r['jnumID'], reportlib.close_accession_anchor()))
	fp.write('%-20s ' % (r['author']))
	fp.write('%s%-25s%s ' % (reportlib.create_accession_anchor(primaryID[key]), primaryID[key], reportlib.close_accession_anchor()))

	if seqIDs.has_key(key):
		fp.write('%-75s ' % (string.join(seqIDs[key], ',')))
	else:
		fp.write('%-75s ' % (''))

	if human.has_key(key):
		fp.write('%-15s ' % (human[key]))
	else:
		fp.write('%-15s ' % (''))

	if otherIDs.has_key(key):
		fp.write(string.join(otherIDs[key], ','))
	fp.write(reportlib.CRT)

	rows = rows + 1

fp.write(reportlib.CRT + '(%d rows affected)' % (rows) + reportlib.CRT)
reportlib.trailer(fp)
reportlib.finish_nonps(fp, isHTML = 1)	# non-postscript file

# re-create a symbolic link between the new file and the current file

os.chdir(os.environ['FTPREPORTDIR2'])
os.symlink(reportName + '.html', currentReport)

