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
# lec	03/28/2012
#	- TR 11027; create_accession_anchor() for MGI 5.0
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
import mgi_html
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

db.useOneConnection(1)

#
# purposely did not add this to lib_py_postgres/pg_db.py/translate_be() (for now)
#

if os.environ['DB_TYPE'] == 'postgres':
	results = db.sql('select to_char(now() - interval \'7 days\', \'YYYY-MM-DD\')', 'auto')
	bdate = results[0]['to_char']

	results = db.sql('select to_char(now(), \'YYYY-MM-DD\')', 'auto')
	edate = results[0]['to_char']

	datequery = '\'%s\' and \'%s\'' \
		% (bdate, edate)
else:
	results = db.sql('select convert(varchar(25), dateadd(day, -7, "%s"))' % (currentDate), 'auto')
	bdate = results[0]['']

	results = db.sql('select convert(varchar(25), dateadd(day, 0, "%s"))' % (currentDate), 'auto')
	edate = results[0]['']

	datequery = 'dateadd(day,-7,"%s") and dateadd(day,1,"%s")' \
		% (currentDate, currentDate)

title = 'Updates to Mouse Nomenclature from %s to %s' % (bdate, edate)
fpHTML = reportlib.init(reportName, title, os.environ['FTPREPORTDIR'], isHTML = 1, printHeading = "MGI")
fpRpt = reportlib.init(reportName, outputdir = os.environ['FTPREPORTDIR'], printHeading = None)

fpHTML.write('J:23000 generally indicates gene family nomenclature revision event.\n\n')

fpHTML.write('%-2s %-30s %-35s %-10s %-30s %-25s %-75s %-15s %-25s\n' \
	% ('Ch', 'Symbol', 'Gene Name', 'J#', 'First Author', 'MGI ID', 'Sequence ID', 'Human Ortholog', 'Other MGI IDs'))

fpHTML.write('%-2s %-30s %-35s %-10s %-30s %-25s %-75s %-15s %-25s\n' \
	% ('--', '------', '---------', '--', '------------', '------', '-----------', '--------------', '-------------'))

db.sql('''
	select h._Marker_key, m.symbol, substring(m.name,1,35) as name, 
		substring(m.chromosome,1,2) as chromosome, 
		c.sequenceNum, 
		b.jnumID, 
		substring(b._primary,1,25) as author
	into #markers 
	from MRK_Marker m, MRK_History h, MRK_Chromosome c, BIB_All_View b 
	where m._Organism_key = 1 
	and m._Marker_key = h._History_key 
	and h.event_date between %s
	and m.chromosome = c.chromosome 
	and m._Organism_key = c._Organism_key 
	and h._Refs_key = b._Refs_key
	''' % (datequery), None)
db.sql('create index idx_key on #markers(_Marker_key)', None)

# get primary MGI ids (2)
results = db.sql('''
	select distinct m._Marker_key, a.accID 
	from #markers m, ACC_Accession a 
	where m._Marker_key = a._Object_key 
	and a._MGIType_key = 2 
	and a.prefixPart = 'MGI:' 
	and a._LogicalDB_key = 1 
	and a.preferred = 1
	''', 'auto')
primaryID = {}
for r in results:
	primaryID[r['_Marker_key']] = r['accID']

# other MGI ids (3)
results = db.sql('''
	select distinct m._Marker_key, a.accID 
	from #markers m, ACC_Accession a 
	where m._Marker_key = a._Object_key 
	and a._MGIType_key = 2 
	and a.prefixPart = 'MGI:' 
	and a._LogicalDB_key = 1 
	and a.preferred = 0
	''', 'auto')
otherIDs = {}
for r in results:
	if not otherIDs.has_key(r['_Marker_key']):
		otherIDs[r['_Marker_key']] = []
	otherIDs[r['_Marker_key']].append(r['accID'])

# get sequence ids (4)
results = db.sql('''
	select distinct m._Marker_key, a.accID 
	from #markers m, ACC_Accession a 
	where m._Marker_key = a._Object_key 
	and a._MGIType_key = 2 
	and a._LogicalDB_key = 9
	''', 'auto')
seqIDs = {}
for r in results:
	if not seqIDs.has_key(r['_Marker_key']):
		seqIDs[r['_Marker_key']] = []
	seqIDs[r['_Marker_key']].append(r['accID'])

# get human ortholog (5)
results = db.sql('''
	select distinct m._Marker_key, humanSymbol = h.symbol 
	from #markers m, MRK_Marker h, MRK_Homology_Cache hm1, MRK_Homology_Cache hm2 
	where m._Marker_key = hm1._Marker_key 
	and hm1._Homology_key = hm2._Homology_key 
	and hm2._Marker_key = h._Marker_key 
	and hm2._Organism_key = 1
	''', 'auto')
human = {}
for r in results:
	human[r['_Marker_key']] = r['humanSymbol']

# retrieve markers, sort (6)
results = db.sql('select distinct * from #markers order by sequenceNum, symbol', 'auto')
rows = 0
for r in results:

	key = r['_Marker_key']
	symbol = mgi_html.escape(r['symbol'])

	fpHTML.write('%-2s ' % (r['chromosome']))

	fpHTML.write('%s%-30s%s ' \
		% (reportlib.create_accession_anchor(primaryID[key], 'marker'), symbol, \
			reportlib.close_accession_anchor()))

	fpHTML.write('%-35s ' % (r['name']))

	fpHTML.write('%s%-10s%s ' \
		% (reportlib.create_accession_anchor(r['jnumID'], 'reference'), r['jnumID'], \
			reportlib.close_accession_anchor()))

	fpHTML.write('%-30s ' % (r['author']))

	fpHTML.write('%s%-25s%s ' \
		% (reportlib.create_accession_anchor(primaryID[key], 'marker'), primaryID[key], \
			reportlib.close_accession_anchor()))

	fpRpt.write('%s\t%s\t%s\t%s\t%s\t%s\t' \
		% (r['chromosome'], symbol, r['name'], r['jnumID'], r['author'], primaryID[key]))

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
db.useOneConnection(0)

# re-create a symbolic link between the new file and the current file
os.chdir(reportDir)
os.symlink(reportName + '.html', currentReport + '.html')
os.symlink(reportName + '.rpt', currentReport + '.rpt')

