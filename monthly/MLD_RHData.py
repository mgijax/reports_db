#!/usr/local/bin/python

'''
#
# TR 1683 - RH Data for Mary Barter (meb@jax.org)
#	  this report is expected to run the first of every month
#
# Notes:
#	- all reports use mgireport directory for output file
#	- all reports use db default of public login
#	- all reports use server/database default of environment
#	- use lowercase for all SQL commands (i.e. select not SELECT)
#	- all public SQL reports require the header and footer
#	- all private SQL reports require the header
#
# History:
#
# lec	06/27/2000
#	- created
#
'''
 
import sys 
import string
import os
import db
import reportlib
import mgi_utils

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

fp = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])

cmds = []

cmds.append('select distinct e._Refs_key, e._Expt_key ' + \
'into #expts ' + \
'from MLD_Expts e ' + \
'where e.exptType = "TEXT-Radiation Hybrid" ' + \
'and datepart(year, e.creation_date) = datepart(year, getdate()) ' + \
'and datepart(month, e.creation_date) = datepart(month, getdate()) - 1 ' + \
'union ' + \
'select distinct e._Refs_key, e._Expt_key ' + \
'from MLD_Expts e, MLD_Expt_Notes n ' + \
'where e.exptType in ("TEXT", "TEXT-Cytogenetic Localization", "TEXT-Physical Mapping", "TEXT-Genetic Cross", "CROSS") ' + \
'and datepart(year, e.creation_date) = datepart(year, getdate()) ' + \
'and datepart(month, e.creation_date) = datepart(month, getdate()) - 1 ' + \
'and e._Expt_key = n._Expt_key ' + \
'and (n.note like "%radiation hybrid%" or n.note like "%T31%" or n.note like "%RH%")')

cmds.append('select e.*, n.note, n.sequenceNum ' + \
'into #expts2 ' + \
'from #expts e, MLD_Expt_Notes n ' + \
'where e._Expt_key *= n._Expt_key ' + \
'order by n._Expt_key, n.sequenceNum ')

cmds.append('select e.*, b.jnumID, b.authors, b.title, b.citation ' + \
'from #expts2 e, BIB_All_View b ' + \
'where e._Refs_key = b._Refs_key ' + \
'order by b.jnumID, e._Expt_key, e.sequenceNum')

results = db.sql(cmds, 'auto')

prevExpt = ''
notes = []

for r in results[2]:

	if prevExpt != r['_Expt_key']:

		if prevExpt != '':
			if len(notes) > 0:
				fp.write(string.join(notes, ''))

			fp.write(CRT)
			notes = []

		prevExpt = r['_Expt_key']

		fp.write(r['jnumID'] + TAB)
		fp.write(mgi_utils.prvalue(r['authors']))
		fp.write(' ' + mgi_utils.prvalue(r['title']))
		fp.write(mgi_utils.prvalue(r['citation'] + TAB))

	if r['note'] != None and not r['note'] in notes:
		notes.append(r['note'])

if len(notes) > 0:
	fp.write(string.join(notes, ' '))

fp.write(CRT)

reportlib.finish_nonps(fp)	# non-postscript file

