#!/usr/local/bin/python

'''
#
# TR 1604 - TEXT/CROSS experiments for JAX BSS/BSB crosses
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
# lec	06/19/2000
#	- added date contraints...this report is expected to be run on the
#	  first of every month
#
# lec	06/09/2000
#	- created
#
'''
 
import sys 
import string
import os
import db
import reportlib

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Main
#

fp = reportlib.init(sys.argv[0], printHeading = 0, outputdir = os.environ['REPORTOUTPUTDIR'])

cmds = []

cmds.append('select _Cross_key ' + \
'into #cross ' + \
'from CRS_Cross ' + \
'where whoseCross like "%jax%bss%" or whoseCross like "%jax%bsb%"')

cmds.append('select distinct e._Refs_key, e._Expt_key ' + \
'into #expts ' + \
'from #cross c, MLD_Matrix m, MLD_Expts e ' + \
'where c._Cross_key = m._Cross_key ' + \
'and m._Expt_key = e._Expt_key ' + \
'and datepart(year, e.creation_date) = datepart(year, getdate()) ' + \
'and datepart(month, e.creation_date) = datepart(month, getdate()) - 1 ' + \
'union ' + \
'select distinct e._Refs_key, e._Expt_key ' + \
'from MLD_Expts e, MLD_Expt_Notes n ' + \
'where e.exptType like "TEXT%" ' + \
'and datepart(year, e.creation_date) = datepart(year, getdate()) ' + \
'and datepart(month, e.creation_date) = datepart(month, getdate()) - 1 ' + \
'and e._Expt_key = n._Expt_key ' + \
'and (n.note like "%jax%bss%" or n.note like "%jax%bsb%" ' + \
'or n.note like "%jackson laboratory%bss%" or n.note like "%jackson laboratory%bsb%")')

cmds.append('select e.*, n.note, n.sequenceNum ' + \
'into #expts2 ' + \
'from #expts e, MLD_Expt_Notes n ' + \
'where e._Expt_key *= n._Expt_key ' + \
'order by n._Expt_key, n.sequenceNum ')

cmds.append('select e.*, m.symbol ' + \
'into #expts3 ' + \
'from #expts2 e, MLD_Expt_Marker em, MRK_Marker m ' + \
'where e._Expt_key = em._Expt_key ' + \
'and em._Marker_key = m._Marker_key')

cmds.append('select e.*, b.jnumID, b.authors, b.title, b.citation ' + \
'from #expts3 e, BIB_All_View b ' + \
'where e._Refs_key = b._Refs_key ' + \
'order by b.jnumID, e._Expt_key, e.sequenceNum')

results = db.sql(cmds, 'auto')

prevExpt = ''
markers = []
notes = []

for r in results[4]:

	if prevExpt != r['_Expt_key']:

		if prevExpt != '':
			fp.write(string.join(markers, ',') + TAB)

			if len(notes) > 0:
				fp.write(string.join(notes, ''))

			fp.write(CRT)
			markers = []
			notes = []

		prevExpt = r['_Expt_key']

		fp.write(r['jnumID'] + TAB)
		fp.write(r['authors'] + SPACE +
			 r['title'] + SPACE +
			 r['citation'] + TAB)

	if not r['symbol'] in markers:
		markers.append(r['symbol'])

	if r['note'] != None and not r['note'] in notes:
		notes.append(r['note'])

if len(markers) > 0:
	fp.write(string.join(markers, ',') + TAB)

if len(notes) > 0:
	fp.write(string.join(notes, ' '))

fp.write(CRT)

reportlib.finish_nonps(fp)	# non-postscript file

