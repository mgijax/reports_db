#!/usr/local/bin/python

'''
#
# PRB_AltTranscripts.py 04/13/99
#
# Report:
#	TR 108
# 	This program generates the Alternative Transcripts public report.
# 	The report is a list of genes that are known to produce alternate
# 	transcripts.  
#
# Usage:
#       PRB_AltTranscripts.py
#
# Notes:
#       - all reports use db default of public login
#       - all reports use server/database default of environment
#       - use lowercase for all SQL commands (i.e. select not SELECT)
#       - all public SQL reports require the header and footer
#       - all private SQL reports require the header
#
# History:
#
# lec   04/13/1999
#       - migrated to new template format
#
# cjd	12/1998
#	- created
#
'''

import sys
import os
import string
import db
import reportlib

CRT = reportlib.CRT
SPACE = reportlib.SPACE

fp = reportlib.init(sys.argv[0], "Alternate Transcripts Report", os.environ['REPORTOUTPUTDIR'])
fp.write(
'''
For the past 10 years, MGI editorial staff have attached a standard
notation to molecular probe records if a probe has been characterized by
researchers as "derived from a gene that produces alternate transcripts".
The Alternate Trancripts Report is based on this information. For further
details, see: Mangan, M.E. and Frazer, K.S.  1999.  An extensive list of genes that
produce alternative transcripts in the mouse.  Bioinformatics 15(2):170-171.
 
Please note: Molecular probes / segments described in one paper to
illustrate alternative transcription may also be used in other papers which
do not mention these characteristics.  Therefore, in cases where multiple
references (jnumbIDs) are associated with one mouse gene symbol, at least
one but not necessarily all references provide evidence for alternative
transcription for that gene.

We provide gene symbols for human orthologs but we do not store or provide
information about alternative transcript data for human genes.  When a
human symbol is unavailable, this is designated as NULL.  Pending symbols
are designated with the suffix "-pending" (mouse) or "-PEN" (human).
Unapproved human symbols are enclosed by asterisks.
''' + CRT
)

fp.write(string.ljust('mouse_symbol', 25))
fp.write(SPACE)
fp.write(string.ljust('chromosome', 10))
fp.write(SPACE)
fp.write(string.ljust('jnumID', 10))
fp.write(SPACE)
fp.write(string.ljust('human_symbol', 25))
fp.write(CRT)

fp.write(string.ljust('-------------------------', 25))
fp.write(SPACE)
fp.write(string.ljust('----------', 10))
fp.write(SPACE)
fp.write(string.ljust('----------', 10))
fp.write(SPACE)
fp.write(string.ljust('-------------------------', 25))
fp.write(CRT)

cmds = []
cmds.append('select distinct  mouse_symbol = mm.symbol, ' + \
            'mm.chromosome, mm._Marker_key, pb._Probe_key, pr._Refs_key ' + \
            'into #csnr ' + \
            'from PRB_Reference pr, PRB_Source ps, PRB_Probe pb, PRB_Notes pn, ' + \
	    'MRK_Marker mm, PRB_Marker pm ' + \
            'where pn.note like "%alternat%" ' + \
            'and ps._ProbeSpecies_key = 1 ' + \
            'and pm.relationship = "E" ' + \
            'and pn._Probe_key = pb._Probe_key ' + \
            'and pb._Probe_key = pr._Probe_key ' + \
            'and pb._Source_key = ps._Source_key ' + \
            'and mm._Marker_key = pm._Marker_key ' + \
            'and pm._Probe_key = pb._Probe_key')

cmds.append('select distinct  c.mouse_symbol, c._Marker_key, human_symbol= m1.symbol ' + \
            'into #homol ' + \
            'from #csnr c, MRK_Marker m1, HMD_Homology r1, ' + \
            'HMD_Homology_Marker h1, HMD_Homology r2, ' + \
            'HMD_Homology_Marker h2 ' + \
            'where m1._Species_key = 2 ' + \
            'and m1._Marker_key = h1._Marker_key ' + \
            'and h1._Homology_key = r1._Homology_key ' + \
            'and r1._Class_key = r2._Class_key ' + \
            'and r2._Homology_key = h2._Homology_key ' + \
            'and h2._Marker_key = c._Marker_key')

cmds.append('select distinct c.mouse_symbol, c.chromosome, a.jnumID, h.human_symbol ' + \
            'from #csnr c, #homol h, BIB_All_View a ' + \
            'where c._Refs_key = a._Refs_key ' + \
            'and c._Marker_key *= h._Marker_key ' + \
            'order by c.mouse_symbol')

results = db.sql(cmds, 'auto')

prevSymbol = ''

for r in results[2]:

	if r['human_symbol'] is None:
		r['human_symbol'] = "NULL"

	if prevSymbol != r['mouse_symbol']:
		fp.write(string.ljust(r['mouse_symbol'], 25))
		fp.write(SPACE)
		fp.write(string.ljust(r['chromosome'], 10))
		fp.write(SPACE)
		fp.write(string.ljust(r['jnumID'], 10))
		fp.write(SPACE)
		fp.write(string.ljust(r['human_symbol'], 25))
		fp.write(CRT)
	else:
		fp.write(string.ljust(SPACE, 25))
		fp.write(SPACE)
		fp.write(string.ljust(SPACE, 10))
		fp.write(SPACE)
		fp.write(string.ljust(r['jnumID'], 10))
		fp.write(CRT)

	prevSymbol = r['mouse_symbol']

reportlib.trailer(fp)
reportlib.finish_nonps(fp)

