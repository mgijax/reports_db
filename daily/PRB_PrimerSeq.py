#!/usr/local/bin/python

'''
#
# PRB_PrimerSeq.py 11/05/99
#
# Report:
#       TR 1047 Primers and Sequences
#
# Usage:
#       PRB_PrimerSeq.py
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
# lec	11/05/1999
#	- created
#
# dbm	02/11/2003
#	- added marker AccID, chromosome and offset for TR 4506
#
# dbm	02/12/2003
#	- added marker name
#
'''
 
import sys 
import os
import db
import reportlib
import string
import regsub

CRT = reportlib.CRT
TAB = reportlib.TAB

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmd = 'select m.symbol, m.name "mname", p.name "pname", ' + \
             'a.accID, p.mgiID, ' + \
             'p.primer1sequence, p.primer2sequence, ' + \
             'p.productSize, m.chromosome, o.offset ' + \
      'from PRB_Primer_View p, PRB_Marker pm, MRK_Marker m, ' + \
           'ACC_Accession a, MRK_Offset o ' + \
      'where p._Probe_key = pm._Probe_key and ' + \
            'pm._Marker_key = m._Marker_key and ' + \
            'm._Marker_key = a._Object_key and ' + \
	    'and a._MGIType_key = 2 ' + \
            'a.prefixPart = "MGI:" and ' + \
            'a.preferred = 1 and ' + \
            'a._LogicalDB_key = 1 and ' + \
            'm._Marker_key = o._Marker_key and ' + \
            'o.source = 0 ' + \
      'order by m.symbol'

results = db.sql(cmd, 'auto')

for r in results:
    mname = r['mname']
    pname = r['pname']
    p1seq = r['primer1sequence']
    p2seq = r['primer2sequence']
    prodSize = r['productSize']

    if (mname == None):
        mname = ""
    else:
        mname = string.strip(regsub.gsub('\n', '', mname))
    if (pname == None):
        pname = ""
    else:
        pname = string.strip(regsub.gsub('\n', '', pname))
    if (p1seq == None):
        p1seq = ""
    else:
        p1seq = string.strip(regsub.gsub('\n', '', p1seq))
    if (p2seq == None):
        p2seq = ""
    else:
        p2seq = string.strip(regsub.gsub('\n', '', p2seq))
    if (prodSize == None):
        prodSize = ""

    fp.write(r['symbol'] + TAB + mname + TAB + pname + TAB +
             r['accID'] + TAB + r['mgiID'] + TAB +
             p1seq + TAB + p2seq + TAB + prodSize + TAB +
             r['chromosome'] + TAB + str(r['offset']) + CRT)

reportlib.finish_nonps(fp)
