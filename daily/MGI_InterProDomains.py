#!/usr/local/bin/python

'''
#
# Report:
#       Tab-delimited file
#       	Interpro Acc Id
#		Domain Name
#
# Usage:
#       MGI_InterProDomains.py
#
#
# Notes:
#
# History:
#
# lec	03/26/2003
#	- TR 3702
#
'''
 
import sys
import os
import db
import mgi_utils
import reportlib

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmd = 'select a.accID, t.term ' + \
      'from VOC_Vocab v, VOC_Term t, VOC_Term_Acc_View a ' + \
      'where v.name = "InterPro Domains" ' + \
      'and v._Vocab_key = t._Vocab_key ' + \
      'and t._Term_key = a._Object_key ' + \
      'order by t.sequenceNum'

results = db.sql(cmd, 'auto')

for r in results:

	fp.write(r['accID'] + reportlib.TAB + \
	         r['term'] + reportlib.CRT)

reportlib.finish_nonps(fp)

