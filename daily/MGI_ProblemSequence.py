#!/usr/local/bin/python

'''
#
# MGI_ProblemSequence.py 09/18/2001
#
# Report:
#	TR 2943
#       Tab-delimited file
#
#	All Molecular Segments that contain the Problem Sequence Note. The file
#       should have the fields indicated below. This is the report for NCBI.
#
#       1. GenBank sequence identifier associated with the Molecular Segment 
#       2. MGI Accession ID(s) for the MGI Marker(s) to which the Molecular Segment 
#          has a Hybridizes relationship. If there is more than one Marker, delimit the MGI:IDs by commas. 
#
# Problem Sequence Note:
#
# "MGI curatorial staff have found evidence of artifact in the sequence of this molecular"
# "segment. The artifact may have been introduced during the cloning process, or may be"
# "due to unreliable sequence data."
#
# Usage:
#       MGI_ProblemSequence.py
#
# Used by:
#       NCBI
#
# Notes:
#
# History:
#
# lec	09/18/2001
#	- new
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

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

cmds = []

/* Select probes w/ problem note */
cmds.append('select _Probe_key ' + \
      'into #probes ' + \
      'from PRB_Notes ' + \
      'where note like "MGI curatorial staff have found evidence of artifact in the sequence of this molecular%"')

/* Select probes w/ Seq IDs */
cmds.append('select distinct p._Probe_key, a.accID ' + \
'into #probeseqs ' + \
'from #probes p, PRB_Acc_View a ' + \
'where p._Probe_key = a._Object_key  ' + \
'and a._LogicalDB_key = 9 ')

/* Select probes w/ only one Seq ID */
cmds.append('select * into #forncbi from #probeseqs group by accID having count(*) = 1')

/* Select probe's markers which hybridizie */
cmds.append('select n.*, markerID = ma.accID ' + \
'from #forncbi n, PRB_Marker_View m, MRK_Acc_View ma ' + \
'where n._Probe_key = m._Probe_key ' + \
'and m.relationship = "H" ' + \
'and m._Marker_key = ma._Object_key ' + \
'and ma.prefixPart = "MGI:" ' + \
'and ma.preferred = 1 ' + \
'order by n._Probe_key')

results = db.sql(cmds, 'auto')

prevProbe = ''
markers = []

for r in results[-1]:

	if prevProbe != r['_Probe_key']:
		if len(markers) > 0:
			fp.write(string.join(markers, ','))
		markers = ''

		if prevProbe != '':
			fp.write(reportlib.CRT)

		fp.write(mgi_utils.prvalue(r['accID']) + reportlib.TAB)

		prevProbe = r['_Probe_key']
		markers = []

        markers.append(r['markerID'])

fp.write(string.join(markers, ','))
reportlib.finish_nonps(fp)

