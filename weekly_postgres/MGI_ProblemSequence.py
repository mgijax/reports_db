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
#	2. Molecular Segment Name
#       3. MGI Accession ID(s) for the MGI Marker(s) to which the Molecular Segment 
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
import mgi_utils
import reportlib
import pg_db
db = pg_db
db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE()

#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

# Select probes w/ problem note
db.sql('''
	select p._Probe_key, p.name 
        into #probes 
        from PRB_Notes n, PRB_Probe p 
        where n.note like '%staff have found evidence of artifact in the sequence of this molecular%'
        and n._Probe_key = p._Probe_key
	''', None)
db.sql('create index probes_idx1 on #probes(_Probe_key)', None)

# Select probes w/ Seq IDs and without Seq IDs
db.sql('''
	(
	select distinct p._Probe_key, p.name, a.accID 
	into #probeseqs 
	from #probes p, ACC_Accession a 
	where p._Probe_key = a._Object_key  
	and a._MGIType_key = 3 
	and a._LogicalDB_key = 9 
	union 
	select distinct p._Probe_key, p.name, null as accID
	from #probes p, ACC_Accession pa 
	where p._Probe_key = pa._Object_key 
	and pa._MGIType_key = 3 
	and pa.prefixPart = 'MGI:' 
	and pa._LogicalDB_key = 1 
	and not exists (select 1 from ACC_Accession a 
	where p._Probe_key = a._Object_key  
	and a._MGIType_key = 3 
	and a._LogicalDB_key = 9) 
	)
	order by _Probe_key
	''', None)
db.sql('create index probeseqs_idx1 on #probeseqs(_Probe_key)', None)

# Select probes w/ only one Seq ID
db.sql('select _Probe_key into #forncbi from #probeseqs group by _Probe_key having count(*) = 1', None)
db.sql('create index forncbi_idx1 on #forncbi(_Probe_key)', None)

# Select probe's markers which hybridizie
results = db.sql('''
	select p._Probe_key, p.name, p.accID, ma.accID as markerID
	from #forncbi n, #probeseqs p, PRB_Marker m, ACC_Accession ma 
	where n._Probe_key = p._Probe_key
	and p._Probe_key = m._Probe_key 
	and m.relationship = 'H' 
	and m._Marker_key = ma._Object_key 
	and ma._MGIType_key = 2 
	and ma.prefixPart = 'MGI:' 
	and ma._LogicalDB_key = 1 
	and ma.preferred = 1 
	order by p.accID
	''', 'auto')

prevProbe = 0
markers = []

for r in results:

	if prevProbe != r['_Probe_key']:
		if len(markers) > 0:
			fp.write(string.join(markers, ','))
		markers = ''

		if prevProbe > 0:
			fp.write(reportlib.CRT)

		fp.write(mgi_utils.prvalue(r['accID']) + reportlib.TAB)
		fp.write(mgi_utils.prvalue(r['name']) + reportlib.TAB)

		prevProbe = r['_Probe_key']
		markers = []

        markers.append(r['markerID'])

fp.write(string.join(markers, ','))
fp.write(reportlib.CRT)
reportlib.finish_nonps(fp)

