#############################################################################
#
# TAL_MRK_Nomen_Change.py
#
# Usage:
#       ${PYTHON} TAL_MRK_Nomen_Change.py
#
# Report:
#       Weekly report from the TAL input file of Marker Nomen changes
#
# History:
#
# sc	10/17/2023
#	- created FL2-564/WTS2-558
#
# Notes:
#
# For requirements see: 
# https://mgi-jira.atlassian.net/browse/FL2-564
#
#############################################################################
 
import sys 
import os
import string
import reportlib
import db

db.setTrace()
db.useOneConnection(1)

CRT = reportlib.CRT
TAB = reportlib.TAB

fpOut = reportlib.init(sys.argv[0], fileExt = '.last30days.rpt', outputdir = os.environ['REPORTOUTPUTDIR'])
fpOut.write('Marker MGI ID%sCurrent Marker Symbol%sPrevious Marker Symbol%sEvent%sEvent Date%s' % (TAB, TAB, TAB, TAB, CRT))

fpIn = open(os.getenv('TAL_FILE'), 'r')

# {mgiID: results, ...}
dbEvents = {}

results = db.sql('''select distinct aa._marker_key, a.accid as markerID, 
        a.numericPart, hv.symbol as currentSymbol, hv.history as oldSymbol, 
        hv.event, to_char(hv.event_date, 'MM/dd/yyyy') as eventDate
    from mrk_history_view hv, all_allele aa, acc_accession a
    where aa._marker_key = hv._marker_key
    and hv._marker_event_key in (
        106563605, 106563606, 106563608, 106563609)
    and hv.event_date >= (now() + interval '-30 day') and hv.event_date <= now()
    and hv._marker_key = a._object_key
    and a._mgitype_key = 2
    and a._logicaldb_key = 1
    and a.preferred = 1
    order by a.numericPart, eventDate''', 'auto')

for r in results:
    markerID = r['numericPart']

    if markerID not in dbEvents:
        dbEvents[markerID] = []
    dbEvents[markerID].append(r)

nomenCt = 0
talMarkersList = []
header = fpIn.readline()
for line in fpIn.readlines():
    tokens = str.split(line, TAB)
    if len(tokens) < 8: # for blank lines
        continue
    
    pid = tokens[4]
    if pid == '':
        continue

    mgiID = tokens[0][4:]
    mgiID = int(mgiID)
    if mgiID not in talMarkersList:
        talMarkersList.append(mgiID)
talMarkersList.sort()

for mgiID in talMarkersList:
    if  mgiID in dbEvents:
        eventList = dbEvents[mgiID]
        for r in eventList:
            nomenCt += 1
            fpOut.write('MGI:%s%s%s%s%s%s%s%s%s%s' % (mgiID, TAB, r['currentSymbol'], TAB, r['oldSymbol'], TAB, r['event'], TAB, r['eventDate'], CRT))

fpOut.write('%sTotal: %s' % (CRT, nomenCt))
fpIn.close()
reportlib.finish_nonps(fpOut)	# non-postscript file

db.useOneConnection(0)

