#!/bin/csh

#
# nightly_reports.sh
#
# Script to generate nightly reports.
#
# Usage: nightly_reports.sh
#

cd `dirname $0` && source ./Configuration

umask 002

setenv LOG	${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

# GOA processing
./goa.sh

foreach i (daily/*.sql)
echo $i, `date`
reportisql.csh $i ${REPORTOUTPUTDIR}/`basename $i`.rpt ${MGD_DBSERVER} ${MGD_DBNAME} >> ${LOG}
echo $i, `date`
end

cd daily

# exclude mgiMarkerFeed and clone reports 

foreach i (*.py)
if ( $i != "mgiMarkerFeed.py" && $i != "mgiMarkerFeed2.py" && $i != "MGI_CloneSet.py" ) then
echo $i, `date`
$i >>& ${LOG}
echo $i, `date`
endif
end

# clone reports

foreach i ("Image", "NIA 15K,NIA 7.4K,NIA", "RIKEN (FANTOM),RIKEN", "RPCI-23", "RPCI-24")
echo $i, `date`
./MGI_CloneSet.py "$i" >> ${LOG}
echo $i, `date`
end

cd ..

cd anatdict
foreach i (*.py)
$i >> ${LOG}
end
cd ..

foreach i (${REPORTOUTPUTDIR}/*)
rcp $i ${FTPREPORTDIR}
end

# this is now being run as a separate task from mgidbutilities/bin/prod/dailytasks.csh
# so that it runs earlier in the evening (per Carolyn Blake)
#
#./mgimarkerfeed_reports.sh
#

