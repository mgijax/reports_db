#!/bin/csh

#
# nightly_reports.csh
#
# Script to generate nightly public reports.
#
# Usage: nightly_reports.csh
#

cd `dirname $0` && source ./Configuration

umask 002

setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start nightly public reports | tee -a ${LOG}

cd ${PUBDAILY}

foreach i (*.sql)
    echo `date`: $i | tee -a ${LOG}
    reportisql.csh $i ${REPORTOUTPUTDIR}/$i.rpt ${MGD_DBSERVER} ${MGD_DBNAME} MGI >> ${LOG}
end

# exclude mgiMarkerFeed and clone reports 

foreach i (*.py)
    if ( $i != "mgiMarkerFeed.py" && $i != "MGI_CloneSet.py" ) then
        echo `date`: $i | tee -a ${LOG}
        $i >>& ${LOG}
    endif
end

# clone reports

foreach i ("Image" "NIA 15K,NIA 7.4K,NIA" "RIKEN (FANTOM),RIKEN" "RPCI-23" "RPCI-24")
    echo `date`: MGI_CloneSet.py $i | tee -a ${LOG}
    ./MGI_CloneSet.py "$i" >>& ${LOG}
end

echo `date`: Copy reports | tee -a ${LOG}
foreach i (${REPORTOUTPUTDIR}/*)
    if ( ! -d $i ) then
        cp $i ${FTPREPORTDIR}
    endif
end

echo `date`: End nightly public reports | tee -a ${LOG}
