#!/bin/csh -f

#
# ondemand_reports.csh
#
# Script to generate on-demand public reports.
#
# Usage: ondemand_reports.csh
#

cd `dirname $0` && source ./Configuration

umask 002

setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start on-demand public reports | tee -a ${LOG}

cd ${PUBWEEKLY}

foreach i (*.py)
    echo `date`: $i | tee -a ${LOG}
    $i >>& ${LOG}
end

echo `date`: End weekly public reports | tee -a ${LOG}
