#!/bin/csh -f

#
# ondemand_reports.csh
#
# Script to generate on-demand public reports.
#
# Usage: ondemand_reports.csh
#

cd `dirname $0` && source ./Configuration

#
# set postgres configuration to use the public instance
#
setenv PG_DBSERVER      ${PG_PUB_DBSERVER}
setenv PG_DBNAME        ${PG_PUB_DBNAME}
setenv PG_DBUSER        ${PG_PUB_DBUSER}

echo "PG_DBSERVER: ${PG_DBSERVER}"
echo "PG_DBNAME: ${PG_DBNAME}"
echo "PG_DBUSER: ${PG_DBUSER}"

umask 002

setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start on-demand public reports | tee -a ${LOG}

cd ${ONDEMAND}

foreach i (*.py)
    echo `date`: $i | tee -a ${LOG}
    $i >>& ${LOG}
end

echo `date`: End on-demand public reports | tee -a ${LOG}
