#!/bin/csh -f

#
# run_weekly.csh
#
# Script to generate weekly public reports.
#
# Usage: run_weekly.csh
#

cd `dirname $0` && source ./Configuration

umask 002

#
# Set postgres configuration to use the public instance.
#
setenv PG_DBSERVER	${PG_PUB_DBSERVER}
setenv PG_DBNAME	${PG_PUB_DBNAME}
setenv PG_DBUSER	${PG_PUB_DBUSER}

#
# Initialize the log file.
#
setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start weekly public reports | tee -a ${LOG}

#
# Generate weekly public reports.
#
cd ${PUBWEEKLY}
foreach i (*.py)
    echo `date`: $i | tee -a ${LOG}
    $i >>& ${LOG}
end

#
# Generate NCBI LinkOut files.
#
echo `date`: ncbilinkout.csh | tee -a ${LOG}
${PUBRPTS}/ncbilinkout/ncbilinkout.csh >>& ${LOG}

#
# Generate inparanoid files.
#
echo `date`: inparanoid.csh | tee -a ${LOG}
${PUBRPTS}/inparanoid/inparanoid.csh >>& ${LOG}

echo `date`: End weekly public reports | tee -a ${LOG}

exit 0
