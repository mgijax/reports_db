#!/bin/csh -f

#
# run_weekly_sybase.csh
#
# Script to generate weekly public sybase reports.
#
# Usage: run_weekly_sybase.csh
#

cd `dirname $0` && source ./Configuration

umask 002

#
# Initialize the log file.
#
setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start weekly public sybase reports | tee -a ${LOG}

#
# Generate weekly public sybase reports.
#
cd ${PUBWEEKLY_SYBASE}
foreach i (*.py)
    echo `date`: $i | tee -a ${LOG}
    $i >>& ${LOG}
end

#
# Generate inparanoid files.
#
echo `date`: inparanoid.csh | tee -a ${LOG}
${PUBRPTS}/inparanoid/inparanoid.csh >>& ${LOG}

echo `date`: End weekly public sybase reports | tee -a ${LOG}

exit 0
