#!/bin/csh -f

#
# weekly_reports.csh
#
# Script to generate weekly public reports.
#
# Usage: weekly_reports.csh
#

cd `dirname $0` && source ./Configuration

umask 002

setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start weekly public reports | tee -a ${LOG}

cd ${PUBWEEKLY}

foreach i (*.py)
    echo `date`: $i | tee -a ${LOG}
    $i >>& ${LOG}
end
exit 0

#
# Generate inparanoid files.
#
echo `date`: inparanoid.csh | tee -a ${LOG}
${PUBRPTS}/inparanoid/inparanoid.csh >>& ${LOG}

#
# Generate NCBI LinkOut files.
#
echo `date`: ncbilinkout.csh | tee -a ${LOG}
${PUBRPTS}/ncbilinkout/ncbilinkout.csh >>& ${LOG}

echo `date`: End weekly public reports | tee -a ${LOG}
