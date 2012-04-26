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

#
# Generate weekly reports
#
cd ${PUBWEEKLY}
foreach i (*.py)
    echo `date`: $i | tee -a ${LOG}
    $i >>& ${LOG}
end

#
# Generate on-demain reports, if set
#
if (${RUN_ONDEMAND} == 1)
	${PUBRPTS}/ondemand.csh >>& ${LOG}
endif

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
