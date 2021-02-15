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
# Initialize the log file.
#
setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start weekly public reports | tee -a ${LOG}

#
# unzip files
#
echo "unzipping input file ${ALLIANCE_HUMAN_FILE_GZ}" >> ${LOG}
gunzip -cf ${ALLIANCE_HUMAN_FILE_GZ} > ${ALLIANCE_HUMAN_FILE}

#
# Generate weekly public reports.
#
cd ${PUBWEEKLY}
foreach i (*.py)
    echo `date`: $i | tee -a ${LOG}
    ${PYTHON} $i >>& ${LOG}
end

#
# Generate NCBI LinkOut files.
#
echo `date`: ncbilinkout.csh | tee -a ${LOG}
${PUBRPTS}/ncbilinkout/ncbilinkout.csh >>& ${LOG}

echo `date`: End weekly public reports | tee -a ${LOG}

exit 0
