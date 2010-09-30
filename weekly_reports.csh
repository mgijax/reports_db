#!/bin/csh

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

echo `date`: Start nightly public reports | tee -a ${LOG}

cd ${PUBWEEKLY}

foreach i (`ls *.py`)
    echo `date`: $i | tee -a ${LOG}
    $i >>& ${LOG}
end

# Generate inparanoid files and copy them to the FTP site
echo `date`: inparanoid.csh | tee -a ${LOG}
cd ${PUBRPTS}
${PUBRPTS}/inparanoid/inparanoid.csh >>& ${LOG}

echo `date`: Copy reports | tee -a ${LOG}
rm -f ${FTPCUSTOM}/inparanoid/Mus-musculus*
cd ${INPARANOIDDIR}
cp Mus-musculus* aaseq* ${FTPCUSTOM}/inparanoid

echo `date`: End nightly public reports | tee -a ${LOG}
