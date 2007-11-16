#!/bin/csh

#
# weekly_reports.sh
#
# Script to generate weekly qc reports.
#
# Usage: weekly_reports.sh
#

cd `dirname $0` && source ./Configuration

umask 002

setenv LOG	${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

cd weekly
foreach i (`ls *.py`)
$i >>& ${LOG}
end

# Generate inparanoid files and copy them to the FTP site
cd ${PUBRPTS}
${PUBRPTS}/inparanoid/inparanoid.sh >>& ${LOG}

cd ${INPARANOIDDIR}
cp Mus-musculus* aaseq* ${FTPCUSTOM}/inparanoid
