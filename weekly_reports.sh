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

# remove old inParanoid reports
rm -rf ${REPORTOUTPUTDIR}/Mus-musculus_MGI*

cd weekly
foreach i (`ls *.py`)
$i >>& ${LOG}
end

