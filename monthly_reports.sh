#!/bin/csh

#
# monthly_reports.sh
#
# Script to generate monthly qc reports.
#
# Usage: monthly_reports.sh
#

cd `dirname $0` && source ./Configuration

umask 002

setenv LOG	${REPORTLOGSDIR}/$0.log
rm -rf ${LOG}
touch ${LOG}

cd monthly
foreach i (`ls *.py`)
$i >>& ${LOG}
end

rcp ${REPORTOUTPUTDIR}/GO_nonmouse.rpt ${FTPCUSTOM}
rm -rf ${REPORTOUTPUTDIR}/GO_nonmouse.rpt
