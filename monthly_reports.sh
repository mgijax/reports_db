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

cd monthly
foreach i (`ls *.py`)
$i
end

rcp ${REPORTOUTPUTDIR}/GO_nonmouse.rpt ${FTPCUSTOM}
rm -rf ${REPORTOUTPUTDIR}/GO_nonmouse.rpt
