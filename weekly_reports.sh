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

cd weekly
foreach i (`ls *.py`)
$i
end

