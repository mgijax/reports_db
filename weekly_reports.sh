#!/bin/csh

#
# weekly_reports.sh
#
# Script to generate weekly qc reports.
#
# Usage: weekly_reports.sh
#

cd `dirname $0` && source Configuration

umask 002

set RPTS="MRK_NomenUpdates.py"

foreach i ("$RPTS")
$i
end

