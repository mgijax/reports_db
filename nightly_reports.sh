#!/bin/csh

#
# nightly_reports.sh
#
# Script to generate nightly reports.
#
# Usage: nightly_reports.sh
#

cd `dirname $0` && source Configuration

umask 002

foreach i (`ls *.sql`)
sql.sh $MGD $i
end

foreach i (`ls *.py`)
$i
end
