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

foreach i (`ls daily/*.sql`)
sql.sh $MGD $i
end

cd daily
foreach i (`ls *.py`)
$i
end

foreach i (`ls $REPORTOUTPUTDIR`)
cp $REPORTOUTPUTDIR/$i $FTPREPORTDIR
end

./mgimarkerfeed.sh
