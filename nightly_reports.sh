#!/bin/csh

#
# nightly_reports.sh
#
# Script to generate nightly reports.
#
# Usage: nightly_reports.sh
#

cd `dirname $0` && source ./Configuration

umask 002

foreach i (daily/*.sql)
reportisql.csh $i $REPORTOUTPUTDIR/`basename $i`.rpt $DSQUERY $MGD
end

cd daily
foreach i (*.py)
$i
end
cd ..

cd anatdict
foreach i (*.py)
$i
end
cd ..

foreach i ($REPORTOUTPUTDIR/*)
rcp $i $FTPREPORTDIR
end

./mgimarkerfeed_reports.sh
