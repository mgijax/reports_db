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
echo $i, `date`
reportisql.csh $i $REPORTOUTPUTDIR/`basename $i`.rpt $DSQUERY $MGD
echo $i, `date`
end

cd daily

# exclude mgiMarkerFeed and clone reports 

foreach i (*.py)
if ( $i != "mgiMarkerFeed.py" && $i != "PRB_CloneSet.py" ) then
echo $i, `date`
$i
echo $i, `date`
endif
end

# clone reports

foreach i ("Image", "NIA 15K,NIA 7.4K,NIA", "RIKEN (FANTOM),RIKEN", "RPCI-23", "RPCI-24")
echo $i, `date`
./PRB_CloneSet.py "$i"
echo $i, `date`
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
