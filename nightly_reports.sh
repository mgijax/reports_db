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
reportisql.csh $i $REPORTOUTPUTDIR/`basename $i`.rpt $DSQUERY $MGD
end

cd daily
foreach i (`ls *.py`)
$i
end
cd ..

mv $REPORTOUTPUTDIR/gene_association.rpt $REPORTOUTPUTDIR/gene_association.mgi
mv $REPORTOUTPUTDIR/go_terms.rpt $REPORTOUTPUTDIR/go_terms.mgi

foreach i (`ls $REPORTOUTPUTDIR`)
cp $REPORTOUTPUTDIR/$i $FTPREPORTDIR
end

./mgimarkerfeed_reports.sh
