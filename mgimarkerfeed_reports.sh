#!/bin/csh

#
# mgimarkerfeed_reports.sh
#
# Script to generate mgi marker feed reports.
#
# Usage: mgimarkerfeed_reports.sh
#

setenv TARFILE mgimarkerfeed.tar

cd `dirname $0` && source ./Configuration

umask 002

cd daily
mgiMarkerFeed.py

cd $REPORTOUTPUTDIR/mgimarkerfeed
tar cvf $TARFILE *.bcp
compress -f $TARFILE
cp $TARFILE.Z $MGIFEEDFTPDIR
