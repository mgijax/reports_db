#!/bin/csh

#
# mgimarkerfeed_reports2.sh
#
# Script to generate mgi marker feed reports.
#
# Usage: mgimarkerfeed_reports2.sh
#

setenv TARFILE mgimarkerfeed2.tar

cd `dirname $0` && source ./Configuration

umask 002

cd daily
mgiMarkerFeed2.py

cd ${REPORTOUTPUTDIR}/mgimarkerfeed2
tar cvf ${TARFILE} *.bcp
compress -f ${TARFILE}
rcp ${TARFILE}.Z ${MGIFEEDFTPDIR}
