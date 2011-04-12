#!/bin/csh

#
# nightly_reports.csh
#
# Script to generate nightly public reports.
#
# Usage: nightly_reports.csh
#

cd `dirname $0` && source ./Configuration

umask 002

setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start nightly public reports | tee -a ${LOG}

cd ${PUBDAILY}

echo `date`: GO_gene_association.py | tee -a ${LOG}
GO_gene_association.py >>& ${LOG}

echo `date`: Copy reports | tee -a ${LOG}
cp ${REPORTOUTPUTDIR}/gene_association.mgi ${FTPREPORTDIR}

echo `date`: End nightly public reports | tee -a ${LOG}
