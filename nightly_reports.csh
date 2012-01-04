#!/bin/csh -f

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

foreach i (*.py)
    echo `date`: $i | tee -a ${LOG}
    $i >>& ${LOG}
end

cd ${REPORTOUTPUTDIR}

echo `date`: Copy reports | tee -a ${LOG}
foreach i (gene_association.mgi)
    echo `date`: $i | tee -a ${LOG}
    cp $i ${FTPREPORTDIR}
end

echo `date`: End nightly public reports | tee -a ${LOG}
