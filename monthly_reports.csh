#!/bin/csh

#
# monthly_reports.csh
#
# Script to generate monthly public reports.
#
# Usage: monthly_reports.csh
#

cd `dirname $0` && source ./Configuration

umask 002

setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start monthly public reports | tee -a ${LOG}

cd ${PUBMONTHLY}

foreach i (`ls *.py`)
    echo `date`: $i | tee -a ${LOG}
    $i >>& ${LOG}
end

# Generate NCBI LinkOut files
echo `date`: ncbilinkout.csh | tee -a ${LOG}
cd ${PUBRPTS}
${PUBRPTS}/ncbilinkout/ncbilinkout.csh >>& ${LOG}

echo `date`: Copy reports | tee -a ${LOG}
cp ${REPORTOUTPUTDIR}/gene_association.mgi_nonmouse ${FTPCUSTOM}
rm -rf ${REPORTOUTPUTDIR}/gene_association.mgi_nonmouse

cp ${REPORTOUTPUTDIR}/ES_CellLine.rpt ${FTPREPORTDIR}
cp ${REPORTOUTPUTDIR}/MRK_GeneTrap.rpt ${FTPREPORTDIR}

echo `date`: End monthly public reports | tee -a ${LOG}
