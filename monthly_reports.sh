#!/bin/csh

#
# monthly_reports.sh
#
# Script to generate monthly qc reports.
#
# Usage: monthly_reports.sh
#

cd `dirname $0` && source ./Configuration

umask 002

setenv LOG	${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

cd monthly
foreach i (`ls *.py`)
$i >>& ${LOG}
end

# Generate NCBI LinkOut files
cd ${PUBRPTS}
${PUBRPTS}/ncbilinkout/ncbilinkout.sh >>& ${LOG}

cp ${REPORTOUTPUTDIR}/gene_association.mgi_nonmouse ${FTPCUSTOM}
rm -rf ${REPORTOUTPUTDIR}/gene_association.mgi_nonmouse

cp ${REPORTOUTPUTDIR}/ES_CellLine.rpt ${FTPREPORTDIR}
cp ${REPORTOUTPUTDIR}/MRK_GeneTrap.rpt ${FTPREPORTDIR}
