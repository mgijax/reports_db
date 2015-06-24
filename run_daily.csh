#!/bin/csh -f

#
# run_daily.csh
#
# Script to generate daily public reports.
#
# Usage: run_daily.csh
#

cd `dirname $0` && source ./Configuration

umask 002

#
# Initialize the log file.
#
setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start daily public reports | tee -a ${LOG}

#
# Generate daily public reports.
#
cd ${PUBDAILY}
foreach i (*.py)
    echo `date`: $i | tee -a ${LOG}
    $i >>& ${LOG}
end

#
# Copy reports to ftp site
#
cd ${REPORTOUTPUTDIR}

echo `date`: Copy reports | tee -a ${LOG}
foreach i (gene_association.mgi)
    echo `date`: $i | tee -a ${LOG}
    cp $i ${FTPREPORTDIR}
end
foreach i (GO_eco_association.rpt)
    echo `date`: $i | tee -a ${LOG}
    cp $i ${FTPREPORTDIR}
end

echo `date`: End daily public reports | tee -a ${LOG}

exit 0
