#!/bin/csh -f

#
# run_daily_sybase.csh
#
# Script to generate daily public sybase reports.
#
# Usage: run_daily_sybase.csh
#

cd `dirname $0` && source ./Configuration

umask 002

#
# Set this variable to either 'sybase' or 'postgres'.
#
setenv DB_TYPE                  sybase

#
# Initialize the log file.
#
setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start daily public sybase reports | tee -a ${LOG}

#
# Generate daily public sybase reports.
#
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

echo `date`: End nightly public sybase reports | tee -a ${LOG}

exit 0
