#!/bin/csh -f

#
# run_daily_sybase.csh
#
# Script to generate daily public sybase reports.
#
# Usage: run_daily_sybase.csh
#

cd `dirname $0` && source ./Configuration

# set this variable to either 'sybase' or 'postgres'
setenv DB_TYPE                  sybase

umask 002

setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start daily public sybase reports >>& ${LOG}

cd ${PUBDAILY}

foreach i (*.py)
    echo `date`: $i >>& ${LOG}
    $i >>& ${LOG}
end

cd ${REPORTOUTPUTDIR}

echo `date`: Copy reports >>& ${LOG}
foreach i (gene_association.mgi)
    echo `date`: $i >>& ${LOG}
    cp $i ${FTPREPORTDIR}
end

echo `date`: End nightly public sybase reports >>& ${LOG}
