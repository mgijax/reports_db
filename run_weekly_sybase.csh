#!/bin/csh -f

#
# run_weekly_sybase.csh
#
# Script to generate weekly public sybase reports.
#
# Usage: run_weekly_sybase.csh
#

cd `dirname $0` && source ./Configuration

umask 002

setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start weekly public sybase reports | tee -a ${LOG}

#
# Generate weekly reports
#
./run_in_parallel.py -p ${SUBPROCESSES} weekly_sybase/*py inparanoid/inparanoid.csh >>& ${LOG}

#cd ${PUBWEEKLY_SYBASE}
#foreach i (*.py)
#    echo `date`: $i | tee -a ${LOG}
#    $i >>& ${LOG}
#end

echo `date`: End weekly public sybase reports | tee -a ${LOG}
