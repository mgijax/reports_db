#!/bin/csh -f

#
# runone_sybase_report.csh
#
# Script to run a single, specified public sybase report. 
#
# Usage: runone_sybase_report.csh <directory>/<report script name>
#
# Example: runone_sybase_report.csh weekly_sybase/MRK_NomenUpdates.py
#

cd `dirname $0` && source ./Configuration

umask 002

set USAGE = "Usage: $0 <directory>/<report script name>"

#
# Check that we received an acceptable argument.
#
if ($#argv != 1) then
    echo "Error: Wrong number of arguments"
    echo $USAGE
    exit 1
else if (! -e $1) then
    echo "Error: File $1 does not exist"
    exit 1
else if (! -r $1) then
    echo "Error: File $1 is not readable"
    exit 1
endif

#
# Get the directory and report names from the argument.
#
set DIRNAME = `echo $1 | cut -d'/' -f1`
set RPTNAME = `echo $1 | cut -d'/' -f2`

#
# Set this variable to either 'sybase' or 'postgres'.
#
setenv DB_TYPE		sybase

#
# Initialize the log file.
#
setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start public sybase report | tee -a ${LOG}

#
# Generate weekly reports.
#
cd ${DIRNAME}
echo `date`: $1 | tee -a ${LOG}
${RPTNAME} >>& ${LOG}

echo `date`: End public sybase report | tee -a ${LOG}

exit 0
