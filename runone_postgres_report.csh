#!/bin/csh -f

#
# runone_postgres_report.csh
#
# Script to run a single, specified public postgres report. 
#
# Usage: runone_postgres_report.csh <directory>/<report script name>
#
# Example: runone_postgres_report.csh weekly_postgres/HOM_MouseHumanSequence.py
#
# History:
#
# sc:    03/21/20 python 3 upgrade
#
# lnh: 12/2014

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
# Initialize the log file.
#
setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start public postgres report | tee -a ${LOG}

#
# Generate weekly reports.
#
cd ${DIRNAME}
echo `date`: $1 | tee -a ${LOG}
${PYTHON} ${RPTNAME} >>& ${LOG}

echo `date`: End public postgres report | tee -a ${LOG}

exit 0
