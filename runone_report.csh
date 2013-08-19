#!/bin/csh -f

#
# runone_report.csh
#
# Script to run a single, specified public report.  Useful for debugging
# without running an entire suite (weekly, daily, etc.) of reports.
#
# Usage: runone_report.csh <directory>/<report script name>
# Example: runone_report.csh weekly/MRK_List.py
#

set USAGE = "$0 <directory>/<report script name>"

cd `dirname $0` && source ./Configuration

#
# set postgres configuration to use the public instance
#
if $?PG_PUB_DBSERVER then
    setenv PG_DBSERVER      ${PG_PUB_DBSERVER}
endif
if $?PG_PUB_DBNAME then
    setenv PG_DBNAME        ${PG_PUB_DBNAME}
endif
if $?PG_PUB_DBUSER then
    setenv PG_DBUSER        ${PG_PUB_DBUSER}
endif
echo $MGD_DBSERVER $MGD_DBNAME

# check that we received an acceptable argument

if ($#argv != 1) then
	echo $USAGE
	echo "Error: Wrong number of arguments"
	exit 1
else if (! -e $1) then
	echo $USAGE
	echo "Error: File $1 does not exist"
	exit 1
else if (! -r $1) then
	echo $USAGE
	echo "Error: File $1 is not readable"
	exit 1
endif

umask 002

setenv LOG ${REPORTLOGSDIR}/`basename $1`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start run of $1 | tee -a ${LOG}

#
# Generate the report
#
echo `date`: $1 | tee -a ${LOG}
cd `dirname $1`
`basename $1` >>& ${LOG}
echo Report output dir: ${REPORTOUTPUTDIR}

echo `date`: End run of $1 | tee -a ${LOG}
