#!/bin/csh -f

#
# runone_postgres_report.csh
#
# Script to run a single, specified public postgres report. 
#
# Usage: runone_postgres_report.csh <directory>/<report script name>
# Example: runone_postgres_report.csh weekly_postgres/HOM_MouseHumanSequence.py
#
# lnh: 12/2014

set USAGE = "Usage: $0 <directory>/<report script name>"

cd `dirname $0` && source ./Configuration
# check that we received an acceptable argument

if ($#argv != 1) then
        echo "Error: Wrong number of arguments"
        echo $USAGE
        exit 1
else if (! -e $1) then
        echo "Error: File $1 does not exist"
        echo $USAGE
        exit 1
else if (! -r $1) then
        echo "Error: File $1 is not readable"
        echo $USAGE
        exit 1
endif
# set this variable to either 'sybase' or 'postgres'
setenv DB_TYPE                 postgres

umask 002

setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Generating $1 public postgres report >>& ${LOG}
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
#
# Generate weekly reports
#
echo "/run_in_parallel.py -p ${SUBPROCESSES} $1 " >>& ${LOG}
./run_in_parallel.py -p ${SUBPROCESSES} $1 >>& ${LOG}
echo `date`: $1 public postgres report generated >>& ${LOG}

exit 0
