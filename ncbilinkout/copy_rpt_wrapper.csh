#!/bin/csh -f
#
#  copy_rpt_wrapper.csh
###########################################################################
#
#  Purpose:
#
#      This script is a wrapper around the process that copies the NCBI
#      linkout reports to the NCBI FTP site.
#
#  Usage:
#
#      copy_rpt_wrapper.csh
#
#  Env Vars:
#
#      - See Configuration file (reports_db product)
#
#      - See master.config.csh (mgiconfig product)
#
#  Inputs:
#
#      - Process control flags
#
#  Outputs:
#
#      - Log file for the script (${LOG})
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  Fatal error occurred
#
#  Assumes:  Nothing
#
#  Implementation:
#
#      This script will perform following steps:
#
#      1) Source the configuration file to establish the environment.
#      2) Wait for the flag to signal that the public webshare has been
#         swapped. This indicates that the public release is proceeding.
#      3) Call the script that copies the reports.
#
#  Notes:  None
#
###########################################################################

cd `dirname $0` && source ../Configuration

setenv SCRIPT_NAME `basename $0`

setenv LOG ${REPORTLOGSDIR}/${SCRIPT_NAME}.log
rm -rf ${LOG}
touch ${LOG}

#
# Wait for the "Webshare Swapped" flag to be set. Stop waiting if the number
# of retries expires or the abort flag is found.
#
date | tee -a ${LOG}
echo 'Wait for the "Webshare Swapped" flag to be set' | tee -a ${LOG}

setenv RETRY ${PROC_CTRL_RETRIES}
while (${RETRY} > 0)
    setenv READY `${PROC_CTRL_CMD_PUB}/getFlag ${NS_PUB_LOAD} ${FLAG_WEBSHR_SWAPPED}`
    setenv ABORT `${PROC_CTRL_CMD_PUB}/getFlag ${NS_PUB_LOAD} ${FLAG_ABORT}`

    if (${READY} == 1 || ${ABORT} == 1) then
        break
    else
        sleep ${PROC_CTRL_WAIT_TIME}
    endif

    setenv RETRY `expr ${RETRY} - 1`
end

#
# Terminate the script if the number of retries expired or the abort flag
# was found.
#
if (${RETRY} == 0) then
   echo "${SCRIPT_NAME} timed out" | tee -a ${LOG}
   date | tee -a ${LOG}
   exit 1
else if (${ABORT} == 1) then
   echo "${SCRIPT_NAME} aborted by process controller" | tee -a ${LOG}
   date | tee -a ${LOG}
   exit 1
endif

#
# Call the script to copy the reports.
#
echo 'Call the script to copy the reports' | tee -a ${LOG}
${PUBRPTS}/ncbilinkout/copy_rpt.csh

echo "${SCRIPT_NAME} completed successfully" | tee -a ${LOG}
date | tee -a ${LOG}
exit 0
