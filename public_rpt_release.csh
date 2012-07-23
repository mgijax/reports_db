#!/bin/csh -f
#
#  public_rpt_release.csh
###########################################################################
#
#  Purpose:
#
#      This script is a wrapper that controls the copying of the public
#      reports to the public FTP site.
#
#  Usage:
#
#      public_rpt_release.csh
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
#      - Process control flags
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
#      3) Copy the non-mouse gene association file to the public FTP site.
#      4) Copy the public reports to the public FTP site.
#      5) Copy the inparanoid files to the public FTP site.
#
#  Notes:  None
#
###########################################################################

cd `dirname $0` && source ./Configuration

setenv SCRIPT_NAME `basename $0`

setenv LOG ${REPORTLOGSDIR}/${SCRIPT_NAME}.log
rm -f ${LOG}
touch ${LOG}

echo "$0" | tee -a ${LOG}
env | sort | tee -a ${LOG}

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
# Copy the reports to the public FTP site.
#
cd ${REPORTOUTPUTDIR}

echo `date`: Copy non-mouse gene association file to ${FTPCUSTOM} | tee -a ${LOG}
echo `date`: gene_association.mgi_nonmouse | tee -a ${LOG}
cp gene_association.mgi_nonmouse ${FTPCUSTOM}
rm -f gene_association.mgi_nonmouse

echo `date`: Copy reports to ${FTPREPORTDIR} | tee -a ${LOG}
foreach i (*)
    if ( ! -d $i ) then
        echo `date`: $i | tee -a ${LOG}
        cp $i ${FTPREPORTDIR}
    endif
end
rm -rf ${FTPREPORTDIR}/iphone*[0-9]*

echo `date`: Copy inparanoid files to ${FTPCUSTOM}/inparanoid | tee -a ${LOG}
rm -f ${FTPCUSTOM}/inparanoid/Mus-musculus*
cd ${INPARANOIDDIR}
foreach i (Mus-musculus* aaseq*)
    echo `date`: $i | tee -a ${LOG}
    cp $i ${FTPCUSTOM}/inparanoid
end

echo "${SCRIPT_NAME} completed successfully" | tee -a ${LOG}
date | tee -a ${LOG}
exit 0
