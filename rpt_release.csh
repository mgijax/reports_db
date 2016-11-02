#!/bin/csh -f
#
#  rpt_release.csh
###########################################################################
#
#  Purpose:
#
#      This script copies the public reports to the public FTP site.
#
#  Usage:
#
#      rpt_release.csh
#
#  Env Vars:
#
#      - See Configuration file (reports_db product)
#
#      - See master.config.csh (mgiconfig product)
#
#  Inputs:  None
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
#      2) Copy the reports.
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

echo `date`: Copy CvDC files to ${FTPREPORTDIR}/cvdc | tee -a ${LOG}
rm -f ${FTPREPORTDIR}/cvdc/*html
cd ${CVDCDIR}
foreach i (*)
    echo `date`: $i | tee -a ${LOG}
    cp $i ${FTPREPORTDIR}/cvdc
end

echo "${SCRIPT_NAME} completed successfully" | tee -a ${LOG}
date | tee -a ${LOG}
exit 0
