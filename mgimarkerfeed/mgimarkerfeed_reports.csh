#!/bin/csh -f

#
# mgimarkerfeed_reports.csh
#
# Script to generate mgi marker feed reports.
#
# Usage: mgimarkerfeed_reports.csh
#

setenv TARFILE mgimarkerfeed.tar

cd `dirname $0` && source ../Configuration

umask 002

setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start MGI marker feed reports | tee -a ${LOG}
rm -rf ${REPORTOUTPUTDIR}/mgimarkerfeed/* >>& ${LOG}

echo `date`: MarkerFeed.py | tee -a ${LOG}
${PYTHON} mgiMarkerFeed.py >>& ${LOG}

echo `date`: Create tar file | tee -a ${LOG}
cd ${REPORTOUTPUTDIR}/mgimarkerfeed
tar cvf ${TARFILE} *.bcp >>& ${LOG}

echo `date`: Compress tar file | tee -a ${LOG}
gzip -f ${TARFILE} >>& ${LOG}

echo `date`: Copy file to FTP site | tee -a ${LOG}
cp ${TARFILE}.gz ${FTPREPORTDIR}

echo `date`: End MGI marker feed reports | tee -a ${LOG}
