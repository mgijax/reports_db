#!/bin/csh -f

#
# run_gxdrnaseq.csh
#
# Run GXD Reg Seq reports
# Run from Jenkins because this takes 8+ hours
# See Jenkins task: Run GXD RNA Seq Report
#

cd `dirname $0` && source ./Configuration

umask 002

#
# Initialize the log file.
#
setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start GXD RNA Seq public reports | tee -a ${LOG}

#
# Generate GXD RNA Seq public reports.
#
cd ${PUBGXDRNASEQ}
${PYTHON} GXD_RnaSeq.py >>& ${LOG}

#
# gzip files
#
#cd ${GXDRNASEQDIR}
#echo `date`: tar and gzip reports | tee -a ${LOG}
#rm -rf gxdrnaseq.tar gxdrnaseq.tar.gz
#tar -cvf gxdrnaseq.tar GXD*.rpt
#gzip gxdrnaseq.tar

echo `date`: End GXD RNA Seq public reports | tee -a ${LOG}

exit 0
