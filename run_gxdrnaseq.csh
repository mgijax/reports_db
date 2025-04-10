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
# remove old files, tar, gzip
#
rm -rf ${GXDRNASEQDIR}/*

#
# Generate GXD RNA Seq public reports.
#
cd ${PUBGXDRNASEQ}
${PYTHON} GXD_RnaSeq.py >>& ${LOG}

#
# tar & gzip files
#
echo `date`: tar and gzip reports | tee -a ${LOG}
cd ${GXDRNASEQDIR}
foreach i (*.rpt)
rm -rf $i.gz
gzip $i
end
rm -rf gxdrnaseq.tar
tar -cvf gxdrnaseq.tar *.gz | tee -a ${LOG}
gzip gxdrnaseq.tar

echo `date`: End GXD RNA Seq public reports | tee -a ${LOG}

exit 0
