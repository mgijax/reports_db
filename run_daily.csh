#!/bin/csh -f

#
# run_daily.csh
#
# Script to generate daily public reports.
#
# Usage: run_daily.csh
#

cd `dirname $0` && source ./Configuration

umask 002

#
# Initialize the log file.
#
setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start daily public reports | tee -a ${LOG}

#
# Generate daily public reports.
#
cd ${PUBDAILY}
#foreach i (*.py)
foreach i (GO_gene_association.py GO_gpi.py)
    echo `date`: $i | tee -a ${LOG}
    ${PYTHON} $i >>& ${LOG}
end

#
# Copy reports to ftp site
#
cd ${REPORTOUTPUTDIR}

echo `date`: Copy reports | tee -a ${LOG}
foreach i (gene_association.mgi gene_association_pro.mgi mgi.gpa mgi.gpi)
    echo `date`: $i | tee -a ${LOG}
    cat $i | gzip -cf9 > $i.gz
    touch $i $i.gz
    cp -p $i $i.gz ${FTPREPORTDIR}
end

echo `date`: End daily public reports | tee -a ${LOG}

exit 0
