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
foreach i (*.py)
    echo `date`: $i | tee -a ${LOG}
    $i >>& ${LOG}
end

#
# After generating the gaf file from step above,
# We need to add step to generate the gpad and gpi files
#
#
cd ${REPORTOUTPUTDIR}

echo `date`: Copy reports | tee -a ${LOG}
foreach i (gene_association.mgi)
    echo `date`: $i | tee -a ${LOG}
    echo "Generating GPAD/GPI files for GAF file $i"
    ${GAF_FPROCESSOR}/gaf2gpad.py --gaf=$i
    cp $i* ${FTPREPORTDIR}
end
foreach i (GO_eco_association.rpt)
    echo `date`: $i | tee -a ${LOG}
    cp $i ${FTPREPORTDIR}
end

echo `date`: End daily public reports | tee -a ${LOG}

exit 0
