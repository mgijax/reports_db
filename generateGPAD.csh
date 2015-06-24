#!/bin/csh -f

#
# generateGPAD.csh  
#
# Script to generate GPAD/GPI files.
#
# Usage: generateGPAD.csh 
#

cd `dirname $0` && source ./Configuration

umask 002

#
# Initialize the log file.
#
setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

#
# After generating the gaf file from public reports generation step,
# We need to add step to generate the gpad and gpi files
#
#
cd ${REPORTOUTPUTDIR}

echo `date`: Copy reports | tee -a ${LOG}
foreach i (gene_association.mgi)
    echo `date`: $i | tee -a ${LOG}
    echo "Generating GPAD/GPI files for GAF file $i"
    ${GAF_FPROCESSOR}/gaf2gpad.py --gaf=$i
    cp $i.* ${FTPREPORTDIR}
end

echo `date`: End GPAD/GPI generation | tee -a ${LOG}

exit 0
