#!/bin/csh

#
# inparanoid.sh
#
# Script to generate InParanoid files
#
# Usage: inparanoid.sh
#

cd `dirname $0` && source ./Configuration

umask 002

setenv LOG	${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

# remove old files
rm -rf ${INPARANOIDDIR}/Mus-musculus_MGI*
rm -rf ${INPARANOIDFASTA}

# derive a FASTA file from the SwissProt load file
/usr/local/wu-blast2.0/sp2fasta ${INPARANOIDDAT} > ${INPARANOIDFASTA}

cd inparanoid
foreach i (`ls *.py`)
$i >>& ${LOG}
end

