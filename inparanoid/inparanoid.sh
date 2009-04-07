#!/bin/csh

#
# inparanoid.sh
#
# Script to generate InParanoid files
#
# Usage: inparanoid.sh
#

cd `dirname $0` && source ../Configuration

umask 002

setenv LOG	${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

# remove old files
rm -rf ${INPARANOIDDIR}/*

# unzip refseq protein fasta file
rm -rf ${REFSEQFASTA}
/usr/bin/gunzip -c ${REFSEQFASTAGZ} > ${REFSEQFASTA}

# derive a FASTA file from the UniProt file
${SP2FASTA} ${UNIPROTDAT} > ${UNIPROTFASTA}

# generate files
foreach i (`ls *.py`)
$i >>& ${LOG}
end

