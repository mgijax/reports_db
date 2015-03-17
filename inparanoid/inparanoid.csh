#!/bin/csh -f

#
# inparanoid.csh
#
# Script to generate InParanoid files
#
# Usage: inparanoid.csh
#

cd `dirname $0` && source ../Configuration

umask 002

setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

# remove old files
rm -rf ${INPARANOIDDIR}/*

# unzip refseq protein fasta file
rm -rf ${REFSEQFASTA}
/usr/bin/gunzip -c ${REFSEQFASTAGZ} > ${REFSEQFASTA}

# derive a FASTA file from the UniProt file
./sp_to_fasta ${UNIPROTDAT} > ${UNIPROTFASTA}

# generate files
foreach i (`ls *.py`)
    $i >>& ${LOG}
end
