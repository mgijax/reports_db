#!/bin/csh -f

#
# inparanoid.csh
#
# Script to generate InParanoid files
#
# Usage: inparanoid.csh
#

cd `dirname $0` && source ../Configuration

#
# set postgres configuration to use the public instance
#
setenv PG_DBSERVER      ${PG_PUB_DBSERVER}
setenv PG_DBNAME        ${PG_PUB_DBNAME}
setenv PG_DBUSER        ${PG_PUB_DBUSER}

echo "PG_DBSERVER: ${PG_DBSERVER}"
echo "PG_DBNAME: ${PG_DBNAME}"
echo "PG_DBUSER: ${PG_DBUSER}"

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
${SP2FASTA} ${UNIPROTDAT} > ${UNIPROTFASTA}

# generate files
foreach i (`ls *.py`)
    $i >>& ${LOG}
end
