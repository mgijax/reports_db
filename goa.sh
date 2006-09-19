#!/bin/csh

#
# goa.sh
#
# Script to generate GOA association file 
# To be appended to the GO_gene_association.py output file
#
# Usage: goa.sh
#

cd `dirname $0` && source ./Configuration

umask 002

setenv LOG	${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

# remove old files
rm -rf ${GOAINPUTFILE1}

# copy new file from /data/downloads and unzip
cd ${GOADIR}
cp ${GOAINPUTFILEGZ1} .
gunzip ${GOAINPUTFILE1}

# parse GOA file
cd ${GOAWORKDIR}
./goa.py

