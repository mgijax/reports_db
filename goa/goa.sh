#!/bin/csh

#
# goa.sh
#
# Script to:
#
#  1.  Take the downloaded GO association file
#  2.  Generate 2 output files
#
#	a.  goa.annot
#	    File of GO annotations that we can load into MGI
#
#	b.  goa.mgi
#	    File of GO annotations that we cannot load into MGI
#	    This file is appended to the gene_association.mgi
#	    file that gets generated from MGI.
#
# Usage: goa.sh
#

cd `dirname $0` && source ./goa.config

umask 002

setenv LOG	${GOALOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

# remove old files
rm -rf ${GOAINPUTFILE}

# copy new file from /data/downloads and unzip
cd ${GOAINPUTDIR}
cp ${GOAINPUTFILEGZ} .
gunzip ${GOAINPUTFILE}
rm -rf ${GOAINPUTFILESRT}
# important to sort the file so we can collapse "duplicate" lines
sort -k2,7 ${GOAINPUTFILE} > ${GOAINPUTFILESRT}

# parse GOA file
cd ${GOAOUTPUTDIR}
./goa_mgi.py

# load the annotations into MGI
${ANNOTLOAD}/annotload.csh ${PUBRPTS}/goa/goa.config

