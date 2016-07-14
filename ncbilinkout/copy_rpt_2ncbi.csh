#!/bin/csh -f

#
# Script to copy NCBI LinkOut files to NCBI FTP site
#
# Usage: copy_rpt_2ncbi.csh
#
# Runs Monday morning with the public swap
#
# NCBI's LinkOut (http://www.ncbi.nlm.nih.gov/projects/linkout/)
# enables us to generate XML files that contain links from NCBI to MGI.
#
# NCBI uses their LinkOut software to take the XML files we provide
# and process them correctly.
#
# NCBI FTP host, etc.:
#        Host: ftp-private.ncbi.nlm.nih.gov
#        Username: mgd
#        Password: keUPdTuD
#
#
# History:
# lnh   - 07/14/2016
#       TR12036 

cd `dirname $0` && source ../Configuration

umask 002

#
# Initialize the log file.
#
setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}
#
# If this is production, FTP the files to the NCBI LinkOut host.
#
if ( ${INSTALL_TYPE} == "prod" ) then

echo `date`: Send data to NCBI LinkOut host | tee -a ${LOG}
cd ${REPORTOUTPUTDIR}
ftp -in ${NCBILINKOUT_HOST} <<END | tee -a ${LOG}
user ${NCBILINKOUT_USER} ${NCBILINKOUT_PASSWORD}
cd holdings
mput protein-mgd.xml
mput providerinfo.xml
mput pubmed-mgd-1.xml
mput nucleotide-mgd-*.xml
mput pubmed-mgd.uid
quit
END

endif

echo `date`: reports copied to NCBI LinkOut host | tee -a ${LOG}
