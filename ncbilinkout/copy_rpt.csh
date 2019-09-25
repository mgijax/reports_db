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
#
# sc	- 09/13/19 - TR13082
#    	updated to add PMC ftp of PubMed linkout files
#
# lnh   - 07/14/2016
#       TR12036 

cd `dirname $0` && source ../Configuration

setenv SCRIPT_NAME `basename $0`

umask 002

setenv LOG ${REPORTLOGSDIR}/${SCRIPT_NAME}.log
rm -f ${LOG}
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

echo `date`: Send data to E_PMC LinkOut host | tee -a ${LOG}
cd ${REPORTOUTPUTDIR}
ftp -in ${E_PMC_HOST} <<END | tee -a ${LOG}
user ${E_PMC_USER} ${E_PMC_PASSWORD}
cd ol9dybmp
mput pmc_providerinfo.xml
mput pubmed-mgd-1.xml
mput pubmed-mgd.uid
quit
END

endif

echo "${SCRIPT_NAME} completed successfully" | tee -a ${LOG}
date | tee -a ${LOG}
exit 0
