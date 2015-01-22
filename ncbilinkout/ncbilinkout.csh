#!/bin/csh -f

#
# Script to generate and copy NCBI LinkOut files to NCBI FTP site
#
# Usage: ncbilinkout.csh
#
# NCBI's LinkOut (http://www.ncbi.nlm.nih.gov/projects/linkout/)
# enables us to generate XML files that contain links from NCBI to MGI.
#
# NCBI uses their LinkOut software to take the XML files we provide
# and process them correctly.
#
# Our job:
#    * provide an identity file (providerinfo.xml) to NCBI
#    * provide one XML resource file for each LinkOut site
#      (PubMed, Protein, Nucleotide, OMIM) we are interested in providing
#    * put the XML files on the NCBI FTP site on a monthly basis
#      in order to keep our links up-to-date
#    * NCBI FTP host, etc.:
#        Host: ftp-private.ncbi.nlm.nih.gov
#        Username: mgd
#        Password: keUPdTuD
#
# NCBI's job:
#    * to pick up the XML files we send them and "make it so" (do the LinkOuts)
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

echo 'begin NCBI LinkOut', `date` | tee -a ${LOG}

foreach i (*.py)
    echo | tee -a ${LOG}
    echo $i, `date` | tee -a ${LOG}
    $i >>& ${LOG}
    echo $i, `date` | tee -a ${LOG}
    echo | tee -a ${LOG}
end

cp providerinfo.xml ${REPORTOUTPUTDIR}

# if running this from production...

if ( ${INSTALL_TYPE} == "prod" ) then

echo 'sending data to NCBI LinkOut host...' | tee -a ${LOG}
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

echo 'finished NCBI LinkOut', `date` | tee -a ${LOG}
