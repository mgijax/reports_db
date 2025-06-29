#!/bin/csh -f

#
# Script to generate and copy NCBI LinkOut files to NCBI FTP site
#
# Usage: ncbilinkout.csh
#
# NCBI's LinkOut (https://www.ncbi.nlm.nih.gov/projects/linkout/)
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

umask 002

#
# Initialize the log file.
#
setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start NCBI LinkOut | tee -a ${LOG}

foreach i (*.py)
    echo `date`: $i | tee -a ${LOG}
    ${PYTHON} $i >>& ${LOG}
end

# for ncbi
cp providerinfo.xml ${REPORTOUTPUTDIR}
# for E PMC
cp pmc_providerinfo.xml ${REPORTOUTPUTDIR}
echo `date`: End NCBI LinkOut | tee -a ${LOG}
