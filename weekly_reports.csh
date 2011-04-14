#!/bin/csh

#
# weekly_reports.csh
#
# Script to generate weekly public reports.
#
# Usage: weekly_reports.csh
#

cd `dirname $0` && source ./Configuration

umask 002

setenv LOG ${REPORTLOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

echo `date`: Start weekly public reports | tee -a ${LOG}

cd ${PUBWEEKLY}

foreach i (*.sql)
    echo `date`: $i | tee -a ${LOG}
#    reportisql.csh $i ${REPORTOUTPUTDIR}/$i.rpt ${MGD_DBSERVER} ${MGD_DBNAME} MGI >> ${LOG}
end

foreach i (*.py)
    if ( $i != "MGI_CloneSet.py" ) then
        echo `date`: $i | tee -a ${LOG}
        $i >>& ${LOG}
    endif
end

#
# Generate clone set reports.
#
foreach i ("Image" "NIA 15K,NIA 7.4K,NIA" "RIKEN (FANTOM),RIKEN" "RPCI-23" "RPCI-24")
    echo `date`: MGI_CloneSet.py $i | tee -a ${LOG}
    ./MGI_CloneSet.py "$i" >>& ${LOG}
end

#
# Generate inparanoid files.
#
echo `date`: inparanoid.csh | tee -a ${LOG}
${PUBRPTS}/inparanoid/inparanoid.csh >>& ${LOG}

#
# Generate NCBI LinkOut files.
#
echo `date`: ncbilinkout.csh | tee -a ${LOG}
${PUBRPTS}/ncbilinkout/ncbilinkout.csh >>& ${LOG}

echo `date`: End weekly public reports | tee -a ${LOG}
