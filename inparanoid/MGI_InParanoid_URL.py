#!/usr/local/bin/python

'''
#
# MGI_InParanoid_URL.py
#
# Report:
#	TR 6500
#	Requested by InParanoid
#
#	Produce a one-line report of:
#
#	Marker Accession ID URL
#
# Usage:
#       MGI_InParanoid_URL.py
#
# History:
#
# lec   01/25/2005
#       - created
#
'''

import sys
import os
import string
import mgi_utils
import reportlib

#
# Main
#

reportName = 'Mus-musculus_MGI_' + mgi_utils.date('%m%d%Y') + '_protein'
fp = reportlib.init(reportName, outputdir = os.environ['INPARANOIDDIR'], printHeading = 0, fileExt = '-URL')
fp.write('http://www.informatics.jax.org/searches/accession_report.cgi?id=\n')
reportlib.finish_nonps(fp)

