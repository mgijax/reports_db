#!/usr/local/bin/python

'''
#
# HMD_HumanPhenotype.py 
#
# Report:
#     A tab-delimited report of all mouse-human orthologous gene pairs, 
#     with all PhenoSlim IDs associated with any allele of the mouse gene.
#
#     Output fields: 
#     
#     Human Gene Symbol 
#     Human LocusLink ID 
#     Mouse Gene Symbol 
#     MGI Accession ID 
#     PhenoSlim IDs associated with any allele of the mouse gene (space- or comma-separated) 
#
#     Sorted by:
#     Human Gene Symbol. 
#    
# Usage:
#       HMD_HumanPhenotype.py
#
# Used by:
#
# Notes:
#
# History:
#
# lec	08/27/2003
#	- TR 5082; make a public report
#
# dm    07/14/03
#       -Inital Creation
# 
'''
 
import sys
import os
import db
import reportlib
import string

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#############
#           #
# Functions #
#           #
#############
'''
# requires: results,    The array of data
#           keyField,   The name of the key to use to 
#                       create the dictionary entries
#           valueField, The name of the row field to use
#                       as the value for the dictionary entry
# effects:
#     Creates a dictionary from a result array from db.sql
#     using the keyfield as the key and 
#     the value field as the key value.
#
# returns:
#     a dictionary
'''
def createDict(results, keyField, valueField):

        d = {}
        for r in results:
                key = r[keyField]

                value = r[valueField]
                if not d.has_key(key):
                        d[key] = ' '
                d[key] = string.join([d[key],value,])
        return d

########
#      #
# MAIN #
#      #
########

cmds = []

#################################################################
# Get mouse to human orthologous marker pair's symbols and keys #
# and the Phenoslim IDs associated to the mouse markers and     #
# place them in a temp table to be used by other queries        #
# needed to produce the desired report                          #
#################################################################
cmds.append(string.join(['select distinct(hpv.markerkey2) as humanKey, ',
                         '                 m2.symbol as humanSym,      ',
                         '                hpv.markerkey1 as mouseKey,  ',
                         '                 m1.symbol as mouseSym       ',
                         'into #homology                               ',
                         'from HMD_Homology_Pairs_View hpv,            ',
                         '     MRK_Marker m1,                          ',
                         '     MRK_Marker m2                           ',
                         'where hpv.specieskey1    = 1                 ',
                         '  and hpv.specieskey2    = 2                 ',
                         '  and  m1._Marker_key    = hpv.markerkey1    ',
                         '  and  m2._Marker_key    = hpv.markerkey2    ',]))

######################################################
# Create some indexes to speed up the report queries #
######################################################
cmds.append('create nonclustered index index_mouseKey on #homology(mouseKey)')
cmds.append('create nonclustered index index_humanKey on #homology(humanKey)')
cmds.append('create nonclustered index index_humanSym on #homology(humanSym)')

#################################################################
# Get the mouse to human orthologous marker pairs from the temp #
# table created in the first query so that it can be used as    #
# data for creating the desired report                          #
#################################################################
cmds.append('select * from #homology order by humanSym')

###################################################################
# Get the MGI IDs for the Mouse markers in the temp table created #
# by the first query                                              #
###################################################################
cmds.append(string.join(['select distinct(a._Object_key),   ',
                         '       a.accID                    ',
                         'from #homology h,                 ',
                         '     ACC_Accession a              ',
                         'where a._Object_key  = h.mouseKey ',
                         '  and a._MGIType_key = 2          ',
                         '  and a._LogicalDB_key = 1        ',
                         '  and a.prefixPart   = "MGI:"     ',
                         '  and a.preferred    = 1          ',]))

#####################################################################
# Get the LocusLink for the Human markers in the temp table created # 
# by the first query                                                #
#####################################################################
cmds.append(string.join(['select distinct(a._Object_key),     ',
                         '       a.accID                      ',
                         'from #homology h,                   ',
                         '     ACC_Accession a                ',
                         'where a._Object_key    = h.humanKey ',
                         '  and a._MGIType_key   = 2          ',
                         '  and a._LogicalDB_key = 24         ',]))

#####################################################################
# Get the Phenoslim for the Mouse markers in the temp table created # 
# by the first query                                                #
#####################################################################
cmds.append(string.join(['select h.mouseKey,                       ',
                         '       a.accID                           ',
                         'from #homology h,                        ',
                         '     GXD_AlleleGenotype g,               ',
                         '     VOC_Annot v,                        ',
                         '     ACC_Accession a                     ',
                         'where g._Marker_key    = h.mouseKey      ',
                         '  and v._Object_key    = g._Genotype_key ',
                         '  and v._AnnotType_key = 1001            ',
                         '  and a._Object_key    = v._Term_key     ',
                         '  and a._MGIType_key   = 13              ',
                         '  and a._LogicalDB_key = 34              ',]))

results = db.sql(cmds, 'auto')

########################################################
# Create keyed dictionaries of the results so they can #
# be called by there associated _Marker_key.           #
########################################################
mmgi   = createDict(results[5], '_Object_key', 'accID')
hlocus = createDict(results[6], '_Object_key', 'accID')
mpheno = createDict(results[7], 'mouseKey', 'accID')

####################
# Open report file #
####################
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = 0)

############################
# Process db query results #
############################
for r in results[4]:
    ################################
    # Output the human gene symbol #
    ################################
    fp.write(r['humanSym'] + TAB)

    #####################################################
    # If the human marker has a locus link ID output it #
    #####################################################
    if hlocus.has_key(r['humanKey']):
        fp.write(hlocus[r['humanKey']] + TAB)
    else:
        fp.write(TAB)

    ################################
    # Output the mouse gene symbol #
    ################################
    fp.write(r['mouseSym'] + TAB)

    #########################################################
    # If the mouse marker has an MGI accession ID output it #
    #########################################################
    if mmgi.has_key(r['mouseKey']):
        fp.write(mmgi[r['mouseKey']] + TAB)

    ##################################################
    # The mouse marker has Phenoslim IDs output them #
    ##################################################
    if mpheno.has_key(r['mouseKey']):
        fp.write(mpheno[r['mouseKey']] + TAB)
    else:
        fp.write(TAB)

    fp.write(CRT)

#########################
# Close the output file #
#########################
reportlib.finish_nonps(fp)

