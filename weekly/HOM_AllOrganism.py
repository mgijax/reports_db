"""
HOM_AllOrganism.py

# Report:
# 
# Produce a homology class report for all organism that includes the following
# tab-delimited fields:
#
# 1)  2/22/21 Replaced HomoloGene ID with database key
# 2)  Common Organism Name
# 3)  NCBI Taxon ID
# 4)  Symbol
# 5)  EntrezGene ID
# 6)  Mouse MGI ID
# 7)  HGNC ID
# 8)  OMIM Gene ID
# 9)  Genetic Location
# 10) Genome Coordinates
# 11) Name
# 12) Synonyms (pipe-delimited)
#
# Usage:
#       HOM_AllOrganism.py
#
# HISTORY
#
# sc    02/22/2021
#       TR13349 - B39 project. Update to use alliance direct homology
#               (was using Homologene and HGNC)
#
"""

import sys
import os
import string
import reportlib
import db

db.setTrace()

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PAGE = reportlib.PAGE

#
# Lookups
#
mgiID = {}
hgncID = {}
omimID = {}
genLoc = {}
genCoord = {}
synonym = {}


#
# Define the sort order of the organisms within a homology class
# based on organism key.
#
# sc 2/22/21 - we only have mouse, human, rat, zebrafish, but leaving
# the others in as the Alliance will be adding organisms to their homology
#
organismOrder = [ 1, 2, 10, 94, 13, 11, 40, 63, 95, 84 ]


#
# Purpose: This function will write one record to the report that represents
#          one results set from the database.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def writeRecord (r):
    markerKey = r['_Marker_key']
    fp.write(str(r['_Cluster_key']) + TAB)
    fp.write(r['commonName'] + TAB)
    fp.write(r['taxonID'] + TAB)
    fp.write(r['symbol'] + TAB)
    fp.write(r['entrezgeneID'] + TAB)
    if markerKey in mgiID:
        fp.write(mgiID[markerKey] + TAB)
    else:
        fp.write(TAB)
    if markerKey in hgncID:
        fp.write(hgncID[markerKey] + TAB)
    else:
        fp.write(TAB)
    if markerKey in omimID:
        fp.write(omimID[markerKey] + TAB)
    else:
        fp.write(TAB)
    if markerKey in genLoc:
        fp.write(genLoc[markerKey] + TAB)
    else:
        fp.write(TAB)
    if markerKey in genCoord:
        fp.write(genCoord[markerKey] + TAB)
    else:
        fp.write(TAB)
    fp.write(r['name'] + TAB)
    if markerKey in synonym:
        fp.write(synonym[markerKey] + CRT)
    else:
        fp.write(CRT)
    return


#
# Purpose: This function will write one homology class to the report. It takes
#          a dictionary as the only argument where the key is an organism key
#          and the value is a list of one or more results sets from the
#          database. Each result set is a list of dictionaries, as returned
#          by db.sql().
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def writeHomologyClass (homology):

    #
    # Process the results in the homology class dictionary in the defined
    # organism sort order.
    #
    for organismKey in organismOrder:

        #
        # If the homology class dictionary has an entry for the current
        # organism, pull out the list of results sets that need to be
        # written to the report.
        #
        if organismKey in homology:
            homologyResults = homology[organismKey]

            #
            # Write each results set to the report.
            #
            for r in homologyResults:
                writeRecord(r)

            #
            # Remove the list of results sets from the homology class
            # dictionary for the current organism.
            #
            homology.pop(organismKey)

    #
    # When we are done processing the organisms that are define in the sort
    # order list, check to see if there are any additional results sets in
    # the homology class dictionary for other organisms that have not been
    # defined in the sort order list. If so, write those to the report.
    #
    for organismKey in homology.keys():
        homologyResults = homology[organismKey]
        for r in homologyResults:
            writeRecord(r)

    return


#
# Main
#

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# Get the mouse and human build versions to display with the genomic location
# in the header record. 
#
results = db.sql('''
        select m._Organism_key, lc.version
        from MRK_Marker m, MRK_Location_Cache lc
        where m._Organism_key in (1,2)
        and m._Marker_key = lc._Marker_key
        and lc.version is not null
        ''','auto')

mouseVersion = ''
humanVersion = ''
for r in results:
    if r['_Organism_key'] == 1:
        mouseVersion =  r['version']
    if r['_Organism_key'] == 2:
        humanVersion =  r['version']


#
# Print the header record for the report.
#
fp.write('DB Class Key' + TAB)
fp.write('Common Organism Name' + TAB)
fp.write('NCBI Taxon ID' + TAB)
fp.write('Symbol' + TAB)
fp.write('EntrezGene ID' + TAB)
fp.write('Mouse MGI ID' + TAB)
fp.write('HGNC ID' + TAB)
fp.write('OMIM Gene ID' + TAB)
fp.write('Genetic Location' + TAB)
fp.write('Genome Coordinates (mouse: ' + mouseVersion + \
                ' human: ' + humanVersion + ')' + TAB)
fp.write('Name' + TAB)
fp.write('Synonyms' + CRT)


#
# Get Cluster Key, organism name, taxon ID, marker symbol/name. and
# EntrezGene ID that will be needed for each row of the report. These are
# saved in a temp table because the marker key for each row will be used
# to query for the remaining info.
#
db.sql('''
        select c._Cluster_key,
               o._Organism_key,
               o.commonName,
               a2.accID as taxonID,
               m._Marker_key,
               m.symbol,
               m.name,
               a3.accID as entrezgeneID
        into temporary table temp1
        from  MRK_Cluster c,
             MRK_ClusterMember cm,
             MRK_Marker m,
             MGI_Organism o,
             ACC_Accession a2,
             ACC_Accession a3
        where c._ClusterType_key = 9272150 and
              c._ClusterSource_key = 75885739 and
              c._Cluster_key = cm._Cluster_key and
              cm._Marker_key = m._Marker_key and
              m._Organism_key = o._Organism_key and
              o._Organism_key = a2._Object_key and
              a2._MGIType_key = 20 and
              a2._LogicalDB_key = 32 and
              cm._Marker_key = a3._Object_key and
              a3._MGIType_key = 2 and
              a3._LogicalDB_key = 55
        ''', None)

db.sql('create index idx1 on temp1 (_Marker_key)', None)


#
# Build a lookup for MGI IDs for each marker key.
#
results = db.sql('''
        select t._Marker_key,
               a.accID
        from temp1 t,
             ACC_Accession a
        where t._Marker_key = a._Object_key and
              a._MGIType_key = 2 and
              a._LogicalDB_key = 1 and
              a.prefixPart = 'MGI:' and
              a.private = 0 and
              a.preferred = 1
        ''','auto')

for r in results:
    mgiID[r['_Marker_key']] = r['accID']


#
# Build a lookup for HGNC IDs for each marker key.
#
results = db.sql('''
        select t._Marker_key,
               a.accID
        from temp1 t,
             ACC_Accession a
        where t._Marker_key = a._Object_key and
              a._MGIType_key = 2 and
              a._LogicalDB_key = 64 and
              a.prefixPart = 'HGNC:' and
              a.private = 0 and
              a.preferred = 1
        ''','auto')

for r in results:
    hgncID[r['_Marker_key']] = r['accID']


#
# Build a lookup for OMIM Gene IDs for each marker key.
#
results = db.sql('''
        select t._Marker_key,
               a.accID
        from temp1 t,
             ACC_Accession a
        where t._Marker_key = a._Object_key and
              a._MGIType_key = 2 and
              a._LogicalDB_key = 15 and
              a.private = 0 and
              a.preferred = 1
        ''','auto')

for r in results:
    omimID[r['_Marker_key']] = r['accID']


#
# Build lookups for genetic location and genome coordinates for each marker key.
# NOTE: The cmoffset is converted to a str.with 2 digits to the right of the
#       decimal point. The start/end coordinates are converted to str.with
#       no decimal positions.
#
# For example: cmoffset 12.349999999999999 becomes 12.35
#              startCoordinate 1234567.0 becomes 1234567
#
results = db.sql('''
        select t._Marker_key,
               lc.chromosome,
               lc.cytogeneticOffset,
               substr(cast(lc.cmoffset as varchar),10,2) as cmoffset,
               lc.genomicChromosome,
               cast(lc.startCoordinate as varchar) as startCoordinate,
               cast(lc.endCoordinate as varchar) as endCoordinate,
               lc.strand
        from temp1 t,
             MRK_Location_Cache lc
        where t._Marker_key = lc._Marker_key and
              lc._Organism_key in (1,2)
        union
        select t._Marker_key,
               m.chromosome,
               m.cytogeneticOffset,
               NULL, NULL, NULL, NULL, NULL
        from temp1 t,
             MRK_Marker m
        where t._Marker_key = m._Marker_key and
              m._Organism_key not in (1,2)
        ''','auto')

for r in results:
    #
    # If there is a genomic chromosome, use it with the start/end coordinates
    # and the strand to build the genomic coordinate str. This should only
    # be applicable for mouse or human.
    #
    genomicChr = r['genomicChromosome']
    if genomicChr != None:

        #
        # The strand can be null, so build that piece separately.
        #
        if r['strand'] != None:
            strand = '(' + r['strand'] + ')'
        else:
            strand = ''

        #
        # Add the genomic coordinate str.to the lookup for the marker.
        #
        genCoord[r['_Marker_key']] = 'Chr' + genomicChr + ':' + \
                                     r['startCoordinate'].strip() + '-' + \
                                     r['endCoordinate'].strip() + strand

    #
    # Use the genetic chromosome with one of the cmoffsets (if available) to
    # build the genetic location str.
    #
    geneticChr = r['chromosome']
    cmoffset = r['cmoffset']
    cytogeneticOffset = r['cytogeneticOffset']

    #
    # If there is an cmoffset, use it to build the genetic location str.
    #
    if cmoffset != None:
        cmoffset = cmoffset.strip()
        if cmoffset == '-999.00':
            location = 'Chr' + geneticChr
        elif cmoffset == '-1.00':
            location = 'Chr' + geneticChr + ' syntenic'
        else:
            location = 'Chr' + geneticChr + ' ' + cmoffset + ' cM'

    #
    # If there is a cytogenetic cmoffset, use it to build the genetic location
    # str.
    #
    elif cytogeneticOffset != None:
        location = 'Chr' + geneticChr + ' ' + cytogeneticOffset

    #
    # Otherwise, just use the chromosome.
    #
    else:
        location = 'Chr' + geneticChr

    #
    # Add the genetic location str.to the lookup for the marker.
    #
    genLoc[r['_Marker_key']] = location


#
# Build a lookup for synonyms for each marker key.
#
results = db.sql('''
        select t._Marker_key,
               s.synonym
        from temp1 t,
             MGI_Synonym s
        where t._Marker_key = s._Object_key and
              s._MGIType_key = 2
        order by t.symbol,
                 s.synonym
        ''','auto')

synList = ''
markerKey = ''

for r in results:
    if r['_Marker_key'] != markerKey and markerKey != '':
        synonym[markerKey] = synList
        synList = ''

    markerKey = r['_Marker_key']
    if synList == '':
        synList = r['synonym']
    else:
        synList = synList + '|' + r['synonym']

synonym[markerKey] = synList


#
# Get the Cluster Key, organism and marker data from the temp table that
# will represent separate lines on the report.
#
results = db.sql('''
        select _Cluster_key,
               _Organism_key,
               commonName,
               taxonID,
               _Marker_key,
               symbol,
               name,
               entrezgeneID
        from temp1 
        order by _Cluster_key,_Organism_key,
                 symbol
        ''','auto')

h = {}
clusterKey = -1

for r in results:
    #
    # If this row is the start of a new homology class and it is not the
    # first one, write out the prior one and reset the homology dictionary.
    #
    if r['_Cluster_key'] != clusterKey and clusterKey != -1:
        writeHomologyClass(h)
        h = {}

    clusterKey = r['_Cluster_key']

    #
    # If there is already a list of results sets in the homology dictionary
    # for the current organism, pull the list and add the new one to it.
    # Otherwise, start a new list.
    #
    organismKey = r['_Organism_key']
    if organismKey in h:
        list = h[organismKey]
    else:
        list = []
    list.append(r)

    #
    # Put the list of results sets in the dictionary for the current
    # organism.
    #
    h[organismKey] = list

#
# Write out the last homology class to the report.
#
writeHomologyClass(h)

reportlib.finish_nonps(fp)
