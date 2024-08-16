'''
#
# Report:
#
# Post on the  MGI Data and Statistical Reports
# area: Vertebrate Homology
# report #: 5
# report title: Mouse Protein Coding Genes having one-to-one Orthology with Human Genes
# 
# Use the same logic as for the MGI Stat referenced above.
# For protein-coding mouse/human orthologs with one-to-one correspondence, report:
# 
# Fields
# MGI Marker Accession ID
# Mouse Gene Symbol
# Mouse NCBI Gene ID
# HGNC ID
# Human Gene Symbol
# Human NCBI Gene ID
# 
# order by Mouse Gene Symbol
# update weekly
#
# Usage:
#       HOM_ProteinCoding.py
#
# History:
#
# lec   08/16/2024
#       - fl2-948/Create a public report for MGI stat: mouse protein-coding genes with 1-to-1 homology to human
#
'''
 
import sys 
import os
import reportlib
import db

db.setTrace()

CRT = reportlib.CRT
TAB = reportlib.TAB

#
# Main
#

fp = reportlib.init(sys.argv[0], printHeading = None, outputdir = os.environ['REPORTOUTPUTDIR'])

results = db.sql('''
        with all_mouse_human_pairs as (
         /* all pairs of mouse/human genes from Alliance Direct clusters */
         select
           m._marker_key as mouse_key,
           m.symbol as mouse_symbol,
           m2._marker_key as human_key,
           m2.symbol as human_symbol
         from  MRK_Marker m, MRK_ClusterMember mcm, MRK_Cluster mc,
               MRK_Marker m2, MRK_ClusterMember mcm2
         where mc._clusterType_key = 9272150
           and mc._clusterSource_key = 75885739
           and mc._cluster_key = mcm._cluster_key
           and mcm._marker_key = m._marker_key
           and m._organism_key = 1
           and m._marker_status_key in (1,3)
           and mc._cluster_key = mcm2._cluster_key
           and mcm2._marker_key = m2._marker_key
           and m2._organism_key = 2
        ),
        mouse_1 as (
         /* Mouse genes with exactly 1 human gene */
         select mouse_key, count(*)
         from all_mouse_human_pairs
         group by mouse_key
         having count(*) = 1
        ),
        human_1 as (
         /* Human genes with exactly 1 mouse gene */
         select human_key, count(*)
         from all_mouse_human_pairs
         group by human_key
         having count(*) = 1
        ),
        mouse_human_1_1 as (
         /* The pairs where both genes are in the _1 tables. */
         select amh.*
         from all_mouse_human_pairs amh, mouse_1 m1, human_1 h1
         where amh.mouse_key = m1.mouse_key
           and amh.human_key = h1.human_key
        ),
        mouse_protein_coding_1_1 as (
         /* Select for mouse protein coding genes */
         select mh11.*, ma.accid as mouse_accid, mncbi.accid as mouse_ncbi, hncbi.accid as human_ncbi, hgnc.accid as hgnc_accid
         from MRK_MCV_Cache mcv, mouse_human_1_1 mh11, 
            ACC_Accession ma, ACC_Accession mncbi, ACC_Accession hncbi, ACC_Accession hgnc
         where mcv.term = 'protein coding gene'
           and mcv.qualifier = 'D'
           and mcv._Marker_key = mh11.mouse_key
           and mcv._Marker_key = ma._Object_key
           and ma._MGIType_key = 2
           and ma._LogicalDB_key = 1
           and ma.preferred = 1
           and mcv._Marker_key = mncbi._Object_key
           and mncbi._MGIType_key = 2
           and mncbi._LogicalDB_key = 55
           and mncbi.preferred = 1
           and mh11.human_key = hncbi._Object_key
           and hncbi._MGIType_key = 2
           and hncbi._LogicalDB_key = 55
           and hncbi.preferred = 1
           and mh11.human_key = hgnc._Object_key
           and hgnc._MGIType_key = 2
           and hgnc._LogicalDB_key = 64
           and hgnc.preferred = 1
        )
        /* And the total is... */
        select *
        from mouse_protein_coding_1_1
        order by mouse_symbol
        ''', 'auto')

for r in results:
        fp.write(r['mouse_accid'] + TAB)
        fp.write(r['mouse_symbol'] + TAB)
        fp.write(r['mouse_ncbi'] + TAB)
        fp.write(r['hgnc_accid'] + TAB)
        fp.write(r['human_symbol'] + TAB)
        fp.write(r['human_ncbi'] + CRT)

reportlib.finish_nonps(fp)      # non-postscript file
