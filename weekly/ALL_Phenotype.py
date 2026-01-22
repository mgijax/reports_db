
'''
# ALL_Phenotype.py
#
# Report:
#       Tab-delimited file of alleles and their associated phenotypes 
#       (as determined by _annotype_key = 1028).
#       WTS2-1792
#
# Usage:
#       ALL_Phenotype.py
#

# Output Format:
#
#   1. Allele ID
#   2. Allele Symbol
#   3. Phenotypes. A pipe-separated list of MPid,term pairs.
#
'''

import sys
import os
import reportlib
import db

db.setTrace()

CRT = reportlib.CRT
TAB = reportlib.TAB

fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

fp.write('Allele ID' + TAB)
fp.write('Allele Symbol' + TAB)
fp.write('Phenotypes' + CRT)

results = db.sql('''
    with allphenos as (
    select aa.accid as allele_id, a.symbol, concat_ws(',', aa2.accid, t.term) as phenotype
    from voc_annot va, all_allele a, acc_accession aa, voc_term t, acc_accession aa2
    where va._annottype_key = 1028
    and va._object_key = a._allele_key
    and a._allele_key = aa._object_key
    and aa._mgitype_key = 11
    and aa._logicaldb_key = 1
    and aa.preferred = 1
    and va._term_key = t._term_key
    and t._term_key = aa2._object_key
    and aa2._mgitype_key = 13
    and aa2.preferred = 1
    )
    select allele_id, symbol, STRING_AGG(distinct phenotype, '|' ORDER BY phenotype) as phenotypes
    from allphenos
    group by allele_id, symbol
    order by symbol
    ''', 'auto')

for r in results:	
    allele_id = r['allele_id']
    symbol = r['symbol']
    phenotypes = r['phenotypes']
    line = TAB.join([allele_id,symbol,phenotypes]) + CRT
    fp.write(line)

reportlib.finish_nonps(fp)
