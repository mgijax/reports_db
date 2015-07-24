#!/usr/local/bin/python

'''
#
# Report:
#       
# Usage:
#       MGI_ElsevierDOI.py
#
# Used by:
# 	TR10980/Implement Elsevier DOI image linking in mgiws product
#
#
# - a simple list of doi ids that are associated with references from >= 1 of the following:
# - markers
# - probes
# - mapping experiments
# - GXD literature entries
# - expresion results
#   - structures with expresion results
#   - expression assays
# - alleles
# - sequences
#
# History:
#
# lec	02/29/2012
#	- created
#
'''
 
import sys
import os
import reportlib
import db

db.setTrace()
db.setAutoTranslate(False)
db.setAutoTranslateBE()

CRT = reportlib.CRT

db.useOneConnection(1)
fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

#
# markers
#
db.sql('''
	select distinct r._Refs_key, a.accID
	into #refs
	from ACC_Accession a, MRK_Reference r
	where a._LogicalDB_key = 65
	and a._MGIType_key = 1
	and a._Object_key = r._Refs_key
	''', None)

db.sql('create index idx on #refs(_Refs_key)', None)

#
# molecular probes & segments
#
db.sql('''
	insert into #refs
	select distinct r._Refs_key, a.accID
	from ACC_Accession a, PRB_Reference r
	where a._LogicalDB_key = 65
	and a._MGIType_key = 1
	and a._Object_key = r._Refs_key
	and not exists (select 1 from #refs rr where r._Refs_key = rr._Refs_key)
	''', None)

#
# mapping
#
db.sql('''
	insert into #refs
	select distinct r._Refs_key, a.accID
	from ACC_Accession a, MLD_Expts r
	where a._LogicalDB_key = 65
	and a._MGIType_key = 1
	and a._Object_key = r._Refs_key
	and not exists (select 1 from #refs rr where r._Refs_key = rr._Refs_key)
	''', None)

#
# gxd index
#
db.sql('''
	insert into #refs
	select distinct r._Refs_key, a.accID
	from ACC_Accession a, GXD_Index r
	where a._LogicalDB_key = 65
	and a._MGIType_key = 1
	and a._Object_key = r._Refs_key
	and not exists (select 1 from #refs rr where r._Refs_key = rr._Refs_key)
	''', None)

#
# gxd expression
#
db.sql('''
	insert into #refs
	select distinct r._Refs_key, a.accID
	from ACC_Accession a, GXD_Expression r
	where a._LogicalDB_key = 65
	and a._MGIType_key = 1
	and a._Object_key = r._Refs_key
	and r.isForGXD = 1
	and not exists (select 1 from #refs rr where r._Refs_key = rr._Refs_key)
	''', None)

#
# alleles
#
db.sql('''
	insert into #refs
	select distinct r._Refs_key, a.accID
	from ACC_Accession a, MGI_Reference_Assoc r, ALL_Allele aa
	where a._LogicalDB_key = 65
	and a._MGIType_key = 1
	and a._Object_key = r._Refs_key
	and r._MGIType_key = 11
	and r._Object_key = aa._Allele_key
	and aa.isWildType = 0
	and not exists (select 1 from #refs rr where r._Refs_key = rr._Refs_key)
	''', None)

#
# sequences
#
db.sql('''
	insert into #refs
	select distinct r._Refs_key, a.accID
	from ACC_Accession a, MGI_Reference_Assoc r
	where a._LogicalDB_key = 65
	and a._MGIType_key = 1
	and a._Object_key = r._Refs_key
	and r._MGIType_key = 19
	and not exists (select 1 from #refs rr where r._Refs_key = rr._Refs_key)
	''', None)


#
# GO annotations
#
db.sql('''
        insert into #refs
	select distinct ve._Refs_key, a.accID
	from ACC_Accession a, VOC_Evidence ve, VOC_Annot va
	where a._LogicalDB_key = 65
        and a._MGIType_key = 1
        and a._Object_key = ve._Refs_key
	and ve._Annot_key = va._Annot_key
	and va._AnnotType_key = 1000
        and not exists (select 1 from #refs rr where ve._Refs_key = rr._Refs_key)
        ''', None)

results = db.sql('select * from #refs', 'auto')
for r in results:
    fp.write(r['accID'] + CRT)

reportlib.finish_nonps(fp)
db.useOneConnection(0)
