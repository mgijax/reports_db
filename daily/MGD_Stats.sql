set nocount on
go

print ""
print "MGI General Statistics"

print ""
print "Number of References"
select count(*) from BIB_Refs
go


print ""
print "Number of Genetic Markers"
select count(*) from MRK_Marker
where _Organism_key = 1 
and _Marker_Status_key in (1,3)
go

/* mutants are genes that have an allele with the same nomenclature */
/* ("foo" and "foo") AND do not have a GenBank, RefSeq, SwissProt or TrEMBL association */

select distinct m._Marker_key 
into #mutants 
from MRK_Marker m, ALL_ALlele a 
where m._Organism_key = 1 
and m._Marker_Type_key = 1 
and m._Marker_key = a._Marker_key 
and m.symbol = a.symbol 
and not exists (select 1 from ACC_Accession a 
where m._Marker_key = a._Object_key 
and a._MGIType_key = 2 
and a._LogicalDB_Key in (9, 13, 27, 41)) 
go

create index idx1 on #mutants(_Marker_key)
go

print ""
print "Number of Genes"
select count(*) from MRK_Marker m1
where m1._Organism_key = 1 
and m1._Marker_Type_key = 1
and m1._Marker_Status_key in (1,3)
and not exists (select 1 from #mutants m2 where m1._Marker_key = m2._Marker_key)
go

select _Marker_key into #markers from MRK_Marker
where _Organism_key = 1 
and _Marker_Type_key = 1
and _Marker_Status_key in (1,3)
go

create index idx1 on #markers(_Marker_key)
go

print ""
print "Number of Genes w/ DNA Sequence information"
select count(distinct a._Object_key)
from #markers m, ACC_Accession a
where m._Marker_key = a._Object_key
and a._MGIType_key = 2
and a._LogicalDB_Key = 9
go

print ""
print "Number of Genes w/ Protein Sequence information"
select count(distinct a._Object_key)
from #markers m, ACC_Accession a
where m._Marker_key = a._Object_key
and a._MGIType_key = 2
and a._LogicalDB_Key = 13
go

print ""
print "Number of Genetic Markers Mapped"
select count(*) from MRK_Marker
where _Organism_key = 1 
and _Marker_Status_key in (1,3)
and chromosome != 'UN'
go

print ""
print "Number of Genes Mapped"
select count(*) from MRK_Marker
where _Organism_key = 1 
and _Marker_Type_key = 1
and _Marker_Status_key in (1,3)
and chromosome != 'UN'
go

select distinct h1 = m1._Marker_key, h2 = m2._Marker_key
into #homology
from HMD_Homology r1, HMD_Homology_Marker h1,
HMD_Homology r2, HMD_Homology_Marker h2, MRK_Marker m1, MRK_Marker m2
where m1._Organism_key = 2
and m1._Marker_key = h1._Marker_key
and h1._Homology_key = r1._Homology_key
and r1._Class_key = r2._Class_key
and r2._Homology_key = h2._Homology_key
and h2._Marker_key = m2._Marker_key
and m2._Organism_key = 1
go

print ""
print "Number of Orthologies Mouse/Human" 
select count(*)
from #homology
go

print ""
print "Genes with Molecular Probes and Segments Data"
select count(distinct m._Marker_key)
from #markers m, PRB_Marker p, PRB_Probe pb
where m._Marker_key = p._Marker_key
and p._Probe_key = pb._Probe_key
and pb._SegmentType_key != 63473
go

print ""
print "Number of Genetic Markers with Molecular Polymorphisms" 
select count(distinct m._Marker_key)
from #markers m, PRB_RFLV p
where m._Marker_key = p._Marker_key
go

print ""
print "Number of Genes with Molecular Polymorphisms" 
select count(distinct m._Marker_key)
from MRK_Marker m, PRB_RFLV p
where m._Organism_key = 1
and m._Marker_Type_key = 1
and m._Marker_Status_key in (1,3)
and m._Marker_key = p._Marker_key
go

