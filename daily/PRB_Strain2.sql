set nocount on
go

select distinct r._Reference_key, s.strain, s._Strain_key
into #rflvs
from PRB_Probe p, PRB_Reference f, PRB_RFLV r, PRB_Allele a,
PRB_Allele_Strain al, PRB_Strain s
where p.DNAtype != 'primer'
and p._Probe_key = f._Probe_key
and f._Reference_key = r._Reference_key
and r._RFLV_key = a._RFLV_key
and a._Allele_key = al._Allele_key
and al._Strain_key = s._Strain_key
go
 
select distinct r._Reference_key, s.strain, s._Strain_key
into #pcrs
from PRB_Probe p, PRB_Reference f, PRB_RFLV r, PRB_Allele a,
PRB_Allele_Strain al, PRB_Strain s
where p.DNAtype = 'primer'
and p._Probe_key = f._Probe_key
and f._Reference_key = r._Reference_key
and r._RFLV_key = a._RFLV_key
and a._Allele_key = al._Allele_key
and al._Strain_key = s._Strain_key
go
 
set nocount off
go

print ""
print "This list represents genetic backgrounds for which polymporphism data have"
print "been curated. These are not standard strains, but mixtures of strains, crosses, "
print "common names, and species designations. See the Standard Strains for a list of "
print "strains with accepted nomenclature."
print ""

print ""
print "PCR Polymorphism Strains and # of References - Sorted by # of References" 
print ""

select distinct strain, num = count(*)
from #pcrs
group by _Strain_key
order by num desc, strain
go

print ""
print "RFLP Polymorphism Strains and # of References - Sorted by # of References" 
print ""

select distinct strain, num = count(*)
from #rflvs
group by _Strain_key
order by num desc, strain
go

