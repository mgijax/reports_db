
set nocount on
go

select distinct p._Probe_key, db = ldb.name
into #riken1
from PRB_Probe p, ACC_Accession pa, ACC_LogicalDB ldb
where p._Probe_key = pa._Object_key
and pa._MGIType_key = 3
and pa._LogicalDB_key = 26
and pa._LogicalDB_key = ldb._LogicalDB_key
go

create index idx1 on #riken1(_Probe_key)
go

select distinct p._Probe_key, db = ldb.name
into #riken2
from PRB_Probe p, ACC_Accession pa, ACC_LogicalDB ldb
where p._Probe_key = pa._Object_key
and pa._MGIType_key = 3
and pa._LogicalDB_key = 51
and pa._LogicalDB_key = ldb._LogicalDB_key
and not exists (select 1 from #riken1 where p._Probe_key = _Probe_key)
go

select _Probe_key, db
into #riken
from #riken1
union
select _Probe_key, db
from #riken2
go

create index idx1 on #riken(_Probe_key)
go

select sm._Object_key
into #rikensource
from MGI_Set s, MGI_SetMember sm
where s._MGIType_key = 5
and s.name = 'RIKEN'
and s._Set_key = sm._Set_key 
go

create index idx1 on #rikensource(_Object_key)
go

select p._Probe_key, segmentName = p.name
into #prbs1
from #rikensource s, PRB_Probe p
where s._Object_key = p._Source_key
go

create index idx1 on #prbs1(_Probe_key)
go

select p._Probe_key, pm._Marker_key, p.segmentName, r.db
into #prbs2
from #prbs1 p, PRB_Marker pm, #riken r
where p._Probe_key *= pm._Probe_key
and p._Probe_key *= r._Probe_key
go

create index idx1 on #prbs2(_Probe_key)
go
create index idx2 on #prbs2(_Marker_key)
go

select segmentID = pa.accID, p.segmentName, p.db, p._Marker_key
into #prbs
from #prbs2 p, ACC_Accession pa
where p._Probe_key = pa._Object_key 
and pa._MGIType_key = 3
and pa.prefixPart = 'MGI:'
and pa.preferred = 1
go

create index idx1 on #prbs(_Marker_key)
go

set nocount off
go

print ""
print "MGI Clone Set: RIKEN"
print ""
print "Description: A row in this report represents a member of the Clone Set:  'RIKEN'."
print "Note a clone in this report may also exist in one of the other Clone Set reports."
print ""

select p.segmentID, p.segmentName, markerID = ma.accID, markerSymbol = m.symbol, p.db
from #prbs p, MRK_Marker m, ACC_Accession ma
where p._Marker_key *= m._Marker_key 
and p._Marker_key *= ma._Object_key
and ma._MGIType_key = 2
and ma.preferred = 1
and ma.prefixPart = 'MGI:'
order by p.db, markerSymbol, p.segmentID
go

