
set nocount on
go

select distinct p._Probe_key, db = ldb.name
into #nia
from PRB_Probe p, PRB_ACC_View pa, ACC_LogicalDB ldb
where p._Probe_key = pa._Object_key
and pa._LogicalDB_key = 46
and pa._LogicalDB_key = ldb._LogicalDB_key
go

insert #nia
select distinct p._Probe_key, db = ldb.name
from PRB_Probe p, PRB_ACC_View pa, ACC_LogicalDB ldb
where p._Probe_key = pa._Object_key
and pa._LogicalDB_key = 50
and pa._LogicalDB_key = ldb._LogicalDB_key
and not exists (select 1 from #nia where p._Probe_key = _Probe_key)
go

insert #nia
select distinct p._Probe_key, db = ldb.name
from PRB_Probe p, PRB_ACC_View pa, ACC_LogicalDB ldb
where p._Probe_key = pa._Object_key
and pa._LogicalDB_key = 49
and pa._LogicalDB_key = ldb._LogicalDB_key
and not exists (select 1 from #nia where p._Probe_key = _Probe_key)
go

select distinct p._Probe_key, pm._Marker_key, n.db
into #prbs
from #nia n, PRB_Source ps, PRB_Probe p, MGI_Set s, MGI_SetMember sm, 
PRB_Marker pm
where s._Set_key = sm._Set_key 
and s._MGIType_key = 5
and s.name = 'NIA'
and sm._Object_key = ps._Source_key
and ps._Source_key = p._Source_key
and p._Probe_key *= n._Probe_key
and p._Probe_key *= pm._Marker_key
go

set nocount off
go

print ""
print "MGI Clone Set: NIA"
print ""
print "Description: A row in this report represents a member of the Clone"
print "Set:'NIA'.  Note a clone in this report may also exist in one"
print "of the other Clone Set reports."
print ""

select segmentID = pa.accID, segmentName = p.name, markerID = ma.accID, 
markerSymbol = m.symbol, tmp.db
from #prbs tmp, PRB_Probe p, PRB_ACC_View pa, MRK_Marker m, MRK_ACC_View ma
where tmp._Probe_key = p._Probe_key 
and p._Probe_key = pa._Object_key 
and pa.preferred = 1
and pa.prefixPart = 'MGI:'
and tmp._Marker_key *= m._Marker_key 
and tmp._Marker_key *= ma._Object_key
and ma.preferred = 1
and ma.prefixPart = 'MGI:'
order by db, markerSymbol, segmentName
go

