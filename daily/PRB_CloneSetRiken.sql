
set nocount on
go

select distinct p._Probe_key, db = ldb.name
into #riken
from PRB_Probe p, ACC_Accession pa, ACC_LogicalDB ldb
where p._Probe_key = pa._Object_key
and pa._MGIType_key = 3
and pa._LogicalDB_key = 26
and pa._LogicalDB_key = ldb._LogicalDB_key
go

insert #riken
select distinct p._Probe_key, db = ldb.name
from PRB_Probe p, ACC_Accession pa, ACC_LogicalDB ldb
where p._Probe_key = pa._Object_key
and pa._MGIType_key = 3
and pa._LogicalDB_key = 51
and pa._LogicalDB_key = ldb._LogicalDB_key
and not exists (select 1 from #riken where p._Probe_key = _Probe_key)
go

select distinct p._Probe_key, pm._Marker_key, r.db
into #prbs
from #riken r, PRB_Source ps, PRB_Probe p, MGI_Set s, MGI_SetMember sm, 
PRB_Marker pm
where s._Set_key = sm._Set_key 
and s._MGIType_key = 5
and s.name = 'RIKEN'
and sm._Object_key = ps._Source_key
and ps._Source_key = p._Source_key
and p._Probe_key *= r._Probe_key
and p._Probe_key *= pm._Marker_key
go

set nocount off
go

print ""
print "MGI Clone Set: RIKEN"
print ""
print "Description: A row in this report represents a member of the Clone"
print "Set:'RIKEN'.  Note a clone in this report may also exist in one"
print "of the other Clone Set reports."
print ""

select segmentID = pa.accID, segmentName = p.name, markerID = ma.accID, 
markerSymbol = m.symbol, tmp.db
from #prbs tmp, PRB_Probe p, ACC_Accession pa, MRK_Marker m, ACC_Accession ma
where tmp._Probe_key = p._Probe_key 
and p._Probe_key = pa._Object_key 
and pa._MGIType_key = 3
and pa.preferred = 1
and pa.prefixPart = 'MGI:'
and tmp._Marker_key *= m._Marker_key 
and tmp._Marker_key *= ma._Object_key
and ma._MGIType_key = 2
and ma.preferred = 1
and ma.prefixPart = 'MGI:'
order by db, markerSymbol, pa.accID
go

