
set nocount on
go

select ps._Source_key
into #source
from PRB_Source ps, MGI_Set s, MGI_SetMember sm
where s.name = 'IMAGE'
and s._Set_key = sm._Set_key
and sm._Object_key = ps._Source_key
and ps._Organism_key = 1
go

create index idx_source_key on #source(_Source_key)
go

select p._Probe_key, pm._Marker_key
into #prbs
from #source s, PRB_Probe p, PRB_Marker pm
where s._Source_key = p._Source_key
and p._Probe_key = pm._Probe_key
union
select p._Probe_key, null
from #source s, PRB_Probe p
where s._Source_key = p._Source_key
and not exists (select 1 from PRB_Marker pm where p._Probe_key = pm._Probe_key)
go

create index idx_prb_key on #prbs(_Probe_key)
go

select segmentID = pa.accID, segmentName = p.name, tmp._Marker_key
into #final
from #prbs tmp, PRB_Probe p, ACC_Accession pa
where tmp._Probe_key = p._Probe_key 
and p._Probe_key = pa._Object_key 
and pa._MGIType_key = 3
and pa._LogicalDB_key = 1
and pa.prefixPart = 'MGI:'
and pa.preferred = 1
go

create index idx_mrk_key on #final(_Marker_key)
go

set nocount off
go

print ""
print "MGI Clone Set: IMAGE"
print ""
print "Description: A row in this report represents a mouse member of the "
print "Clone Set:'IMAGE'.  Note a clone in this report may also exist in one"
print "of the other Clone Set reports."
print ""

select f.segmentID, f.segmentName, markerID = ma.accID, markerSymbol = m.symbol
from #final f, MRK_Marker m, ACC_Accession ma
where f._Marker_key is not null
and f._Marker_key = m._Marker_key 
and f._Marker_key = ma._Object_key
and ma._MGIType_key = 2
and ma._LogicalDB_key = 1
and ma.prefixPart = 'MGI:'
and ma.preferred = 1
union
select f.segmentID, f.segmentName, null, null
from #final f
where f._Marker_key is null
order by markerSymbol, f.segmentID
go

