
set nocount on
go

select distinct p._Probe_key, pm._Marker_key
into #prbs
from PRB_Source ps, PRB_Probe p, MGI_Set s, MGI_SetMember sm, PRB_Marker pm
where s._Set_key = sm._Set_key 
and s._MGIType_key = 5
and s.name = 'RPCI-23'
and sm._Object_key = ps._Source_key
and ps._Source_key = p._Source_key
and p._Probe_key *= pm._Probe_key
go

set nocount off
go

print ""
print "MGI Clone Set: RPCI-23"
print ""
print "Description: A row in this report represents a member of the Clone"
print "Set:'RPCI-23'.  Note a clone in this report may also exist in one"
print "of the other Clone Set reports."
print ""

select segmentID = pa.accID, segmentName = p.name, markerID = ma.accID, 
markerSymbol = m.symbol
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
order by markerSymbol, pa.accID 
go

