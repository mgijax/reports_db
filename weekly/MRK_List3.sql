set nocount on
go

select m._Marker_key, m._Marker_Status_key, m._Marker_Type_key, 
	m.symbol, name = substring(name,1,150), m.chromosome, c.sequenceNum
into #markers
from MRK_Marker m, MRK_Chromosome c
where m._Organism_key = 1
and m._Organism_key = c._Organism_key
and m.chromosome = c.chromosome
go

create index idx1 on #markers(_Marker_key)
go

select m.*, 
markerStatus = upper(substring(s.status, 1, 1)),
markerType = substring(t.name,1,25),
cmPosition =
        case
        when o.offset >= 0 then str(o.offset, 10, 2)
        when o.offset = -999.0 then "       N/A"
        when o.offset = -1.0 then "  syntenic"
        end
into #markersAll
from #markers m, MRK_Status s, MRK_Types t, MRK_Offset o
where m._Marker_key = o._Marker_key
and o.source = 0
and m._Marker_Type_key = t._Marker_Type_key
and m._Marker_Status_key = s._Marker_Status_key
go

create index idx1 on #markersAll(_Marker_key)
create index idx2 on #markersAll(symbol)
create index idx3 on #markersAll(sequenceNum)
go

select a.accID, m.chromosome, m.cmPosition, m.symbol, m.markerStatus, m.name, m.markerType, m.sequenceNum
into #markersPrint
from #markersAll m, ACC_Accession a
where m._Marker_Status_key in (1, 3)
and m._Marker_key = a._Object_key
and a._MGIType_key = 2
and a.prefixPart = "MGI:"
and a._LogicalDB_key = 1
and a.preferred = 1
union
select null, m.chromosome, m.cmPosition, m.symbol, m.markerStatus, m.name, m.markerType, m.sequenceNum
from #markersAll m
where m._Marker_Status_key = 2
order by m.sequenceNum, m.symbol
go

print ""
print "Genetic Marker List (sorted by chromosome/includes withdrawns)"
print ""

select accID "MGI Accession ID", chromosome "Chr", cmPosition "cM Position", 
	symbol "Symbol", markerStatus "Status", name "Name", markerType "Type"
from #markersPrint
go

