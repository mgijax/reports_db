print ""
print "Genetic Marker List (sorted alphabetically/includes withdrawns)"
print ""

select a.accID "MGI Accession ID", m.chromosome "Chr", 
"cM Position" =
        case
        when o.offset >= 0 then str(o.offset, 10, 2)
        when o.offset = -999.0 then "       N/A"
        when o.offset = -1.0 then "  syntenic"
        end
,m.symbol "Symbol", upper(substring(s.status, 1, 1)) "Status", substring(m.name,1,150) "Name", substring(t.name,1,25) "Type"
from MRK_Marker m, ACC_Accession a, MRK_Offset o, MRK_Types t, MRK_Status s
where m._Organism_key = 1
and m._Marker_Status_key in (1,3)
and m._Marker_key = a._Object_key
and a._MGIType_key = 2
and a.prefixPart = "MGI:"
and a._LogicalDB_key = 1
and a.preferred = 1
and m._Marker_key = o._Marker_key
and o.source = 0
and m._Marker_Type_key = t._Marker_Type_key
and m._Marker_Status_key = s._Marker_Status_key
union
select null, m.chromosome,
"cM Position" =
        case
        when o.offset >= 0 then str(o.offset, 10, 2)
        when o.offset = -999.0 then "       N/A"
        when o.offset = -1.0 then "  syntenic"
        end
,m.symbol "Symbol", upper(substring(s.status, 1, 1)) "Status", substring(m.name,1,150) "Name", substring(t.name,1,25) "Type"
from MRK_Marker m, MRK_Offset o, MRK_Types t, MRK_Status s
where m._Organism_key = 1
and m._Marker_Status_key = 2
and m._Marker_key = o._Marker_key
and o.source = 0
and m._Marker_Type_key = t._Marker_Type_key
and m._Marker_Status_key = s._Marker_Status_key
order by m.symbol
go

