print ""
print "Genetic Marker List (sorted by chromosome/excludes withdrawns)"
print ""

select a.accID "MGI Accession ID", m.chromosome "Chr", 
"cM Position" =
        case
        when o.offset >= 0 then str(o.offset, 10, 2)
        when o.offset = -999.0 then "       N/A"
        when o.offset = -1.0 then "  syntenic"
        end
, m.symbol "Symbol", substring(m.name,1,150) "Name", substring(t.name,1,25) "Type"
from MRK_Marker m, MRK_Chromosome c, MRK_Acc_View a, MRK_Offset o, MRK_Types t
where m._Species_key = 1
and m._Marker_Status_key in (1,3)
and m._Species_key = c._Species_key
and m.chromosome = c.chromosome
and m._Marker_key = a._Object_key
and a.prefixPart = "MGI:"
and a._LogicalDB_key = 1
and a.preferred = 1
and m._Marker_key = o._Marker_key
and o.source = 0
and m._Marker_Type_key = t._Marker_Type_key
order by c.sequenceNum, m.symbol
go

