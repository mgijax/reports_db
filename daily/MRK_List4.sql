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
, m.symbol "Symbol", upper(substring(s.status, 1, 1)) "Status", substring(m.name,1,150) "Name", substring(t.name,1,25) "Type"
from MRK_Marker m, MRK_Chromosome c, ACC_Accession a, MRK_Offset o, MRK_Types t, MRK_Status s
where m._Organism_key = 1
and m._Marker_Status_key in (1,3)
and m._Organism_key = c._Organism_key
and m.chromosome = c.chromosome
and m._Marker_key = a._Object_key
and a._MGIType_key = 2
and a.prefixPart = "MGI:"
and a._LogicalDB_key = 1
and a.preferred = 1
and m._Marker_key = o._Marker_key
and o.source = 0
and m._Marker_Type_key = t._Marker_Type_key
and m._Marker_Status_key = s._Marker_Status_key
order by c.sequenceNum, m.symbol
go

