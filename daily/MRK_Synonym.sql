print ""
print "Marker Symbols and their Synonyms (ordered by Marker Symbol)"
print ""

select a.accID "MGI Accession ID", m.chromosome "Chr", "cM Position" =
        case
        when o.offset >= 0 then str(o.offset, 10, 2)
        when o.offset = -999.0 then "       N/A"
        when o.offset = -1.0 then "  syntenic"
        end,
"Marker Symbol" = m.symbol, "Synonym" = substring(s.name,1,90)
from MRK_Marker m, ACC_Accession a, MRK_Other s, MRK_Offset o 
where m._Organism_key = 1 
and m._Marker_key = a._Object_key 
and a._MGIType_key = 2
and a.prefixPart = "MGI:" 
and a._LogicalDB_key = 1
and a.preferred = 1 
and m._Marker_key = o._Marker_key 
and o.source = 0 
and m._Marker_key = s._Marker_key
order by m.symbol
go
