print ""
print "Genetic Marker List (sorted alphabetically/excludes withdrawns)"
print ""

select chromosome "Chr", symbol "Symbol", substring(name, 1, 40) "Name"
from MRK_Marker
where _Species_key = 1
and _Marker_Status_key = 1
order by symbol
go

