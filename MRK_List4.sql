print ""
print "Genetic Marker List (sorted by chromosome/excludes withdrawns)"
print ""

select chromosome "Chr", symbol "Symbol", substring(name, 1, 40) "Name"
from MRK_Marker
where _Species_key = 1
and chromosome like '[0-9]'
and _Marker_Status_key = 1
order by chromosome, symbol
go

select chromosome "Chr", symbol "Symbol", substring(name, 1, 40) "Name"
from MRK_Marker
where _Species_key = 1
and chromosome not like '[0-9]'
and _Marker_Status_key = 1
order by chromosome, symbol
go

