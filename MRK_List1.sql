print ""
print "Genetic Marker List (sorted alphabetically/includes withdrawns)"
print ""

select chromosome "Chr", symbol "Symbol", substring(name, 1, 40) "Name"
from MRK_Marker
where _Species_key = 1
order by symbol
go

