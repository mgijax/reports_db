print ""
print "Genetic Marker List (sorted by chromosome/includes withdrawns)"
print ""

select chromosome "Chr", symbol "Symbol", substring(name, 1, 40) "Name"
from MRK_Marker
where _Species_key = 1
and chromosome like '[0-9]'
order by chromosome, symbol
go

select chromosome "Chr", symbol "Symbol", substring(name, 1, 40) "Name"
from MRK_Marker
where _Species_key = 1
and chromosome not like '[0-9]'
order by chromosome, symbol
go

