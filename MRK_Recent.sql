print ""
print "Genetic Marker List Updates Within the Last Month (sorted alphabetically)"
print ""

select l.symbol "Symbol", l.chromosome "Chr", substring(l.name, 1, 40) "Name"
from MRK_Mouse_View l
where datepart(year, modification_date) = datepart(year, getdate())
and datepart(month, modification_date) >= datepart(month, getdate()) - 1
order by symbol
go

