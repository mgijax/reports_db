print ""
print "Genetic Marker List (sorted alphabetically/includes withdrawns)"
print ""

select a.accID "MGI Accession ID", m.chromosome "Chr", m.symbol "Symbol", substring(m.name, 1, 40) "Name"
from MRK_Marker m, MRK_Acc_View a
where m._Species_key = 1
and m._Marker_key = a._Object_key
and a.prefixPart = "MGI:"
and a.preferred = 1
order by m.symbol
go

