print ""
print "Genetic Marker List (sorted by chromosome/includes withdrawns)"
print ""

select a.accID "MGI Accession ID", m.chromosome "Chr", m.symbol "Symbol", substring(m.name, 1, 40) "Name"
from MRK_Marker m, MRK_Chromosome c, MRK_Acc_View a
where m._Species_key = 1
and m._Species_key = c._Species_key
and m.chromosome = c.chromosome
and m._Marker_key = a._Object_key
and a.prefixPart = "MGI:"
and a.preferred = 1
order by c.sequenceNum, m.symbol
go

