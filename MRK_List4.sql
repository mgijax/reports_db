print ""
print "Genetic Marker List (sorted by chromosome/excludes withdrawns)"
print ""

select m.chromosome "Chr", m.symbol "Symbol", substring(m.name, 1, 40) "Name"
from MRK_Marker m, MRK_Chromosome c
where m._Species_key = 1
and m._Marker_Status_key = 1
and m._Species_key = c._Species_key
and m.chromosome = c.chromosome
order by c.sequenceNum, m.symbol
go

