print ""
print "Genetic Marker List (sorted by chromosome/excludes withdrawns)"
print ""

select a.accID "MGI Accession ID", m.chromosome "Chr", m.symbol "Symbol", m.name "Name"
from MRK_Marker m, MRK_Chromosome c, MRK_Acc_View a
where m._Species_key = 1
and m._Marker_Status_key = 1
and m._Species_key = c._Species_key
and m.chromosome = c.chromosome
and m._Marker_key = a._Object_key
and a.prefixPart = "MGI:"
and a.preferred = 1
order by c.sequenceNum, m.symbol
go

