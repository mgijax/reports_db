print ""
print "Withdrawn Mutant Allele Markers (sorted by chromosome)"
print ""

select m2.chromosome "Chr ", m2.symbol "Marker Symbol", m1.symbol "Mutant Symbol",
a.symbol "Allele", substring(a.name, 1, 30) "Allele Name"
from MRK_Current c, MRK_Marker m1, MRK_Marker m2, ALL_Allele a, MRK_Chromosome h
where m1._Marker_Status_key = 2
and m1.name like '%allele of%'
and m1._Marker_key = c._Marker_key
and c._Current_key = m2._Marker_key
and m2._Marker_key = a._Marker_key
and a.symbol like m2.symbol + "%<" + m1.symbol + "%>"
and m2._Species_key = h._Species_key
and m2.chromosome = h.chromosome
order by h.sequenceNum, m2.symbol
go

