print ""
print "Homology - Human vs. Mouse (Sorted by Mouse Chromosome)"
print ""

select distinct 
m2.chromosome "Mouse-Chr", 

"Mouse-cM Position" =
        case
        when o.offset >= 0 then str(o.offset, 10, 2)
        when o.offset = -999.0 then "       N/A"
        when o.offset = -1.0 then "  syntenic"
        end

,m2.symbol "Mouse-Symbol", 
m1.chromosome + m1.cytogeneticOffset "Human-Chr", 
m1.symbol "Human-Symbol",
substring(m2.name, 1, 75) "Name"

from HMD_Homology r1, HMD_Homology_Marker h1,
HMD_Homology r2, HMD_Homology_Marker h2,
MRK_Marker m1, MRK_Marker m2, MRK_Offset o, MRK_Chromosome c

where m1._Species_key = 2
and m1._Marker_key = h1._Marker_key
and h1._Homology_key = r1._Homology_key
and r1._Class_key = r2._Class_key
and r2._Homology_key = h2._Homology_key
and h2._Marker_key = m2._Marker_key
and m2._Species_key = 1
and m2._Marker_key = o._Marker_key
and o.source = 0
and m2._Species_key = c._Species_key
and m2.chromosome = c.chromosome
order by c.sequenceNum, o.offset
go

