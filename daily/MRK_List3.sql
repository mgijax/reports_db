set nocount on
go

select c.sequenceNum, a.accID, m.chromosome, o.offset, m.symbol, m.name
into #output
from MRK_Marker m, MRK_Chromosome c, MRK_Acc_View a, MRK_Offset o
where m._Species_key = 1
and m._Species_key = c._Species_key
and m.chromosome = c.chromosome
and m._Marker_key = a._Object_key
and a.prefixPart = "MGI:"
and a.preferred = 1
and m._Marker_key = o._Marker_key
and o.source = 0
union
select c.sequenceNum, null, m.chromosome, o.offset, m.symbol, m.name
from MRK_Marker m, MRK_Chromosome c, MRK_Offset o
where m._Species_key = 1
and m._Marker_Status_key = 2
and m._Species_key = c._Species_key
and m.chromosome = c.chromosome
and m._Marker_key = o._Marker_key
order by c.sequenceNum, m.symbol
go

set nocount off
go

print ""
print "Genetic Marker List (sorted by chromosome/includes withdrawns)"
print ""

select o.accID "MGI Accession ID", o.chromosome "Chr", 
"cM Position" =
        case
        when o.offset >= 0 then str(o.offset, 10, 2)
        when o.offset = -999.0 then "       N/A"
        when o.offset = -1.0 then "  syntenic"
        end
, o.symbol "Symbol", o.name "Name"
from #output o
go

