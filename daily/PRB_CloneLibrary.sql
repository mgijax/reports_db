
print ""
print "MGI Clone Libraries"
print ""

select name
from PRB_Source
where name is not null
order by name
go

