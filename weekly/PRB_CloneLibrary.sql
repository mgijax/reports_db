
print ""
print "MGI Clone Libraries"
print ""

select name = substring(name, 1, 100)
from PRB_Source
where name is not null
order by name
go

