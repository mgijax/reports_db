
print ""
print "MGI cDNA Clone Libraries"
print ""

select name = substring(ps.name, 1, 100), cloneCollection = s.name
from MGI_Set s, MGI_SetMember sm, PRB_Source ps
where s._MGIType_key = 5
and s.name in ('IMAGE', 'NIA', 'RIKEN', 'RIKEN (FANTOM)')
and s._Set_key = sm._Set_key
and sm._Object_key = ps._Source_key
order by name
go

