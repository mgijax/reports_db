print ""
print "InterPro Domains"
print ""

select a.accID "ID", t.term "InterPro Domain Name"
from VOC_Vocab v, VOC_Term t, VOC_Term_Acc_View a
where v.name = "InterPro Domains"
and v._Vocab_key = t._Vocab_key
and t._Term_key = a._Object_key
order by t.sequenceNum
go

