select a.accID "ID", t.term "InterPro Domain Name"
from VOC_Vocab v, VOC_Term t, ACC_Accession a 
where v.name = "InterPro Domains" 
and v._Vocab_key = t._Vocab_key 
and t._Term_key = a._Object_key 
and a._MGIType_key = 13
order by t.sequenceNum
go
