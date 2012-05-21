#!/usr/local/bin/python
"""
# tr9848.py
# 
# ORIGINAL REQUEST:
# Need a file of the following 7 fields (I added #8)
# 
# 1. MGI ID of the transgene or cre knock-in allele
# 2. Official name of line (get line name from IMSR w/ MGI correct name)
# 3. Promoter (driver) for transgene or gene for knock-in
# 4. Pubmed ID
# 5. Inducible or not
# 6. Holding site - either (i) repository name and ideally an ID so 
#    we can link or (ii) contact email of original publishing author
# 7. Anatomical site(s) of expression/specificity.
# 8. High level term (what we will use for the Cre release)
#
# We will need to generate this download regularly for the CREATE project.
# Thanks Janan 
# 
# USAGE:
#       tr9848.py <server> <database>
#
# HISTORY:
#
# lec   02/08/2012
#	- add 'allele['mname'] != None' to parseAllele
#	  because marker may be null
#
# lec	12/28/2011
#	- changed non-ansi-standard query to left outer join
#
# 2009-09-17	Joel Richardson
#	- created
# 2010-02-22	Joel Richardson
#	- modified according to updates to the TR. The new list of fields
#	is as follows (copied from TR):
#
#   1.MGI ID of the transgene allele or cre knock-in allele
#   *2.transgene or cre knock-in allele symbol
#   *3.transgene or cre knock-in allele name -- this should be either the transgene name 
#    (for transgenes) or a concatenation of the gene and allele name for knock-in alleles 
#    (example: Chattm1(cre)Lowl should have a name of: 
#	choline acetyltransferase; targeted mutation 1, Bradford B Lowell
#   4.Official name of line (get line name from IMSR w/ MGI correct name)
#   5.Promoter (driver) for transgene or gene for knock-in
#   * 6.Pubmed ID. For those that do not have a Pubmed ID, use "Journal, Author". Examples: 
#         MGI Direct Data Submission, G Schutz 
#         Unpublished, URL: http://jaxmice.jax.org/
#   7.Inducible or not
#   8.Holding site - abbreviation for repository name
#   *9.Repository strain ID or link (used by IMSR)
#   10.Anatomical site(s) of expression/specificity (specific terms).
#   11.High level term (what we will use for the Cre release)
#
# 2010-04-19	Joel Richardson
#	- modified again:
#	  1. Remove email addresses from col 8 (holder site)
#	  2. Make this a regular report, run weekly:
#		- do not set server/database variables. These are
#		already set in the env
#		- set the report output directory 
#		to os.environ['REPORTOUTPUTDIR']
#		- change script name to ALL_CRE.py
#
"""

import sys
import os
import types
import urllib
import reportlib

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
	db.setAutoTranslate(False)
        db.setAutoTranslateBE()
    else:
        import db
except:
    import db


eUtilUrl = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

CRT = reportlib.CRT
SPACE = reportlib.SPACE
TAB = reportlib.TAB
PIPE = "|"

REF_MGI_TYPE_KEY	= 1
ALLELE_MGI_TYPE_KEY	= 11
PUBMED_LDB_KEY		= 29
MGI_LDB_KEY		= 1
DRIVER_NOTE_TYPE_KEY	= 1034
INDUCIBLE_NOTE_TYPE_KEY	= 1032
SYSTEM_VOCAB_KEY	= 75
ALLELE_ORIGREF_TYPE_KEY	= 1011
IMSR_STRAIN_TYPE	= 1
IMSR_GENE_TYPE		= 2
IMSR_ALLELE_TYPE	= 3
IMSR_FACIL_TYPE		= 4

#
# Main
#

alleles = {}
allelesById = {}
facilities = {}

def parseAllele( r ):
    global alleles
    ak = r['_Allele_key']
    if not alleles.has_key(ak):
	allele = r.copy()
	if not allele['alleleType'].startswith('Transgen') and allele['mname'] != None:
	    allele['name'] = allele['mname'] + '; ' + allele['name']
	allele['structure'] = {}
	allele['system'] = {}
	allele['pmid'] = ''
	allele['inducible'] = 'N'
	allele['strains'] = []
	alleles[ak] = allele
	allelesById[ r['accID'] ] = allele
    else:
	allele = alleles[ak]
    if r['structure'] is not None:
	allele['structure'][r['structure']] = 1
    if r['system'] is not None:
	allele['system'][r['system']] = 1

def parseRef( r ):
    global alleles
    ak = r['_Allele_key']
    a = alleles.get(ak,None)
    if a :
	a['journal'] = r['journal'] and r['journal'] or ''
	a['authors'] = r['authors'] and r['authors'] or ''

def parseInducible( r ):
    global alleles
    ak = r['_Object_key']
    if alleles.has_key(ak):
	alleles[ak]['inducible'] = 'Y'

def parsePMID( r ):
    global alleles
    ak = r['_Allele_key']
    if alleles.has_key(ak):
	alleles[ak]['pmid'] = r['accID']

def parseFacility(r):
    global facilities
    facilities[ r['_facility_key'] ] = r

def parseStrain(r):
    id = r['accID']
    allele = allelesById[id]
    fk = r['_facility_key']
    facility = facilities[fk]
    r['facility'] = facility
    allele['strains'].append(r)

def fetchPubmed(ids):
    if type(ids) is types.StringType:
        ids = ids.replace(',',' ').split()
    params = urllib.urlencode(
       {'db': 'pubmed',
        'tool' : 'MGI',
        'email' : 'mgi-help@informatics.jax.org',
        'retmode' : 'text',
        'rettype' : 'medline',
        'id' : ','.join(ids)
        })
    f = urllib.urlopen(eUtilUrl, params)
    rec = {}
    tag = None
    for line in f:
        if len(line.strip()) == 0:
            if len(rec) > 0:
                yield rec
                rec = {}
                tag = None
                continue
        if line[0] != " " and line[4:6] == "- ":
            tag = line[0:4].strip()
            data = line[6:].strip()
            rec.setdefault(tag, []).append(data)
            continue
        if line.startswith("      "):
            rec[tag][-1] += " " + line.strip()

def main():
    db.useOneConnection(1)
    fp = reportlib.init(sys.argv[0], outputdir = os.environ['REPORTOUTPUTDIR'], printHeading = None)

    # get basic info from Cre cache table and accession table
    cmd = '''
	select a._Allele_key, a.symbol, a.name, mm.name as mname, a.alleleType,
	       a.driverNote, a.structure, a.system, aa.accID
	from ALL_Cre_Cache a, ACC_Accession aa, 
	     ALL_Allele alle LEFT OUTER JOIN MRK_Marker mm on (alle._Marker_key = mm._Marker_key)
	where a._Allele_key = aa._Object_key
	and aa._MGIType_key = %d
	and aa._LogicalDB_key = %d
	and aa.preferred = 1
	and a._Allele_key = alle._Allele_key
	''' % (ALLELE_MGI_TYPE_KEY, MGI_LDB_KEY)
    db.sql( cmd, parseAllele )

    # set 'inducible' field by querying for inducible-type notes
    cmd = '''
	select distinct _Object_key
	from MGI_Note
	where _MGIType_key = %d
	and _NoteType_key = %d
	''' % (ALLELE_MGI_TYPE_KEY,INDUCIBLE_NOTE_TYPE_KEY)
    db.sql( cmd, parseInducible )

    # get original reference PMID
    cmd = '''
	select distinct a._Allele_key, aa.accID
	from ALL_Cre_Cache a, ACC_Accession aa, MGI_Reference_Assoc ra
	where a._Allele_key = ra._Object_key
	and ra._MGIType_key = %d
	and ra._RefAssocType_key = %d
	and ra._Refs_key = aa._Object_key
	and aa._MGIType_key = %d
	and aa._logicalDB_key = %d
	''' % (ALLELE_MGI_TYPE_KEY,ALLELE_ORIGREF_TYPE_KEY,REF_MGI_TYPE_KEY,PUBMED_LDB_KEY)
    db.sql( cmd, parsePMID )

    # get Journal, Author for cases when there's no PMID
    cmd = '''
	select distinct a._Allele_key, br.journal, br.authors
	from ALL_Cre_Cache a, MGI_Reference_Assoc ra, BIB_Refs br
	where a._Allele_key = ra._Object_key
	and ra._MGIType_key = %d
	and ra._RefAssocType_key = %d
	and ra._Refs_key = br._Refs_key
	''' % (ALLELE_MGI_TYPE_KEY,ALLELE_ORIGREF_TYPE_KEY)
    db.sql( cmd, parseRef )

    # get facility data from IMSR
    cmd = '''
	select _facility_key, baseURL, siteName, abbrevName
	from imsr..Facility
	'''
    db.sql(cmd, parseFacility)

    # get IMSR strain data for all MGI Cre alleles
    cmds = [
	# Select Cre alleles into a temp table
	'''
	select distinct a._Allele_key, aa.accID
	into #ids
	from ALL_Cre_Cache a, ACC_Accession aa
	where a._Allele_key = aa._Object_key
	and aa._MGIType_key = %d
	and aa._LogicalDB_key = %d
	and aa.preferred = 1
	''' % (ALLELE_MGI_TYPE_KEY, MGI_LDB_KEY),

	# create an index on it
	'''
	create index ids_ix
	on #ids(accID)
	''',

	# Join on MGI id of the Cre allele to IMSR strain.
	'''
	select distinct i.accID, lbl.label, sfa._facility_key, sfa.strainURL
	from #ids i, imsr..Accession a, imsr..SGAAssoc sga, 
	     imsr..Label lbl, imsr..StrainFacilityAssoc  sfa
	where i.accID = a.accID
	and a._IMSRType_key = %d
	and a._Object_key = sga._Allele_key
	and sga._Strain_key = lbl._Object_key
	and lbl._IMSRType_key = %d
	and lbl.labelType = "N"
	and sga._Strain_key = sfa._Strain_key
	order by i.accID
	''' % (IMSR_ALLELE_TYPE,IMSR_STRAIN_TYPE)
	]
    db.sql(cmds, [None,None,parseStrain])

    aids = allelesById.keys()
    aids.sort()

    #------------------------------------------
    pmid2email = {}
    ''' # disable getting email addresses
    pmids = set()
    for id in aids:
	allele = allelesById[id]
	if len(allele['strains']) == 0:
	    pmids.add(allele['pmid'])

    for pmRec in fetchPubmed(pmids): 
	if pmRec.has_key('AD'):
	    affil = pmRec['AD'][0]
	    email = affil.split()[-1]
	    if email.find('@') > 0:
		pmid2email[ pmRec['PMID'][0] ] = email
    '''
    #------------------------------------------

    for id in aids:
	allele = allelesById[id]
	if len(allele['strains']) == 0:
	    email = pmid2email.get(allele['pmid'],'')
	    allele['strains'] = [{
		'label' : '',
		'facility' : {'abbrevName':email},
		'strainURL' : ''
		}]
	if allele['pmid'] == '':
	    allele['pmid']=allele['journal']+', '+allele['authors']
	for srec in allele['strains']:
	    orec = [
		allele['accID'],
		allele['symbol'],
		allele['name'],
		srec['label'],
		allele['driverNote'],
		allele['pmid'],
		allele['inducible'],
		srec['facility']['abbrevName'],
		srec['strainURL'],
		PIPE.join(allele['structure'].keys()),
		PIPE.join(allele['system'].keys()),
		]
	    fp.write( TAB.join(orec) + CRT)

    reportlib.finish_nonps(fp)
    db.useOneConnection(0)


main()
