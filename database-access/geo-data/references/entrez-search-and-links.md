# Entrez search and GEO links

### Search GEO for studies matching a query

**Goal:** Find GSE accessions matching keywords + organism + study type.

**Approach:** ESearch on `gds` db with field-qualified terms; filter to `gse[Entry Type]`; summarize with ESummary.

**Reference (BioPython 1.83+):**
```python
from Bio import Entrez
import time

Entrez.email = 'researcher@institution.edu'


def search_geo(term, study_type='gse', organism=None, max_results=50):
    full_term = f'{term} AND {study_type}[Entry Type]'
    if organism:
        full_term += f' AND {organism}[Organism]'
    h = Entrez.esearch(db='gds', term=full_term, retmax=max_results)
    s = Entrez.read(h); h.close()
    if not s['IdList']:
        return []
    h = Entrez.esummary(db='gds', id=','.join(s['IdList']))
    summaries = Entrez.read(h); h.close()
    return summaries


for s in search_geo('breast cancer RNA-seq', organism='Homo sapiens', max_results=10):
    # Surface SuperSeries
    relation = s.get('summary', '')
    is_super = 'SuperSeries' in str(relation)
    print(f"  {s['Accession']:12} {s['n_samples']:>4} samples  {'[SuperSeries]' if is_super else '':12}  {s['title'][:60]}")
```

### Link GEO Series to SRA runs (preferred path: pysradb)

```python
from pysradb import SRAweb


def gse_to_srr(gse):
    db = SRAweb()
    srp_df = db.gse_to_srp(gse)
    if srp_df.empty:
        return []
    srp = srp_df['study_accession'].iloc[0]
    srr_df = db.srp_to_srr(srp)
    return srr_df['run_accession'].tolist()


srrs = gse_to_srr('GSE123456')
print(f'GSE123456 -> {len(srrs)} SRR runs')
```

### Find datasets by PubMed citation

```python
def geo_from_pubmed(pmid):
    h = Entrez.elink(dbfrom='pubmed', db='gds', id=pmid)
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDb']:
        return []
    gds_ids = [l['Id'] for l in r[0]['LinkSetDb'][0]['Link']]
    h = Entrez.esummary(db='gds', id=','.join(gds_ids))
    summaries = Entrez.read(h); h.close()
    return summaries
```
