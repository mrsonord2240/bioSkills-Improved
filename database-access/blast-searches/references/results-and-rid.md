# Results handling and RID polling

Saving XML, extracting hits, and polling a long job by RID. Read when parsing results or when a search will outlive one blocking `qblast()` call. RID lifecycle table: `SKILL.md`.

### Save XML for re-parsing

```python
handle = NCBIWWW.qblast('blastn', 'refseq_select_rna', query)
with open('blast.xml', 'w') as f:
    f.write(handle.read())
handle.close()

with open('blast.xml') as f:
    record = NCBIXML.read(f)
```

### Hit extraction with identity + coverage filtering

**Goal:** Return structured top hits with biological metrics, not just E-values.

**Approach:** Walk alignments + first HSP; compute identity and query coverage as fractions; sort by bit-score (database-size invariant) not E-value.

**Reference (BioPython 1.83+):**
```python
def top_hits(record, min_identity=0.5, min_coverage=0.7, top_n=10):
    qlen = record.query_length
    hits = []
    for aln in record.alignments:
        hsp = aln.hsps[0]
        ident = hsp.identities / hsp.align_length
        cov = hsp.align_length / qlen
        if ident >= min_identity and cov >= min_coverage:
            hits.append({
                'accession': aln.accession,
                'title': aln.title,
                'evalue': hsp.expect,
                'bits': hsp.bits,
                'identity': ident,
                'coverage': cov,
            })
    return sorted(hits, key=lambda h: -h['bits'])[:top_n]
```

### Programmatic RID polling for long jobs

**Goal:** Submit a long job, keep the RID, poll without blocking, resume later.

**Approach:** `qblast()` never exposes the RID, so use the NCBI BLAST URL API directly (`scripts/blast_rid.py`, stdlib only): `Put` returns RID + RTOE, `SearchInfo` returns `Status=WAITING|READY|FAILED|UNKNOWN`, `Get` returns the XML. It waits RTOE, polls at most once per 60 s, refuses queries without a defline, and exits with the NCBI error text on a rejected submit.
```bash
python scripts/blast_rid.py run --query q.fa --program blastn --db refseq_select_rna --hitlist 500 --expect 1e-10 --out hits.xml
python scripts/blast_rid.py submit --query q.fa --program tblastn --db nr      # prints RID; later: status RID / fetch RID --out hits.xml
```
The saved XML parses with `NCBIXML.read()` (see `examples/save_and_parse.py`).
