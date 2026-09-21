# ClinVar and dbSNP retrieval

Read when fetching `clinvar` (VCV) or `snp` records; both return XML that `Entrez.read()` cannot parse.

Assumes the Required Setup block in `SKILL.md` (`Entrez.email`, imports, `expect_start`).

### clinvar

| rettype | retmode | Returns | Use when |
|---|---|---|---|
| `vcv` | `xml` | ClinVar Variation Archive XML | Clinical significance, variant name, accession (checked live 2026-09-17) |

The VCV XML ships with **no DTD or XSD**, so `Entrez.read()` cannot parse it (Biopython raises `ValueError` and recommends `xml.etree.ElementTree` instead — see `clinvar_record()` under Code patterns). ESummary also works (`Entrez.esummary(db='clinvar', id=uid)`) but its docsum references a `common_name` tag missing from Biopython's cached DTD; pass `validate=False` to `Entrez.read()` or use EFetch instead.

ClinVar records report clinically sensitive information (pathogenicity classifications). Report only the record's own stated classification — never infer a diagnosis or a treatment recommendation from it; defer those to a clinician or genetic counselor.

### snp

| rettype | retmode | Returns | Use when |
|---|---|---|---|
| `xml` | `xml` | dbSNP docsum XML (alleles, gene, clinical significance, frequencies) | Variant-level metadata for an rs# / numeric SNP UID (checked live 2026-09-17) |

Response is a namespaced (`https://www.ncbi.nlm.nih.gov/SNP/docsum`) `ExchangeSet`/`DocumentSummary` document; parse with `xml.etree.ElementTree` and the namespace map, same as clinvar — see `snp_record()` under Code patterns.

### ClinVar record by UID

**Goal:** Get a variant's clinical significance from a ClinVar UID.

**Approach:** `rettype='vcv', retmode='xml'`. The response has no DTD/XSD, so `Entrez.read()` cannot parse it (Biopython raises `ValueError` and recommends ElementTree) — parse with `xml.etree.ElementTree` instead.

```python
import xml.etree.ElementTree as ET

def clinvar_record(uid):
    h = Entrez.efetch(db='clinvar', id=uid, rettype='vcv', retmode='xml')
    raw = h.read(); h.close()
    text = raw.decode() if isinstance(raw, bytes) else raw
    archive = ET.fromstring(text).find('.//VariationArchive')
    sig = archive.find('.//Classifications/GermlineClassification/Description')
    return {
        'accession': archive.get('Accession'),
        'variation_name': archive.get('VariationName'),
        'clinical_significance': sig.text if sig is not None else None,
    }
```

Report only the record's own stated classification — never a diagnosis or treatment recommendation.

### dbSNP record by UID

**Goal:** Get alleles, gene, and clinical significance for a numeric SNP UID (e.g. `rs429358` -> UID `429358`).

**Approach:** `rettype='xml', retmode='xml'`; response is namespaced, parse with `xml.etree.ElementTree`.

```python
def snp_record(uid):
    h = Entrez.efetch(db='snp', id=uid, rettype='xml', retmode='xml')
    raw = h.read(); h.close()
    text = raw.decode() if isinstance(raw, bytes) else raw
    ns = {'s': 'https://www.ncbi.nlm.nih.gov/SNP/docsum'}
    doc = ET.fromstring(text).find('s:DocumentSummary', ns)
    gene = doc.find('.//s:GENE_E/s:NAME', ns)
    sig = doc.find('s:CLINICAL_SIGNIFICANCE', ns)
    return {
        'chr': doc.findtext('s:CHR', default=None, namespaces=ns),
        'gene': gene.text if gene is not None else None,
        'clinical_significance': sig.text if sig is not None else None,
    }
```
