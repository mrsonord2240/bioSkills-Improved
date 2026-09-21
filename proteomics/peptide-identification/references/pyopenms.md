# Database search and FDR with pyOpenMS

Read when searching or filtering inside Python with pyOpenMS 3.5 instead of an external engine.

### Database Search with pyOpenMS

**Goal:** Match tandem mass spectra in an mzML file against a protein FASTA and produce scored PSMs as idXML.

**Approach:** `SimpleSearchEngineAlgorithm` actually scores spectra (the hand-rolled `ProteaseDigestion` loop only digests, it never matches a spectrum). The FASTA must already contain target + decoy sequences concatenated for downstream FDR; decoys carry a recognizable prefix (or set `decoys` to `'true'` to let the engine generate them). Defaults are 10 ppm fragment tolerance and 1 missed cleavage, so set parameters explicitly. The search already annotates `target_decoy` on each hit.

```python
from pyopenms import SimpleSearchEngineAlgorithm, IdXMLFile, PeptideIdentificationList

protein_ids = []
peptide_ids = PeptideIdentificationList()   # pyOpenMS 3.5+: a plain [] raises TypeError
search = SimpleSearchEngineAlgorithm()
p = search.getParameters()
p.setValue('precursor:mass_tolerance', 10.0)
p.setValue('precursor:mass_tolerance_unit', 'ppm')
p.setValue('fragment:mass_tolerance', 0.02)         # HCD Orbitrap; default is 10 ppm
p.setValue('fragment:mass_tolerance_unit', 'Da')
p.setValue('peptide:missed_cleavages', 2)           # default is 1
p.setValue('modifications:fixed', [b'Carbamidomethyl (C)'])
p.setValue('modifications:variable', [b'Oxidation (M)'])
search.setParameters(p)
# spectra are scored against in-silico fragment ions of every candidate peptide
search.search('sample.mzML', 'human_target_decoy.fasta', protein_ids, peptide_ids)

# protein_ids FIRST in load/store -- the OpenMS argument order is fixed
IdXMLFile().store('search_results.idXML', protein_ids, peptide_ids)
```

### Annotate Target/Decoy and Estimate FDR with pyOpenMS

**Goal:** Convert raw PSM scores into q-values and keep only PSMs at 1% FDR.

**Approach:** `PeptideIndexing` maps each PSM back to proteins and flags target vs decoy from the decoy prefix -- needed for idXML from other engines or after changing the FASTA; `SimpleSearchEngineAlgorithm` output above is already annotated and can go straight to `FalseDiscoveryRate`. `FalseDiscoveryRate.apply` runs the concatenated competition; `IDFilter` keeps q <= 0.01. This is the real pyOpenMS path -- not a hand-rolled decoy/target ratio of unknown provenance.

```python
from pyopenms import PeptideIndexing, FalseDiscoveryRate, IDFilter, FASTAFile

fasta = []
FASTAFile().load('human_target_decoy.fasta', fasta)
indexer = PeptideIndexing()
params = indexer.getParameters()
params.setValue('decoy_string', 'DECOY_')      # must match the decoy prefix in the FASTA
params.setValue('decoy_string_position', 'prefix')
indexer.setParameters(params)
indexer.run(fasta, protein_ids, peptide_ids)   # sets target/decoy flags on every hit

FalseDiscoveryRate().apply(peptide_ids)         # concatenated competition -> per-PSM q-value as the new score
IDFilter().filterHitsByScore(peptide_ids, 0.01) # 0.01 = 1% FDR, the community list-level standard
IDFilter().removeDecoyHits(peptide_ids)
```
