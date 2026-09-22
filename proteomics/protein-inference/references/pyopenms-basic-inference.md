# Group Proteins with pyOpenMS: Basic Inference + Greedy Resolution + Picked-Group FDR

### Group Proteins with pyOpenMS (aggregation + greedy resolution)

**Goal:** Turn an FDR-filtered peptide identification list into protein groups with a leading protein, resolving shared-peptide ambiguity.

**Approach:** Load the idXML from peptide identification, run `BasicProteinInferenceAlgorithm` (score aggregation per protein) with indistinguishable-group annotation AND greedy group resolution on -- without resolution, subsumable and shared-only proteins stay as their own groups -- apply picked protein-group FDR with the built-in, and read the groups off the protein identification run as records: `leading_protein`, `accessions`, `n_peptides`, `n_unique_peptides` (unique to the GROUP), `is_decoy`, `qvalue` -- the contract to bind downstream code to, and the same record `examples/protein_groups.py` returns.

```python
from pyopenms import (IdXMLFile, BasicProteinInferenceAlgorithm, PeptideIdentificationList,
                      FalseDiscoveryRate, String)

DECOY_PREFIX = 'DECOY_'   # tool-specific: DECOY_ OpenMS/Comet, rev_ Philosopher, REV__ MaxQuant

protein_ids = []
peptide_ids = PeptideIdentificationList()   # pyOpenMS 3.5+: a plain [] fails
# protein_ids is FIRST in both load() and store() for IdXMLFile
IdXMLFile().load('peptides_1pct_fdr.idXML', protein_ids, peptide_ids)

inference = BasicProteinInferenceAlgorithm()
params = inference.getParameters()
# annotate_indistinguishable_groups reports indistinguishable proteins as ONE group
params.setValue('annotate_indistinguishable_groups', 'true')
# greedy_group_resolution assigns shared peptides to the best group, so subsumable proteins drop out
params.setValue('greedy_group_resolution', 'true')
inference.setParameters(params)
inference.run(peptide_ids, protein_ids)

# picked protein-group FDR: (decoy string, is prefix, groups too); group.probability becomes a q-value
FalseDiscoveryRate().applyPickedProteinFDR(protein_ids[0], String(DECOY_PREFIX), True, True)

# peptide sequences behind each protein, from the best hit of every PSM
peps_of = {}
for pid in peptide_ids:
    hit = pid.getHits()[0]
    for ev in hit.getPeptideEvidences():
        peps_of.setdefault(ev.getProteinAccession(), set()).add(hit.getSequence().toString())

# group record: leading_protein, accessions, n_peptides, n_unique_peptides, is_decoy, qvalue
# (the same record examples/protein_groups.py returns)
groups = []
for g in protein_ids[0].getIndistinguishableProteins():
    accs = [a.decode() for a in g.accessions]   # bytes; pyOpenMS sorts them alphabetically
    peptides = set().union(*(peps_of.get(a, set()) for a in accs))
    others = set().union(*(s for a, s in peps_of.items() if a not in accs))   # unique = not in any other protein
    groups.append({
        # members share the same evidence; choose the lead explicitly: canonical (no -N isoform suffix) first
        'leading_protein': sorted(accs, key=lambda a: ('-' in a, a))[0],
        'accessions': accs,
        'n_peptides': len(peptides),
        'n_unique_peptides': len(peptides - others),
        'is_decoy': all(a.startswith(DECOY_PREFIX) for a in accs),
        'qvalue': g.probability,
    })
passing = [g for g in groups if not g['is_decoy'] and g['qvalue'] <= 0.01]
print(len(passing), 'groups at 1% picked-group FDR')
for g in passing[:5]:
    print(g)
```
