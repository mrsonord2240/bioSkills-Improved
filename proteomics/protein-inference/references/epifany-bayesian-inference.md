# Bayesian Inference with EPIFANY

### Bayesian Inference + Group FDR with EPIFANY

**Goal:** Assign calibrated protein/group posteriors and control protein-group FDR with a probability model rather than greedy parsimony.

**Approach:** EPIFANY consumes idXML whose PSMs already carry posterior error probabilities (from Percolator or IDPosteriorErrorProbability), then propagates belief over the peptide-protein graph. The TOPP tool is `Epifany`; the pyOpenMS class is `BayesianProteinInferenceAlgorithm`.

```python
from pyopenms import IdXMLFile, BayesianProteinInferenceAlgorithm, PeptideIdentificationList

protein_ids = []
peptide_ids = PeptideIdentificationList()
IdXMLFile().load('peptides_with_pep.idXML', protein_ids, peptide_ids)

algo = BayesianProteinInferenceAlgorithm()
# EPIFANY expects PSM posteriors as input. The third POSITIONAL argument is
# greedy_group_resolution. Unlike BasicProteinInferenceAlgorithm's parameter of the same
# name, flipping it here did not change the result on the reference idXML (583 groups at 1%
# picked FDR, 3.43% true FDP with either value) -- so do NOT assume it removed subsumable
# proteins. Check the surviving groups, or use Basic + greedy when FDP is what you care about.
algo.inferPosteriorProbabilities(protein_ids, peptide_ids, False)

for group in protein_ids[0].getIndistinguishableProteins():
    print([a.decode() for a in group.accessions], group.probability)   # higher = more likely present
```
