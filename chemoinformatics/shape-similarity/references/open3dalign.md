## Open3DAlign (RDKit)

Open3DAlign uses MMFF atom types and partial charges to find an atom-based 3D alignment:

**Goal:** Align a target molecule onto a query in 3D and score volume overlap with Open3DAlign.

**Approach:** Build 3D structures for query and target, run `GetO3A`, and call `Align()` to transform the probe in place. `Score()` is the unnormalized O3A objective, not a shape Tanimoto or ROCS TanimotoCombo. If a normalized shape similarity is required, compute `1 - rdShapeHelpers.ShapeTanimotoDist(...)` after alignment.

```python
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdShapeHelpers

query = Chem.MolFromSmiles('CCC(=O)Nc1ccccc1')
query = Chem.AddHs(query)
AllChem.EmbedMolecule(query, AllChem.ETKDGv3())

target = Chem.MolFromSmiles('CCC(=O)Nc1ccc(F)cc1')
target = Chem.AddHs(target)
AllChem.EmbedMolecule(target, AllChem.ETKDGv3())

O3A = rdMolAlign.GetO3A(target, query)
rmsd = O3A.Align()  # aligns target to query in place
o3a_score = O3A.Score()
shape_tanimoto = 1.0 - rdShapeHelpers.ShapeTanimotoDist(target, query)
```

`GetO3A` finds an alignment between conformers; `Align()` applies it and returns RMSD. Keep `o3a_score` and normalized `shape_tanimoto` distinct in outputs. `Align()` mutates the probe's coordinates: copy the probe (`Chem.Mol(target)`) or store the best-aligned conformer if the coordinates are part of the output.

**Open3DAlign vs ROCS:** Open3DAlign is open-source and competitive on small benchmarks; slower than ROCS at scale.

## Conformer-Ensemble Shape Searching

For each library molecule, generate ensemble of conformers; pick best-shape conformer:

**Goal:** Run shape-similarity search over a conformer ensemble per library molecule so bound-conformer-like shapes are recovered.

**Approach:** For each library molecule, reject disconnected fragments (salts have no single shape to compare), add hydrogens, embed n_conf conformers with ETKDGv3, MMFF-optimize, drop only the conformers that fail to converge (not the whole molecule), score the surviving conformers against the query with Open3DAlign, and keep the best score per molecule. Report every molecule that is dropped and why -- do not let a molecule silently disappear from the results.

```python
def shape_search_ensemble(query_mol, library_mols, n_conf=20):
    hits = []
    dropped = []  # (smiles, reason) for every molecule that produced no usable score
    for target in library_mols:
        smi = Chem.MolToSmiles(target)
        if len(Chem.GetMolFrags(target)) > 1:
            dropped.append((smi, 'disconnected fragments (salt/multi-component); '
                                  'no single shape to compare'))
            continue

        target = Chem.AddHs(target)
        ids = list(AllChem.EmbedMultipleConfs(target, numConfs=n_conf,
                                               params=AllChem.ETKDGv3()))
        if not ids:
            dropped.append((smi, 'embedding failed for all requested conformers'))
            continue
        if not AllChem.MMFFHasAllMoleculeParams(target):
            dropped.append((smi, 'MMFF parameters unavailable'))
            continue

        optimization = AllChem.MMFFOptimizeMoleculeConfs(target)
        # Keep only the conformer ids that converged (status == 0); a molecule with
        # 19 good conformers and 1 non-convergent one must not be discarded outright.
        converged_ids = [cid for cid, (status, _) in zip(ids, optimization) if status == 0]
        n_failed = len(ids) - len(converged_ids)
        if n_failed:
            print(f'WARNING: {smi}: {n_failed}/{len(ids)} conformers failed MMFF '
                  f'convergence; scoring the remaining {len(converged_ids)}')
        if not converged_ids:
            dropped.append((smi, f'all {len(ids)} conformers failed MMFF convergence'))
            continue

        scores = []
        for c in converged_ids:
            O3A = rdMolAlign.GetO3A(target, query_mol, prbCid=c)
            O3A.Align()
            scores.append(1.0 - rdShapeHelpers.ShapeTanimotoDist(
                target, query_mol, confId1=c,
            ))
        hits.append((target, max(scores)))

    if dropped:
        print(f'WARNING: {len(dropped)} library molecule(s) produced no usable '
              f'conformer and were dropped:')
        for smi, reason in dropped:
            print(f'  {smi}: {reason}')
    return sorted(hits, key=lambda x: x[1], reverse=True)
```

**Critical:** Results depend on conformer coverage. Use an ensemble sized and validated for the library and query (20 conformers is a repository starting budget, not a universal minimum) rather than assuming one conformer is representative. A molecule that loses some but not all of its conformers to non-convergence is still scored on the survivors; only a molecule with zero usable conformers is dropped, and every drop is printed with its SMILES and reason -- never assume "fewer hits than input molecules" means the missing ones failed cleanly.
