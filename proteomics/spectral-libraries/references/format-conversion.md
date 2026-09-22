# Library format conversion and OpenSWATH decoys

Moved from SKILL.md. Read when converting a library between DIA-NN, OpenSWATH and Spectronaut, or building an OpenSWATH TSV/TraML and generating decoys.

## Convert Library Formats

**Goal:** Move a library between DIA-NN, OpenSWATH, and Spectronaut conventions without silently corrupting RT, intensity, modification, or decoy content.

**Approach:** Conversion is renaming columns AND reconciling units, not a copy. Check RT units (iRT ~ -25..150 vs normalized 0-1 vs minutes), intensity scaling (relative vs absolute), and modification notation (UniMod:35 vs +15.9949 vs Oxidation). For OpenSWATH, generate decoys with OpenSwathDecoyGenerator -- a target-only library has no null.

```python
import pandas as pd

# Spectronaut -> DIA-NN column mapping; iRT and RelativeIntensity are renamed, not recomputed.
SPECTRONAUT_TO_DIANN = {'ModifiedPeptide': 'ModifiedPeptide', 'iRT': 'iRT',
                        'RelativeIntensity': 'LibraryIntensity', 'FragmentMz': 'ProductMz',
                        'FragmentNumber': 'FragmentSeriesNumber', 'PrecursorMz': 'PrecursorMz',
                        'PrecursorCharge': 'PrecursorCharge', 'FragmentCharge': 'FragmentCharge',
                        'FragmentType': 'FragmentType', 'Genes': 'Genes'}

def spectronaut_to_diann(lib):
    out = lib.rename(columns=SPECTRONAUT_TO_DIANN)
    assert out['iRT'].between(-50, 200).all(), 'RT not in iRT units; check column before converting'
    return out
```

**OpenSwathDecoyGenerator's real input requirements (verified on OpenMS 3.5.0):**
`TargetedFileConverter` converts a TSV without complaint even when it is unusable. The transition
list needs ALL of: a literal `Annotation` column (e.g. `y3^1`), chemically real theoretical fragment
m/z (not placeholders), and a per-precursor grouping column -- `transition_group_id`, unique per
PeptideSequence + PrecursorCharge (`FullUniModPeptideName` also works). Without the grouping
column, distinct peptides sharing a charge collapse into one `<Peptide id="_2">` group with no error
(2 targets become 1), and the decoy step then fails or reports wrong counts. Missing `Annotation` or
placeholder m/z gives `Number of decoy peptides: 0`. Check the peptide count after conversion:

```python
import pandas as pd
from pyteomics import mass

def build_openswath_tsv(peptides, path, n_frag=6):
    """peptides: [(sequence, charge, protein, iRT)] -> OpenSWATH TSV of y-ion transitions."""
    rows = []
    for seq, z, prot, irt in peptides:
        group = f'{seq}_{z}'  # unique per PeptideSequence + PrecursorCharge; required
        for i in range(1, n_frag + 1):
            rows.append({'PrecursorMz': mass.fast_mass(seq, charge=z),
                         'ProductMz': mass.fast_mass(seq[-i:], ion_type='y', charge=1),  # real m/z
                         'Tr_recalibrated': irt, 'transition_name': f'{group}_y{i}',
                         'transition_group_id': group, 'decoy': 0, 'LibraryIntensity': 1000.0 / i,
                         'PeptideSequence': seq, 'FullUniModPeptideName': seq,
                         'ProteinName': prot, 'PrecursorCharge': z, 'FragmentType': 'y',
                         'FragmentSeriesNumber': i, 'FragmentCharge': 1,
                         'Annotation': f'y{i}^1'})  # literal Annotation; required
    pd.DataFrame(rows).to_csv(path, sep='	', index=False)
```

```bash
TargetedFileConverter -in library.tsv -in_type tsv -out library.TraML -out_type TraML
grep -c "<Peptide " library.TraML   # must equal the number of distinct sequence+charge precursors
OpenSwathDecoyGenerator -in library.TraML -out library_decoy.TraML -method pseudo-reverse
```

The default `-method shuffle` has no seed flag and is NOT reproducible: two runs on identical
input produce different decoy peptide sequences. Use `-method pseudo-reverse` when decoys must be
reproducible (5/5 repeated runs byte-identical). `-method reverse` gives reproducible decoy
sequences and fragments, but repeated runs are not always byte-identical (3 distinct hashes in 8
runs): only the last digits of the isolation-window target m/z metadata float vary. `-method shift`
is listed by `--helphelp` but rejected every peptide as a duplicate in testing (OpenMS 3.5.0)
because it leaves the amino-acid sequence unchanged; do not rely on it.
