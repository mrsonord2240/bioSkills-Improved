# Screen Analysis: Efficiency Filtering, Bystander Attribution, Hit Calling

Read when parsing CRISPResso2 output after a BE screen: filtering sgRNAs by editing efficiency, attributing signal to target vs bystander edits, and aggregating sgRNA scores to variants.

## Editing Efficiency Filtering (Critical Pre-Hit-Calling)

**Goal:** Drop sgRNAs that do not edit efficiently, since unedited reads represent no biological perturbation.

**Approach:** From CRISPResso2 output, compute target-base-conversion percentage per sgRNA; filter library to sgRNAs with >50% target editing in a pilot or co-screened control.

```python
def filter_by_editing_efficiency(crispresso_outputs_dir, target_pos, target_base, efficiency_threshold=0.5):
    '''Drop sgRNAs that edit <efficiency_threshold of reads at target position.
    crispresso_outputs_dir: directory containing CRISPResso per-sample outputs.
    target_pos: 1-indexed position WITHIN THE QUANTIFICATION WINDOW, in column
    order -- CRISPResso2 does not emit a literal "Position" column.'''
    from pathlib import Path
    results = []
    for sample_dir in Path(crispresso_outputs_dir).glob('CRISPResso_on_*'):
        sgrna_id = sample_dir.name.replace('CRISPResso_on_', '')
        quant_file = sample_dir / 'Quantification_window_nucleotide_percentage_table.txt'
        if not quant_file.exists():
            continue
        # Real CRISPResso2 2.3.4 file layout (verified against actual output, not
        # assumed): rows are nucleotide identity (A/C/G/T/N/-, the index column);
        # columns are one per window position, header-labeled with the REFERENCE
        # base at that position (so headers repeat -- pandas suffixes duplicates
        # .1/.2/... -- select columns positionally, not by label). Values are
        # FRACTIONS in [0, 1], not 0-100, despite the filename.
        df = pd.read_csv(quant_file, sep='\t', index_col=0)
        # Schema check: fail loudly and specifically on drift instead of a bare
        # KeyError deep in a groupby/indexing call.
        if target_base not in df.index:
            raise ValueError(f"target_base={target_base!r} not in table rows {list(df.index)} ({quant_file}); "
                              "unexpected CRISPResso2 Quantification_window_nucleotide_percentage_table.txt schema")
        if not (1 <= target_pos <= df.shape[1]):
            raise ValueError(f"target_pos={target_pos} out of range for a {df.shape[1]}-position "
                              f"quantification window in {quant_file}")
        original_frac = df.loc[target_base].iloc[target_pos - 1]
        editing_pct = 1 - original_frac
        results.append({'sgrna_id': sgrna_id, 'editing_pct': editing_pct,
                         'pass_filter': editing_pct >= efficiency_threshold})
    return pd.DataFrame(results)
```

**Convention:** Drop sgRNAs below 50% editing for variant-function screens. A common working split is a 30% editing floor for primary screening and a 50% floor for confirmed hits. Below 30%, the screen has insufficient power; above 70%, results approach saturation editing.

## Bystander Edit Attribution

**Why this matters:** When a sgRNA's editing window contains the target base AND a bystander base, the screen scores the combination. To attribute screen signal to the target variant alone, either (a) include sgRNAs that edit only the target (no bystander) -- often impossible -- or (b) deconvolute via parallel measurements.

**Strategies for variant-by-variant attribution:**

1. **Tile multiple sgRNAs with different bystander patterns:** If 5 different sgRNAs all hit the target base but have different bystanders, common signal across them is target-attributable (Hanna 2021 approach).

2. **Use orthogonal chemistry:** Run the same variant scan with prime editor (no bystanders); cross-validate. See [[prime-editing-screens]].

3. **Bystander stratification:** From CRISPResso2 allele table, partition reads by exact edit pattern (target only, target+bystander_1, target+bystander_2, etc.); separately score each pattern's contribution to the phenotype.

4. **Restrict library:** Use only sgRNAs with zero bystanders in the editing window (rare; may exclude most candidate spacers).

```python
def deconvolute_bystander(allele_table_path, target_pos, bystander_pos_list):
    '''From CRISPResso2 allele table, partition reads by edit pattern at target + bystanders.
    Returns: per-pattern frequency for each combination of target/bystander edits.'''
    alleles = pd.read_csv(allele_table_path, sep='\t', compression='zip')
    required_cols = {'Aligned_Sequence', 'Reference_Sequence', '%Reads'}
    missing = required_cols - set(alleles.columns)
    if missing:
        raise ValueError(f"Unexpected Alleles_frequency_table schema: missing {missing}; "
                          f"got columns {list(alleles.columns)}")
    # Mark target_edited and per-bystander_edited
    alleles['target_edited'] = alleles['Aligned_Sequence'].str[target_pos-1] != alleles['Reference_Sequence'].str[target_pos-1]
    for bp in bystander_pos_list:
        alleles[f'bystander_{bp}_edited'] = alleles['Aligned_Sequence'].str[bp-1] != alleles['Reference_Sequence'].str[bp-1]
    # Real Alleles_frequency_table.zip has no 'Reference_pct' column -- the
    # per-allele read-fraction column is '%Reads' (verified against actual
    # CRISPResso2 2.3.4 output).
    return alleles.groupby(['target_edited'] + [f'bystander_{bp}_edited' for bp in bystander_pos_list])['%Reads'].sum().reset_index()
```

## Hit Calling for Variant-Function Screens

**Goal:** Score per-variant fitness from a base-editor screen.

**Approach:** Filter library to efficiency-passing sgRNAs (>50% editing), then run MAGeCK MLE or drugZ on the sgRNA-level counts; map each significant sgRNA to its predicted variant + bystander pattern; aggregate to per-variant scores.

```python
def aggregate_variant_scores(mageck_sgrna_summary, variant_annotation_df):
    '''Aggregate sgRNA-level scores to per-variant scores.
    variant_annotation_df: per-sgRNA -> predicted variants (target + bystanders),
    keyed on the same sgRNA-identifier column name as mageck_sgrna_summary.
    MAGeCK's real sgrna_summary.txt column is lowercase 'sgrna' (not 'sgRNA') --
    build variant_annotation_df with that same column name.'''
    for _name, _frame in (('mageck_sgrna_summary', mageck_sgrna_summary), ('variant_annotation_df', variant_annotation_df)):
        if 'sgrna' not in _frame.columns:
            raise ValueError(f"{_name} is missing the 'sgrna' merge column (got {list(_frame.columns)})")
    df = mageck_sgrna_summary.merge(variant_annotation_df, on='sgrna')
    # Target-only contribution: sgRNAs with no bystanders
    target_only = df[df['n_bystanders'] == 0]
    target_only_scores = target_only.groupby('target_variant')['LFC'].agg(['mean', 'std', 'count'])
    # Mixed signal: sgRNAs with bystanders
    mixed = df[df['n_bystanders'] > 0]
    return target_only_scores, mixed
```
