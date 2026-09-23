---
name: bio-proteomics-ptm-analysis
description: Frames PTM/phosphoproteomics analysis as three stacked inference layers on a biased enrichment - chemistry selection, site localization (FLR), and protein-level-adjusted quantification with MSstatsPTM - plus kinase-activity and functional triage. Covers MaxQuant Phospho (STY)Sites multiplicity expansion, localization-probability filtering (class I, Ascore, ptmRS, DIA EG.PTMLocalizationProbabilities, DIA-NN PTM.Site.Confidence), false localization rate (LuciPHOr/DeepFLR), motif analysis with experiment-matched backgrounds, diGly/K-GG ubiquitin specificity, acetyl/glyco traps, and KSEA/PTM-SEA. Use when localizing and quantifying phosphorylation, acetylation, ubiquitination, or glycosylation sites from enrichment-based runs and deciding whether an apparent site change is real after subtracting protein abundance. Peptide ID and open/variable-mod search is peptide-identification; underlying protein-level quant is quantification and differential-abundance; DIA acquisition mechanics is dia-analysis.
tool_type: mixed
primary_tool: MSstatsPTM
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: MSstatsPTM 2.8.1, MSstatsTMT 2.14.2 (the TMT route), KSEAapp 2.0, pandas 2.2+, numpy 1.26+, scipy 1.12+ (checked 2026-09-15 on R 4.4.3)

Install: `pip install pandas numpy scipy`; `BiocManager::install('MSstatsPTM')` (plus `MSstatsTMT` for the TMT route); `install.packages('KSEAapp')`; PTM-SEA is `broadinstitute/ssGSEA2.0` on GitHub (R scripts, not a package).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# PTM and Phosphoproteomics Analysis -- Three Inference Layers Stacked on a Biased Extraction

**"Find the regulated phosphosites in my enriched samples"** -> Localize each modification, then test whether its abundance change survives subtracting the protein-level change -- because a PTM result is three separate inferences (enrichment, localization, quantification) and each fails silently if the layer below is treated as solved.
- R: `MSstatsPTM::groupComparisonPTM()` for protein-adjusted site testing (the load-bearing tool)
- Python: `pandas` to expand MaxQuant `Phospho (STY)Sites` multiplicity and filter localization probability
- R: `KSEAapp` / PTM-SEA (`ssGSEA2.0`) for kinase-activity inference from the site fold-changes

Scope: this skill OWNS enrichment-chemistry framing, site localization and FLR, multiplicity-resolved site quant, protein-level adjustment (label-free AND TMT/isobaric), motif analysis, and kinase-activity inference. Peptide identification and open/variable-mod search route to peptide-identification; the underlying protein-level (unenriched) quant routes to quantification and differential-abundance; DIA acquisition mechanics route to dia-analysis. OUT OF SCOPE: intact-glycopeptide glycan-composition search (pGlyco3/MSFragger-Glyco) and absolute occupancy from three-ratio SILAC are noted but not implemented here.

## The Single Most Important Modern Insight -- A PTM Result Is Three Inferences, Not One

1. **Enrichment IS the experiment, and the chemistry is a filter confounded with biology.** The data only contain what the chemistry captured. TiO2 and Fe-IMAC give partially-overlapping phosphoproteomes; anti-K-GG enriches ubiquitin, NEDD8, and ISG15 indistinguishably; a lectin reports only its cognate glycoforms. A between-method or between-lab "biological difference" must FIRST be excluded as a chemistry artifact before it is called biology.
2. **Identifying a peptide is NOT localizing the modification.** A phosphopeptide with two S/T and one phosphate has isobaric positional isomers of identical precursor mass and identical peptide-level score; the localization is a SECOND inference decided only by site-determining fragment ions, with its own error rate (false localization rate, FLR). Target-decoy peptide FDR cannot estimate FLR: a wrong localization is the correct sequence with the mod one residue over, not a decoy sequence (Fermin 2013). A 1% peptide FDR does NOT yield a 1% site FDR -- report them separately.
3. **A change in phosphopeptide abundance is NOT a change in phosphorylation (the biggest quant trap).** Observed PTM signal ~ (site occupancy) x (protein abundance) x (enrichment/ionization factor), so `log2FC(PTM_observed) = log2FC(occupancy) + log2FC(protein)`. Without a paired global (unenriched) proteome run on the SAME samples to subtract `log2FC(protein)`, every protein-abundance change masquerades as a regulated site. Because co-regulated proteins move together, the false positives are pathway-coherent and look biologically convincing -- the worst kind of artifact. This is the entire reason MSstatsPTM exists (Kohler 2023).
4. **Most identified sites have no known function.** Fewer than ~5% of phosphosites are functionally annotated; a fold-change alone says nothing about regulatory relevance. Functional triage (conservation, stoichiometry, Ochoa functional score, confident kinase assignment) is a separate fourth layer on top of the quant (Ochoa 2020).

Bottom line: report THREE numbers, not one -- peptide/PSM FDR, per-site localization probability with its threshold, and an empirically estimated global FLR -- and never call a site "regulated" from a phospho-only run without protein-level adjustment.

## Tool Taxonomy

### Enrichment chemistry (phospho)

| Method | Citation | Mechanism / bias | When |
|---|---|---|---|
| TiO2 (MOAC) | Larsen 2005 | Metal-oxide Lewis-acid surface; skews mono-phospho; needs hydroxy-acid additive | General single-method depth; EasyPhos basis |
| Fe(III)-IMAC | Ruprecht 2015 | Chelated Fe3+ coordinates phosphate; multidentate avidity skews MULTI-phospho | Hierarchical/processive signaling; the mono/multi divergence is STRONGEST here vs TiO2 |
| Ti4+/Zr4+-IMAC | Matheron 2014 | Chelated metal ION on immobilized phosphonate (NOT bulk oxide); bias vs TiO2 is SMALL | Modern automated workflows; metal identity matters more than IMAC-vs-MOAC |
| SIMAC (sequential) | Thingholm 2008 | IMAC acidic elution = mono, basic = multi, then TiO2 on mono fraction | Recovering both populations IMAC alone biases |

Naming trap: Ti4+/Zr4+-IMAC (chelated ions) is DIFFERENT chemistry from TiO2/ZrO2 (bulk oxide). Glycolic acid is the modern additive standard (load 80% ACN / 5% TFA / 0.1 M glycolic acid).

### Other-PTM enrichment and identity traps

| PTM | Reagent / mass | Citation | Headline trap |
|---|---|---|---|
| Ubiquitin (diGly, K-GG, +114.0429) | Anti-K-GG antibody | Xu 2010; Kim 2011 | NOT ubiquitin-specific: K-GG = ubiquitin + NEDD8 (~6% at basal) + ISG15 (rises under interferon). UbiSite (Akimov 2018) is the ubiquitin-specific alternative |
| Ubiquitin alkylation artifact | use chloroacetamide | Nielsen 2008 | Iodoacetamide creates a +114.0429 lysine adduct mimicking ubiquitination; chloroacetamide does not |
| Acetyl-K (+42.0106) | Anti-acetyllysine cocktail | Svinkina 2015 | Isobaric with trimethyl +42.0470 (0.0364 Da, needs high-res); acetyl blocks trypsin -> allow >=4 missed cleavages |
| Glyco N-linked | PNGase F (released) or intact | Riley 2021 | Released loses the glycan; N->D tag +0.984 is isobaric with deamidation -- use PNGase F in H2-18O (+2.988) to disambiguate; N-X-S/T (X!=Pro) sequon is necessary not sufficient |

### Localization scoring

| Tool | Citation | Mechanism | Note |
|---|---|---|---|
| Ascore | Beausoleil 2006 | Cumulative binomial of site-determining ions; DIFFERENCE between best and 2nd-best localization, peak-depth sweep | Ascore >=19 ~ p 0.01 PAIRWISE per-PSM, NOT a dataset FLR |
| PhosphoRS / ptmRS | Taus 2011 | Per-isomer cumulative binomial, tolerance-aware (correct for high-res), per-site probs sum to 100% | In Proteome Discoverer |
| PTMProphet | Shteynberg (TPP) | EM/Bayesian mixture; per-site probs combinable across PSMs to a global FLR | TPP/FragPipe |
| MaxQuant Localization prob | Cox/Mann (Andromeda) | Normalized posterior on the site (fixed peak depth) | column `Localization prob`; >=0.75 = class I |
| DIA localization | Bekker-Jensen 2020 | XIC peak-shape correlation substitutes for missing precursor isolation | Spectronaut `EG.PTMLocalizationProbabilities`; DIA-NN `PTM.Site.Confidence` |
| DeepFLR | Zong 2023 | Deep-learning spectrum predictor + target-decoy FLR | SOTA direction; DDA + DIA |

### Site-FDR / FLR (the layer most pipelines skip)

| Tool | Citation | Mechanism |
|---|---|---|
| LuciPHOr | Fermin 2013 (MCP) | Decoy localizations on non-modifiable residues; rate decoys win = empirical FLR |
| Decoy amino-acid FLR | Ramsbottom/Jones 2022 | Add a non-modifiable residue to the candidate set; global FLR ~ decoy-site-hits / target-site-hits, frequency-corrected |

### Kinase-activity inference

| Tool | Citation | Mechanism | Limitation |
|---|---|---|---|
| KSEA | Casado 2013 | z-score of a kinase's substrate fold-changes | sqrt(m) favors well-annotated kinases; inherits PhosphoSitePlus bias |
| PTM-SEA / PTMsigDB | Krug 2019 | ssGSEA2.0 on site-level +/-7 flanking-sequence signatures | robust to isoform drift; PERT signatures score "looks like EGF stim" |
| RoKAI | Yilmaz 2021 | Network-smooth profiles before z-score so unobserved sites borrow neighbor signal | attacks missingness; feeds KSEA |

Benchmark result (Mueller-Dott 2025): across ~19 methods, simple z-score (KSEA/RoKAI) matched or beat sophisticated methods. Performance is PRIOR-limited, not algorithm-limited; all methods inherit PhosphoSitePlus curation bias toward CK2/CDK1/PKA/MAPK, and the dark kinome is structurally invisible. Spend effort on the substrate prior, not the estimator.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|---|---|---|
| Phospho-only run, want regulated sites | Acquire a PAIRED global proteome -> MSstatsPTM `groupComparisonPTM` -> require significance in `ADJUSTED.Model` | Unadjusted site changes are confounded with protein abundance |
| No global proteome available | Report site changes as UNADJUSTED (`PTM.Model`; `ADJUSTED.Model` is absent) and flag the confound explicitly; optionally `use_unmod=TRUE` on `scripts/msstatsptm_labelfree.R` (passed as `use_unmod_peptides`) uses unmodified peptides co-enriched in the PTM runs as a weak protein proxy, which needs those unmodified rows kept in the evidence -- label those results proxy-adjusted | Cannot separate occupancy from abundance; do not claim "regulation" |
| TMT / isobaric plexes (enriched + global) | Same MSstatsPTM route with `labeling_type = 'TMT'`, `dataSummarizationPTM_TMT`, `groupComparisonPTM(data.type = 'TMT')`, and a `Condition = 'Norm'` reference channel in EVERY plex (`references/tmt-isobaric.md`) | The adjustment is identical; only the converter, the summarization function and the test branch differ |
| More than one TMT plex | Refuse to compare until each plex carries the pooled reference channel; `reference_norm = TRUE` bridges them | Reporter ratios are meaningful only within a plex |
| Between-method phospho difference | Suspect chemistry (TiO2 vs Fe-IMAC mono/multi bias) BEFORE biology | Enrichment is a confounded filter |
| Multiply-phospho peptides present | Localize per-site (Ascore/ptmRS) AND report empirical global FLR | Peptide FDR != site FDR |
| "Ubiquitination" sites | Confirm chloroacetamide alkylation; treat K-GG as ub + NEDD8 + ISG15; consider UbiSite | Iodoacetamide artifact + NEDD8/ISG15 confound |
| Which kinases moved? | KSEA or PTM-SEA (`references/kinase-activity-ksea.md`, `references/ptm-sea.md`) with a curated prior; do not over-interpret dark-kinome silence | Prior-limited; simple z-score suffices |
| Motif logo from the hits | Background = experiment-matched S/T/Y from the identified proteins (NOT whole proteome); code in `references/motif-analysis.md` | Whole-proteome background rediscovers disordered-region composition bias |

Default when uncertain: localize with the search engine's probability (class I >=0.75), expand MaxQuant multiplicity, run MSstatsPTM with a paired global proteome, and call only `ADJUSTED.Model` hits regulated.

## Reference Files

| File | Read it when |
|---|---|
| `references/tmt-isobaric.md` | The enriched and global runs are TMT/isobaric plexes: the three calls that differ from label-free, the annotation columns, `Channel` 0-indexing, the `Norm` reference channel, ratio-compression traps |
| `references/motif-analysis.md` | Building a kinase/writer motif from the site windows with an experiment-matched background (Fisher + BH code) |
| `references/kinase-activity-ksea.md` | Scoring kinase activity from the adjusted site fold-changes with KSEAapp: PX layout, FC direction, infinite-fold-change and prior-coverage guards |
| `references/ptm-sea.md` | Scoring PTMsigDB kinase/perturbation/pathway signatures with ssGSEA2.0: GCT construction, CLI call, reading results, what was checked |

Runnable code lives in `scripts/`, each with a header (purpose, inputs, usage):

| Script | Does |
|---|---|
| `scripts/msstatsptm_labelfree.R` | Class-I filter, MSstatsPTM label-free adjustment, TREAT-style calls; writes `adjusted_sites.csv` |
| `scripts/msstatsptm_tmt.R` | The same for TMT/isobaric plexes; writes `adjusted_sites_tmt.csv` |
| `scripts/ksea_scores.R` | KSEA z-scores per kinase from `adjusted_sites.csv` and a kinase-substrate prior, with the failure guards |
| `scripts/motif_enrichment.py` | Fisher + BH motif enrichment against an experiment-matched background |
| `scripts/ptmsea.py` | Build the PTM-SEA GCT from `Sequence window` and read ssGSEA2.0 scores back |

## Expand the MaxQuant Site Table Before Any Quant

**Goal:** Produce a long, multiplicity-resolved, class-I-filtered phosphosite intensity matrix from `Phospho (STY)Sites.txt`.

**Approach:** Each site row spreads its quant across `Intensity___1/___2/___3` (singly/doubly/triply-phospho forms, THREE underscores); the collapsed base `Intensity` mixes phospho-states and can mask a switch between forms. Per-run columns are `Intensity <experiment>___<n>`; the aggregated `Intensity___<n>` (no experiment name) is the sum over runs and must NOT be melted as a run. Drop Reverse/contaminant, filter `Localization prob`, then melt the per-run, per-multiplicity columns into rows.

The runnable recipe is `examples/phospho_analysis.py`: it reads the sites table (spaced or no-space filename), drops Reverse/contaminant rows, keeps `Localization prob` >= 0.75 (`CLASS_I_PROB`), builds `site_id` from gene (or protein) + residue + position, and melts the per-run columns `Intensity <run>___[123]` into a long table with `multiplicity`, `run` and `log2_intensity`. It builds a synthetic table in a temp dir; to run it on your data point `path` at your own `Phospho (STY)Sites.txt`.

## Protein-Level Adjustment with MSstatsPTM

**Goal:** Decide whether each site change is real after subtracting the matched protein-abundance change.

**Approach:** MSstatsPTM carries TWO datasets -- a PTM dataset (enriched) and a PROTEIN dataset (global/unenriched). `groupComparisonPTM` fits independent linear models to each and returns a list of THREE: `PTM.Model` (unadjusted), `PROTEIN.Model`, and `ADJUSTED.Model`. The adjustment is `dFC_adj = dFC_PTM - dFC_protein` with `SE_adj = sqrt(SE_PTM^2 + SE_protein^2)`, so adjustment ADDS uncertainty -- a site can be significant unadjusted yet lose significance after adjustment. A confident regulation call requires significance in `ADJUSTED.Model`.

```bash
Rscript scripts/msstatsptm_labelfree.R dir=<data_dir> out=<out_dir>   # writes <out_dir>/adjusted_sites.csv
```

`scripts/msstatsptm_labelfree.R` runs the whole route, with the traps in its comments: class-I pre-filter on the enriched evidence (best `Phospho (STY) Probabilities` >= 0.75; unmodified rows kept only with `use_unmod=TRUE`); `MaxQtoMSstatsPTMFormat` with the global run's `evidence_prot`, `proteinGroups` and annotation, then `stopifnot('PROTEIN' %in% names(input))`; `dataSummarizationPTM(use_log_file = FALSE, append = FALSE)`; an explicit `Treatment vs Control` contrast; `groupComparisonPTM(data.type = 'LabelFree')`; site rows of `ADJUSTED.Model` only (`Protein_<residue><position>`); and the TREAT-style test `pt((|log2FC| - lfc) / SE, DF)` with BH, in place of a post-hoc `adj.pvalue < 0.05 & |log2FC| > 1` double filter (whose FDR refers to log2FC != 0). Compare `PTM.Model` against `ADJUSTED.Model` to see how much of each call was protein-driven. For TMT plexes see `references/tmt-isobaric.md`.

## A Note on Home-Grown Ascore

The function below is an ILLUSTRATIVE approximation, NOT real Ascore. Real Ascore (Beausoleil 2006) competes the best localization against the second-best, sweeps peak depth 1-10 per 100 Th, and restricts to site-determining ions -- none of which this captures. Use the search engine's own localization probability (MaxQuant `Localization prob`, ptmRS, PTMProphet) for real work, or pyOpenMS `AScore` (introspect the exact API before relying on it). The home-grown form is here only to show the binomial intuition.

```python
import numpy as np
from scipy.stats import binom

def illustrative_localization_score(matched_site_ions, total_ions, depth_p=0.04):
    '''Binomial intuition only; NOT Ascore (no best-vs-second competition or depth sweep).'''
    if total_ions == 0 or matched_site_ions == 0:
        return 0.0
    p_random = 1 - binom.cdf(matched_site_ions - 1, total_ions, depth_p)
    return -10 * np.log10(p_random) if p_random > 0 else 100.0
```

## Per-Method Failure Modes

### Skipping protein-level adjustment
**Trigger:** Differential testing on a phospho-only run with no paired global proteome. **Mechanism:** `log2FC(PTM_observed) = log2FC(occupancy) + log2FC(protein)`; the two terms are inseparable. **Symptom:** Pathway-coherent "regulated sites" that are pure protein-abundance changes (cyclins/histones in cell cycle, stabilized substrates under drug). **Fix:** Run a matched global proteome and adjust via MSstatsPTM; route the protein-level quant to quantification.

### Collapsing the MaxQuant multiplicity
**Trigger:** Quantifying on base `Intensity` instead of `Intensity___1/___2/___3`. **Mechanism:** The collapsed column mixes singly/doubly/triply-phospho forms of the same site. **Symptom:** A switch between forms is masked: the singly-phospho form dropping while the doubly-phospho form rises reads as no change (synthetic test: collapsed +0.16 log2, ___1 -1.40, ___2 +1.56). **Fix:** Expand multiplicity to long form (Perseus "Expand site table" or the melt in `examples/phospho_analysis.py`) before any stats.

### Treating identification as localization
**Trigger:** Reporting sites at peptide FDR without a localization threshold. **Mechanism:** Isobaric positional isomers share precursor mass and peptide score; CID/ion-trap neutral loss (-98 Da) starves site-determining ions. **Symptom:** A 1% peptide FDR result with a much higher true site error. **Fix:** Filter localization probability (class I >=0.75) in EVERY route, including the evidence fed to MSstatsPTM; report an empirical global FLR (LuciPHOr2/DeepFLR), prefer HCD/EThcD. Empirical FLR needs the spectra (mzML); from MaxQuant txt tables alone, report only the model-based expected FLR `mean(1 - Localization prob)` over the accepted sites and label it as such -- it misses mis-assignment sources the search model does not know about.

### diGly read as ubiquitin
**Trigger:** Calling the K-GG proteome "ubiquitination". **Mechanism:** NEDD8 and ISG15 share the LRLRGG C-terminus and leave the identical +114.0429 remnant; iodoacetamide adds a fourth source. **Symptom:** Inflated/false ubiquitin sites, worst under interferon (ISG15) or with iodoacetamide. **Fix:** Ask which alkylation reagent was used; chloroacetamide for future experiments; treat K-GG as ub+NEDD8+ISG15; use UbiSite for ubiquitin-specific mapping. For data already acquired, flag K-GG on the peptide C-terminal lysine (a GG-modified K blocks trypsin, so a C-terminal GG-K is suspect): `ev['Modified sequence'].str.contains(r'K\(GlyGly \(K\)\)_$')`.

### Motif logo against the wrong background
**Trigger:** Whole-proteome or IUPAC-random background. **Mechanism:** Phosphosites sit in disordered, Ser/Pro/acidic-rich regions; that composition dominates the enrichment. **Symptom:** "Enriched" proline/serine motifs that are region bias, not kinase preference. **Fix:** Experiment-matched S/T/Y background or central-residue-preserving shuffle (MoMo default).

### Over-reading kinase-activity output
**Trigger:** Naming the top KSEA/atlas kinase as the responsible enzyme. **Mechanism:** Substrate priors are PhosphoSitePlus-curated (CK2/CDK1/PKA/MAPK heavy); atlas hits are biochemical preference ignoring expression/localization/timing. **Symptom:** Always-the-usual-suspects kinase lists; dark-kinome activity invisible. **Fix:** Use a curated prior, report z-scores with their substrate counts, do not infer absence from silence.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|---|---|---|
| Localization prob class I >= 0.75 | Olsen 2006 | Best site holds >3x the posterior of all alternatives (single-phospho); a comparability standard, not a calibrated error rate |
| Class II 0.5-0.75; class III 0.25-0.5 | Olsen 2006 | Partial / poor localization |
| Ascore >= 19 (p~0.01); loose >13 | Beausoleil 2006 | Pairwise per-PSM best-vs-next confidence; NOT a dataset FLR |
| DIA directDIA localization >= 0.99 | Bekker-Jensen 2020 | Stricter than library-based (0.75) to match DDA error rates |
| DIA-NN site matrix 0.90 / 0.99 | DIA-NN docs | phosphosites_90/99.tsv; class-I-equivalent stringency is HIGHER than MaxQuant 0.75 |
| Acetyl missed cleavages >= 4 | -- | Acetyl-K blocks trypsin; pair LysC + trypsin |
| Kinase atlas motif match >= 90th percentile | Johnson 2023 | Strong motif preference, NOT proof the kinase acted |
| Report peptide FDR and site FLR separately | Fermin 2013 | 1% peptide FDR != 1% site FDR; true site error is typically several-fold higher |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| FileNotFoundError on the sites table | Filename has a SPACE: `Phospho (STY)Sites.txt` | Accept either spaced or no-space form |
| Form switch hidden (no change on a site that switched singly/doubly) | Quantified base `Intensity`, mixing multiplicities | Use `Intensity <run>___1/___2/___3` (three underscores) |
| Extra "run" named `Intensity` in the long table, ~2.7 log2 above the others | Melt also matched the aggregated `Intensity___1/2/3` columns | Select `Intensity <run>___n` (note the space) |
| `Assertion on '!(append & !use_log_file)' failed` | `dataSummarizationPTM` defaults `append = TRUE` | `dataSummarizationPTM(input, use_log_file = FALSE, append = FALSE)` |
| `object 'ptm_model' not found` in `groupComparisonPTM` | `data.type = 'LF'` matches neither branch | `data.type = 'LabelFree'` (or `'TMT'`) |
| `Extra columns included in the annotation file that are not required ... Run, Raw.file, Fraction, TechRepMixture, Channel, Condition, Mixture, BioReplicate` | a label-free annotation passed with `labeling_type = 'TMT'` | Use the TMT annotation columns the message itself lists |
| One of several loud errors from `MaxQtoMSstatsPTMFormat`, depending on which half is mismatched: `A non-empty vector of column names for 'by' is required`; `Extra columns included in the annotation file ... Run, Raw.file, Condition, BioReplicate, IsotopeLabelType`; `Each MS run (Raw.file) can't have multiple conditions or BioReplicates` | TMT evidence left on the default `labeling_type = 'LF'` (or the reverse); none of the messages names the labeling type | Check `labeling_type` first whenever the converter rejects a TMT input: `'TMT'`, then `dataSummarizationPTM_TMT` and `data.type = 'TMT'` |
| `the channel name must be matched with that in input data` | Annotation `Channel` counted from 1 (`channel.1` .. `channel.10`) while MaxQuant's reporter columns start at 0 | Name the channels from the evidence header: `channel.0` .. `channel.9` for a 10-plex |
| `names(input)` is only `PTM`; `ADJUSTED.Model` missing | `evidence_prot` not passed to `MaxQtoMSstatsPTMFormat` | Pass the global run's evidence as `evidence_prot`; `stopifnot('PROTEIN' %in% names(input))` |
| `z.score` NaN for every kinase from `KSEA.Scores` | a site with log2FC -Inf/Inf (missing in one condition) entered PX; FC = 0 passes `is.finite(FC)` | drop rows with non-finite `log2FC` before building PX |
| `arguments imply differing number of rows: 0, 1` while building PX | the non-finite-`log2FC` filter emptied the table, and `Peptide = 'NULL'` is a length-1 literal against zero-length columns | guard `nrow(ks) == 0` with a message about one-condition sites, and use `rep('NULL', nrow(ks))` |
| `no rows to aggregate` from `KSEA.Scores` | the prior overlaps the site list in 0 or 1 place -- usually the abbreviated `data(KSData)`, a non-HUGO `Gene` column, or residues written `pS473`/`Ser473` | count `paste(Gene, Residue.Both)` against `paste(SUB_GENE, SUB_MOD_RSD)` before calling and stop with the coverage number |
| `Can't assign 4 names to a 0-column data.table` with `use_unmod_peptides = TRUE` | the evidence pre-filter removed every unmodified row | pass `use_unmod=TRUE` to `scripts/msstatsptm_labelfree.R` so unmodified rows are kept; the class-I rule applies to modified rows only |
| Site signs inverted in KSEA / up-down calls | Default pairwise Label `Control vs Treatment` means log2FC = Control - Treatment | Check `adjusted$Label` or pass an explicit `contrast.matrix` |
| Protein names without a site suffix in `ADJUSTED.Model` | Unmodified peptides from the enriched runs entered `$PTM` | Pre-filter evidence to modified rows; drop rows without `_<residue><position>` |
| KeyError / NaN on `Gene names` | Column is FASTA-dependent, absent without gene annotation | Guard with `.notna()` and fall back to `Protein` |
| All sites "regulated" and pathway-coherent | No protein-level adjustment | Require significance in MSstatsPTM `ADJUSTED.Model` |
| `PTM.Q.Value` / `PhosphoSite` not found (DIA-NN) | Those columns do not exist | Use `PTM.Site.Confidence` and `Site.Occupancy.Probabilities` |
| False "ubiquitination" sites | Iodoacetamide +114.0429 lysine artifact | Alkylate with chloroacetamide |
| Acetyl confused with trimethyl | +42.0106 vs +42.0470 isobaric at nominal mass | Require high-res MS; check 0.0364 Da split |

## References

- Beausoleil SA, Villen J, Gerber SA, Rush J, Gygi SP. A probability-based approach for high-throughput protein phosphorylation analysis and site localization. *Nat Biotechnol* 2006;24(10):1285-1292.
- Taus T, Kocher T, Pichler P, et al. Universal and confident phosphorylation site localization using phosphoRS. *J Proteome Res* 2011;10(12):5354-5362.
- Olsen JV, Blagoev B, Gnad F, et al. Global, in vivo, and site-specific phosphorylation dynamics in signaling networks. *Cell* 2006;127(3):635-648.
- Fermin D, Walmsley SJ, Gingras AC, Choi H, Nesvizhskii AI. LuciPHOr: algorithm for phosphorylation site localization with false localization rate estimation using modified target-decoy approach. *Mol Cell Proteomics* 2013;12(11):3409-3419.
- Bekker-Jensen DB, Bernhardt OM, Hogrebe A, et al. Rapid and site-specific deep phosphoproteome profiling by data-independent acquisition without the need for spectral libraries. *Nat Commun* 2020;11:787.
- Kohler D, Tsai TH, Verschueren E, et al. MSstatsPTM: Statistical Relative Quantification of Posttranslational Modifications in Bottom-Up Mass Spectrometry-Based Proteomics. *Mol Cell Proteomics* 2023;22(1):100477.
- Ochoa D, Jarnuczak AF, Vieitez C, et al. The functional landscape of the human phosphoproteome. *Nat Biotechnol* 2020;38(3):365-373.
- Casado P, Rodriguez-Prados JC, Cosulich SC, et al. Kinase-Substrate Enrichment Analysis Provides Insights into the Heterogeneity of Signaling Pathway Activation in Leukemia Cells. *Sci Signal* 2013;6(268):rs6.
- Krug K, Mertins P, Zhang B, et al. A Curated Resource for Phosphosite-specific Signature Analysis. *Mol Cell Proteomics* 2019;18(3):576-593.
- Yilmaz S, Ayati M, Schlatzer D, et al. Robust inference of kinase activity using functional networks. *Nat Commun* 2021;12:1177.
- Larsen MR, Thingholm TE, Jensen ON, Roepstorff P, Jorgensen TJD. Highly selective enrichment of phosphorylated peptides from peptide mixtures using titanium dioxide microcolumns. *Mol Cell Proteomics* 2005;4(7):873-886.
- Ruprecht B, Koch H, Medard G, et al. Comprehensive and reproducible phosphopeptide enrichment using iron immobilized metal ion affinity chromatography (Fe-IMAC) columns. *Mol Cell Proteomics* 2015;14(1):205-215.
- Matheron L, van den Toorn H, Heck AJR, Mohammed S. Characterization of biases in phosphopeptide enrichment by Ti(IV)-IMAC and TiO2 using a massive synthetic library and human cell digests. *Anal Chem* 2014;86(16):8312-8320.
- Thingholm TE, Jensen ON, Robinson PJ, Larsen MR. SIMAC (sequential elution from IMAC), a phosphoproteomics strategy for the rapid separation of monophosphorylated from multiply phosphorylated peptides. *Mol Cell Proteomics* 2008;7(4):661-671.
- Svinkina T, Gu H, Silva JC, et al. Deep, Quantitative Coverage of the Lysine Acetylome Using Novel Anti-acetyl-lysine Antibodies and an Optimized Proteomic Workflow. *Mol Cell Proteomics* 2015;14(9):2429-2440.
- Xu G, Paige JS, Jaffrey SR. Global analysis of lysine ubiquitination by ubiquitin remnant immunoaffinity profiling. *Nat Biotechnol* 2010;28(8):868-873.
- Kim W, Bennett EJ, Huttlin EL, et al. Systematic and Quantitative Assessment of the Ubiquitin-Modified Proteome. *Mol Cell* 2011;44(2):325-340.
- Nielsen ML, Vermeulen M, Bonaldi T, Cox J, Moroder L, Mann M. Iodoacetamide-induced artifact mimics ubiquitination in mass spectrometry. *Nat Methods* 2008;5(6):459-460.
- Akimov V, Barrio-Hernandez I, Hansen SVF, et al. UbiSite approach for comprehensive mapping of lysine and N-terminal ubiquitination sites. *Nat Struct Mol Biol* 2018;25(7):631-640.
- Riley NM, Bertozzi CR, Pitteri SJ. A Pragmatic Guide to Enrichment Strategies for Mass Spectrometry-Based Glycoproteomics. *Mol Cell Proteomics* 2021;20:100029.
- Johnson JL, Yaron TM, Huntsman EM, et al. An atlas of substrate specificities for the human serine/threonine kinome. *Nature* 2023;613(7945):759-766.
- Mueller-Dott S, Jaehnig EJ, et al. Comprehensive evaluation of phosphoproteomic-based kinase activity inference. *Nat Commun* 2025;16:4771.
- Zong Y, Wang Y, Yang Y, et al. DeepFLR facilitates false localization rate control in phosphoproteomics. *Nat Commun* 2023;14:2269.
- Ramsbottom KA, Prakash A, Riverol YP, et al. Method for independent estimation of the false localization rate for phosphoproteomics. *J Proteome Res* 2022;21(7):1603-1615.

## Related Skills

- peptide-identification - Identify modified peptides and run open/variable-mod search
- quantification - Underlying protein-level quant feeding the MSstatsPTM PROTEIN dataset
- differential-abundance - Moderated testing on the protein-level intensity matrix
- pathway-analysis/gsea - Enrichment scoring of regulated-site protein lists and PTM-SEA-style signatures
- data-visualization/sequence-logos - Render motif logos from the foreground/background windows
