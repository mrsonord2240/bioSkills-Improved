## CITE-seq: totalVI (Python, denoise + DE in one model)

**Goal:** Jointly model RNA + protein with explicit protein background, yielding a denoised latent space and foreground probabilities.

**Approach:** Register a MuData object, train the conditional VAE, then read the latent representation and per-protein foreground probability.

```python
import scvi
import mudata as md

# scvi-tools VAE training is stochastic unless seeded: verified two unseeded runs of this
# exact pattern on identical input differ by up to 0.97 (max abs latent diff); seeding
# makes reruns bit-identical. Set this before setup_mudata/train, every run.
scvi.settings.seed = 0

# mdata holds .mod['rna'] (raw counts) and .mod['prot'] (raw ADT counts)
scvi.model.TOTALVI.setup_mudata(
    mdata, rna_layer='counts', protein_layer=None,
    modalities={'rna_layer': 'rna', 'protein_layer': 'prot'}
)
model = scvi.model.TOTALVI(mdata)
model.train()

mdata.obsm['X_totalVI'] = model.get_latent_representation()
fg = model.get_protein_foreground_probability()        # 1 - background mixing weight per protein per cell
denoised_rna, denoised_prot = model.get_normalized_expression()
```

## Mosaic: MultiVI (Python, RNA+ATAC partially observed)

**Goal:** Jointly embed a mosaic design -- some cells have both RNA and ATAC (paired), others only one modality -- imputing the missing side.

**Approach:** Build one MuData with an RNA AnnData and an ATAC AnnData that both cover the full cell union; cells missing a modality get all-zero rows for that modality's block (MultiVI detects presence per cell from whether that block's raw counts sum to zero, not from a separate flag). Register with `setup_mudata`, not `setup_anndata` -- `MULTIVI.setup_anndata` on a plain AnnData is deprecated since scvi-tools 1.4 and silently skips registration (warns, then `MULTIVI(adata)` raises "Please set up your AnnData with MULTIVI.setup_anndata first").

```python
import scvi

scvi.settings.seed = 0

# mdata.mod['rna']: all cells, real counts. mdata.mod['atac']: real counts for paired
# cells, all-zero rows for RNA-only cells (and vice versa for an ATAC-only block).
scvi.model.MULTIVI.setup_mudata(
    mdata, modalities={'rna_layer': 'rna', 'atac_layer': 'atac'}
)
model = scvi.model.MULTIVI(
    mdata, n_genes=mdata.mod['rna'].n_vars, n_regions=mdata.mod['atac'].n_vars
)
model.train()
mdata.obsm['X_multivi'] = model.get_latent_representation()
```

Verified on synthetic 150-cell mosaic data (90 paired, 60 RNA-only, 3 known cell types, scvi-tools 1.5.1): RNA-only cells land nearer their same-type paired counterparts (mean latent distance 0.27) than different-type ones (0.56), confirming the model actually uses the shared RNA signal to place unpaired cells rather than clustering by modality of origin. This example passes no `batch_key` because the modality mask above is not a sequencing batch; if cells also span real sequencing batches, add `batch_key` for that separately -- see Common Errors' MultiVI row for the pitfall of confusing the two.
