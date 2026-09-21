# S-PrediXcan + S-MultiXcan Pipeline

**Goal:** Run TWAS across all GTEx tissues using pre-trained PredictDB models, then combine via S-MultiXcan for a joint multi-tissue test.

**Approach:** For each tissue, run S-PrediXcan with the matched model DB and covariance file; collect per-tissue outputs into a folder; run S-MultiXcan with the same model folder and GWAS to produce a joint multi-tissue Z and per-tissue significance.

```bash
# PredictDB models: predictdb.org
# GTEx v8 MASHR-EUR is the standard EUR panel
# Each tissue has a .db (model) and .txt.gz (covariance) file pair

mkdir -p spredixcan_out
for tissue in Whole_Blood Liver Brain_Frontal_Cortex_BA9 Adipose_Subcutaneous; do
    python SPrediXcan.py \
        --model_db_path mashr_models/mashr_${tissue}.db \
        --covariance mashr_models/mashr_${tissue}.txt.gz \
        --gwas_file gwas.txt \
        --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 \
        --beta_column BETA --pvalue_column P \
        --output_file spredixcan_out/${tissue}.csv
done

# Joint multi-tissue
python SMulTiXcan.py \
    --models_folder mashr_models/ \
    --models_name_pattern "mashr_(.*)\.db" \
    --snp_covariance gtex_v8_expression_mashr_snp_covariance.txt.gz \
    --metaxcan_folder spredixcan_out/ \
    --metaxcan_filter "(.*)\.csv" \
    --metaxcan_file_name_parse_pattern "(.*)\.csv" \
    --gwas_file gwas.txt \
    --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 \
    --beta_column BETA --pvalue_column P \
    --cutoff_condition_number 30 \
    --output joint_multitissue.csv
```

`--cutoff_condition_number 30` (the canonical MetaXcan setting) is required, not optional; omitting it raises the `InvalidArguments` error in SKILL.md Common Errors (confirmed 2026-09-19 against a real run).

S-MultiXcan applies PCA regularisation on the inter-tissue correlation matrix: `--cutoff_condition_number 30` drops near-collinear components, and `--regularization 0.1` (off unless passed explicitly) adds a ridge. Tissues that are nearly collinear with another (e.g. multiple brain sub-regions) are absorbed into shared components and do not contribute independent power.
