# PTM-SEA With ssGSEA2.0

**Goal:** Score PTMsigDB kinase, perturbation and pathway signatures from one signed site-level statistic.

**Approach:** PTM-SEA is ssGSEA run on site identifiers. Each row is ONE localized site, named by its +/-7 flanking sequence plus `-p` (`AALRQLRSPRRAQAP-p`). PTMsigDB scores each signature site with its reported direction (`;u` / `;d`, which never appear in your ids). Feed a SIGNED per-site statistic (log2FC or the moderated t from the protein-adjusted model), not intensities. Get ssGSEA2.0 from `github.com/broadinstitute/ssGSEA2.0` (checked on `cae7bed`); the human flanking database is `db/ptmsigdb/v1.9.1/ptm.sig.db.all.flanking.human.v1.9.1.gmt`, with mouse and rat files beside it.

```bash
python scripts/ptmsea.py write-gct --sites sites.tsv --out sites.gct
```

`sites.tsv` has one row per LOCALIZED site with `seq_window` (MaxQuant "Sequence window": 31 aa, 15 each side) and `stat` (the signed statistic). The script keeps the central 15 (+/-7) plus `-p` as the id; ids must be unique, so duplicates are averaged (a choice: report it).

```bash
SSG=ssGSEA2.0
DB=$SSG/db/ptmsigdb/v1.9.1/ptm.sig.db.all.flanking.human.v1.9.1.gmt
mkdir -p ptmsea
Rscript $SSG/ssgsea-cli.R -i sites.gct -o ptmsea/run -d $DB -z $SSG \
  -n rank -w 0.75 -c z.score -t area.under.RES -s NES -p 1000 -m 10 -x TRUE -e FALSE -l FALSE
```

```bash
python scripts/ptmsea.py read --prefix ptmsea/run --out ptmsea_scores.csv   # NES and FDR per signature, sorted by NES
```

The run also writes `-pvalues.gct` and `-combined.gct` (signature size and overlap).

Checked on ssGSEA2.0 `cae7bed` with PTMsigDB v1.9.1, running the three steps above (`write-gct`, the ssgsea-cli call, `read`; originally as inline blocks, re-run 2026-09-21 as `scripts/ptmsea.py`) on a planted input: the 588 human CDK1 substrate sites of `KINASE-PSP_CDK1` given a +2 shift among 2,582 sites in all (real PTMsigDB flanking sequences, 2,000 of them from other signatures). `KINASE-PSP_CDK1` ranked 1 of 100 scored signatures (NES 44.0, FDR 0.0056, the floor at 1,000 permutations, which at least four other signatures also reached); CDK2, CDK5, CDK6 and ERK2/MAPK1 scored too because they share substrates. With the same statistics shuffled across sites, CDK1 fell to rank 37 (NES 0.47, FDR 0.95) and no signature had FDR below 0.05 (SERUM and PKCZ sat at exactly 0.05).

Only 100 of the 495 signatures had the `-m 10` overlap with those sites, so read the overlap column (`-combined.gct`) before calling a signature absent, and raise `-p` for finer p-values. PTMsigDB also holds non-phospho ids (`-ac`, `-m2`, `-m3`); `write-gct` builds `-p` ids only. Signatures that share substrates co-score: name the substrate set, not an independent kinase (see "Over-reading kinase-activity output" in SKILL.md).
