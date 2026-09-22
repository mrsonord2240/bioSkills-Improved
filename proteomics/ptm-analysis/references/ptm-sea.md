# PTM-SEA With ssGSEA2.0

## PTM-SEA With ssGSEA2.0

**Goal:** Score PTMsigDB kinase, perturbation and pathway signatures from one signed site-level statistic.

**Approach:** PTM-SEA is ssGSEA run on site identifiers. Each row is ONE localized site, named by its +/-7 flanking sequence plus `-p` (`AALRQLRSPRRAQAP-p`). PTMsigDB scores each signature site with its reported direction (`;u` / `;d`, which never appear in your ids). Feed a SIGNED per-site statistic (log2FC or the moderated t from the protein-adjusted model), not intensities. Get ssGSEA2.0 from `github.com/broadinstitute/ssGSEA2.0` (checked on `cae7bed`); the human flanking database is `db/ptmsigdb/v1.9.1/ptm.sig.db.all.flanking.human.v1.9.1.gmt`, with mouse and rat files beside it.

```python
import pandas as pd

def write_ptmsea_gct(sites, path):
    # sites: one row per LOCALIZED site, with `seq_window` (MaxQuant "Sequence window": 31 aa, 15 each side)
    # and `stat` (signed statistic). PTMsigDB wants the central 15 (+/-7) plus '-p'.
    t = sites[sites['seq_window'].str.len() == 31].copy()
    t['id'] = t['seq_window'].str[8:23] + '-p'
    t = t.groupby('id', as_index=False)['stat'].mean()   # ids must be unique; averaging duplicates is a choice, report it
    with open(path, 'w', newline='\n') as f:
        f.write(f'#1.2\n{len(t)}\t1\nName\tDescription\tsample1\n')
        for i, v in zip(t['id'], t['stat']):
            f.write(f'{i}\t{i}\t{v:.6f}\n')
    return t
```

```bash
SSG=ssGSEA2.0
DB=$SSG/db/ptmsigdb/v1.9.1/ptm.sig.db.all.flanking.human.v1.9.1.gmt
mkdir -p ptmsea
Rscript $SSG/ssgsea-cli.R -i sites.gct -o ptmsea/run -d $DB -z $SSG \
  -n rank -w 0.75 -c z.score -t area.under.RES -s NES -p 1000 -m 10 -x TRUE -e FALSE -l FALSE
```

```python
def read_ptmsea(prefix='ptmsea/run'):
    gct = lambda kind: pd.read_csv(f'{prefix}-{kind}.gct', sep='\t', skiprows=2, index_col=0).iloc[:, -1]
    res = pd.DataFrame({'NES': gct('scores'), 'FDR': gct('fdr-pvalues')})
    return res.sort_values('NES', ascending=False)   # the run also writes -pvalues.gct and -combined.gct (signature size and overlap)
```

Checked on ssGSEA2.0 `cae7bed` with PTMsigDB v1.9.1, running the three blocks above as written, on a planted input: the 588 human CDK1 substrate sites of `KINASE-PSP_CDK1` given a +2 shift among 2,582 sites in all (real PTMsigDB flanking sequences, 2,000 of them from other signatures). `KINASE-PSP_CDK1` ranked 1 of 100 scored signatures (NES 44.0, FDR 0.0056, the floor at 1,000 permutations, which at least four other signatures also reached); CDK2, CDK5, CDK6 and ERK2/MAPK1 scored too because they share substrates. With the same statistics shuffled across sites, CDK1 fell to rank 37 (NES 0.47, FDR 0.95) and no signature had FDR below 0.05 (SERUM and PKCZ sat at exactly 0.05).

Only 100 of the 495 signatures had the `-m 10` overlap with those sites, so read the overlap column (`-combined.gct`) before calling a signature absent, and raise `-p` for finer p-values. PTMsigDB also holds non-phospho ids (`-ac`, `-m2`, `-m3`); the block above builds `-p` ids only. Signatures that share substrates co-score: name the substrate set, not an independent kinase (see "Over-reading kinase-activity output").
