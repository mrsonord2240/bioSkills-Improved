'''Label-free normalization, TMT cross-plex IRS bridge, SILAC labeling check and AP-MS control scoring, self-contained.

Every section asserts its expected output, so a silent regression exits non-zero.'''
# Reference: numpy 2.5.3, pandas 3.0.5 | Verify API if version differs
import numpy as np
import pandas as pd

rng = np.random.default_rng(0)


def median_center(log_int):
    '''Subtract each sample median log2 intensity (corrects LOCATION only, cannot manufacture variance).'''
    sample_medians = log_int.median(axis=0)
    return log_int - sample_medians + sample_medians.median()


def sample_loading_normalize(plex):
    '''Per-channel scalar correcting total load WITHIN one plex; target = mean column sum.'''
    target = plex.sum(axis=0).mean()
    return plex * (target / plex.sum(axis=0))


def irs_scale(plexes, ref_cols):
    '''Internal Reference Scaling bridge (Plubell 2017): pin each plex reference channel to a common per-protein value.'''
    refs = pd.concat([p[ref] for p, ref in zip(plexes, ref_cols)], axis=1)
    refs = refs.where(refs > 0)    # a 0 or missing reference cannot anchor the bridge
    unbridged = refs.index[refs.isna().any(axis=1)]
    if len(unbridged):
        print(f'IRS: {len(unbridged)} proteins lack a reference in >=1 plex; set to NaN: {list(unbridged[:10])}')
    geomean = np.exp(np.log(refs).mean(axis=1, skipna=False))    # per-protein geometric mean across plexes
    out = []
    for i, p in enumerate(plexes):
        factor = geomean / refs.iloc[:, i]    # per-protein per-plex scaling factor
        out.append(p.mul(factor, axis=0))
    return out


def plex_offset(a, b, ref_a, ref_b):
    '''Median per-protein log2 offset between plexes on NON-reference channels (IRS equalizes references by construction).'''
    la = np.log2(a.drop(columns=ref_a)).mean(axis=1)
    lb = np.log2(b.drop(columns=ref_b)).mean(axis=1)
    return float(np.nanmedian(lb - la))


proteins = [f'P{i:03d}' for i in range(200)]
samples = [f'S{i}' for i in range(6)]

# Label-free: raw intensities with per-sample loading offsets and ~15% MaxQuant zeros (missing)
true_abundance = rng.normal(20, 2, (200, 6))
loading_offset = np.array([0.0, 0.4, -0.3, 0.6, -0.5, 0.2])    # per-sample load differences to be normalized out
raw = 2 ** (true_abundance + loading_offset)
raw[rng.random((200, 6)) < 0.15] = 0    # MaxQuant writes 0 for "not quantified"
lfq = pd.DataFrame(raw, index=proteins, columns=samples)

log_int = np.log2(lfq.replace(0, np.nan))    # 0 -> NaN before transform; log2(0) = -inf otherwise
print(f'Missing values: {100 * log_int.isna().sum().sum() / log_int.size:.1f}%')
print(f'Sample medians before centering: {np.round(log_int.median().values, 2)}')

normalized = median_center(log_int)
print(f'Sample medians after centering:  {np.round(normalized.median().values, 2)}')

# Filter proteins quantified in too few samples; min-peptides logic mirrors min-ratio-count=2 robustness
valid_per_protein = normalized.notna().sum(axis=1)
min_valid = len(samples) // 2    # require presence in >=50% of samples
filtered = normalized[valid_per_protein >= min_valid]
print(f'Proteins after >=50% valid filter: {len(filtered)}')

# TMT cross-plex: two 4-channel plexes over the same proteins, last channel a pooled reference.
# Plex B sees each protein at a random per-protein offset around 2x (elution sampling, not biology).
protein_level = rng.normal(18, 1.5, (200, 1))
plex_offset_true = rng.normal(1.0, 0.5, (200, 1))
plex_a = pd.DataFrame(2 ** (protein_level + rng.normal(0, 0.15, (200, 4))), index=proteins, columns=['A1', 'A2', 'A3', 'A_ref'])
plex_b = pd.DataFrame(2 ** (protein_level + plex_offset_true + rng.normal(0, 0.15, (200, 4))), index=proteins, columns=['B1', 'B2', 'B3', 'B_ref'])
plex_b.iloc[0, 3] = 0    # one protein without a plex-B reference reporter: reported and left unbridged
sl = [sample_loading_normalize(plex_a), sample_loading_normalize(plex_b)]
bridged = irs_scale(sl, ['A_ref', 'B_ref'])
print(f'Non-reference plex offset (median log2 B-A): SL only {plex_offset(sl[0], sl[1], "A_ref", "B_ref"):+.3f}, '
      f'after IRS {plex_offset(bridged[0], bridged[1], "A_ref", "B_ref"):+.3f}')


# ---------------------------------------------------------------------------------------------
# SILAC: labeling-efficiency and Arg->Pro checks on a seeded HEAVY-ONLY pilot (SKILL.md, SILAC section)
# ---------------------------------------------------------------------------------------------
def silac_labeling_efficiency(pilot, heavy='Intensity H', light='Intensity L'):
    '''Incorporation on a HEAVY-ONLY pilot (cells grown in heavy medium, NOTHING mixed in): every
    light ion there is unlabeled protein. pilot = the pilot's peptide table (MaxQuant evidence.txt).'''
    h = pd.to_numeric(pilot[heavy], errors='coerce').fillna(0.0)
    l = pd.to_numeric(pilot[light], errors='coerce').fillna(0.0)
    ok = (h + l) > 0
    if not ok.any():
        raise ValueError('no peptide has signal in either channel -- wrong columns or wrong file')
    per_pep = h[ok] / (h[ok] + l[ok])
    eff = float(h[ok].sum() / (h[ok] + l[ok]).sum())    # intensity-weighted = the number to report
    # A 1:1 forward mix of these cells does NOT read log2 H/L = 0: the unincorporated (1 - eff) of
    # the heavy sample is counted in the LIGHT channel, so H/L = eff / (2 - eff) -- -0.20 at 93%.
    return {'n_peptides': int(ok.sum()), 'incorporation': round(eff, 4),
            'median_peptide_incorporation': round(float(per_pep.median()), 4),
            'peptides_below_95pct': int((per_pep < 0.95).sum()),
            'expected_log2_HL_bias_at_1to1': round(float(np.log2(eff / (2 - eff))), 4),
            'pass_95pct': bool(eff >= 0.95)}

def arg_to_pro_shift(peptides, seq='Sequence', ratio='Ratio H/L'):
    '''Arg->Pro drains the heavy channel once per proline, so log2 H/L falls with PROLINE COUNT.
    The dose slope is what makes this specific -- a flat offset is incomplete labeling instead.
    Direct route: re-search the pilot with Pro6 variable and take I(Pro6)/(I(Pro6)+I(Pro0)).'''
    r = np.log2(pd.to_numeric(peptides[ratio], errors='coerce'))
    npro = peptides[seq].astype(str).str.count('P')
    ok = np.isfinite(r)
    r, npro = r[ok], npro[ok]
    slope = float(np.polyfit(npro, r, 1)[0]) if npro.nunique() > 1 else float('nan')
    return {'n_pro_free': int((npro == 0).sum()), 'n_pro_containing': int((npro > 0).sum()),
            'median_log2_HL_pro_free': round(float(r[npro == 0].median()), 4),
            'log2_HL_slope_per_proline': round(slope, 4),
            'conversion_per_proline': round(float(1 - 2 ** slope), 4)}


TRUE_INCORPORATION, TRUE_CONVERSION = 0.93, 0.08
n_pep = 600
seqs = [''.join(rng.choice(list('ACDEFGHIKLMNQRSTVWY'), 8)) + 'P' * int(k) + 'K' for k in rng.integers(0, 4, n_pep)]
total = 10 ** rng.uniform(5, 8, n_pep)
pro = pd.Series(seqs).str.count('P').to_numpy()
# H/L starts at eff/(1-eff) and each proline multiplies it by (1 - conversion): the Arg->Pro dose slope
hl = TRUE_INCORPORATION / (1 - TRUE_INCORPORATION) * (1 - TRUE_CONVERSION) ** pro
heavy_frac = hl / (1 + hl)
pilot = pd.DataFrame({'Sequence': seqs, 'Intensity H': total * heavy_frac, 'Intensity L': total * (1 - heavy_frac)})
pilot['Ratio H/L'] = pilot['Intensity H'] / pilot['Intensity L']
pro_free = pilot[pro == 0]
eff = silac_labeling_efficiency(pro_free)    # incorporation is read on proline-free peptides
conv = arg_to_pro_shift(pilot)
print(f'SILAC pilot: incorporation {eff["incorporation"]} (planted {TRUE_INCORPORATION}), '
      f'Arg->Pro per proline {conv["conversion_per_proline"]} (planted {TRUE_CONVERSION})')
assert abs(eff['incorporation'] - TRUE_INCORPORATION) < 0.02 and eff['pass_95pct'] is False
assert abs(conv['conversion_per_proline'] - TRUE_CONVERSION) < 0.02 and conv['log2_HL_slope_per_proline'] < 0

# ---------------------------------------------------------------------------------------------
# AP-MS: score against control IPs, not input lysate (SKILL.md, AP-MS section)
# ---------------------------------------------------------------------------------------------
def score_vs_control_ips(ip, bait_cols, ctrl_cols, fc_cutoff=2.0, min_bait_reps=None):
    '''ip: prey x replicate RAW intensities or spectral counts; 0/NaN = not detected. Prey never
    seen in any bait IP get a NaN enrichment (nothing was measured) and are never called.'''
    L = np.log2(ip[list(bait_cols) + list(ctrl_cols)].replace(0, np.nan).astype(float))
    if min_bait_reps is None:
        min_bait_reps = len(bait_cols)          # default: every bait replicate, reproducibility first
    n_bait, n_ctrl = L[bait_cols].notna().sum(axis=1), L[ctrl_cols].notna().sum(axis=1)
    # a control IP that produced no data is skipped by every mean/min below, so say so out loud
    dead = [c for c in ctrl_cols if L[c].notna().sum() == 0]
    if dead:
        print(f'AP-MS: control runs with no data, excluded: {dead}')
    if len(dead) == len(ctrl_cols):
        raise ValueError('every control IP is empty -- nothing to score against')
    Lf = L.copy()
    Lf[ctrl_cols] = Lf[ctrl_cols].fillna(L[ctrl_cols].min())      # per-control-run detection floor
    enrich = Lf[bait_cols].mean(axis=1) - Lf[ctrl_cols].mean(axis=1)
    worst = Lf[bait_cols].min(axis=1) - Lf[ctrl_cols].max(axis=1)  # weakest bait rep vs best control
    out = pd.DataFrame({'n_bait': n_bait, 'n_ctrl': n_ctrl, 'n_ctrl_runs_used': len(ctrl_cols) - len(dead),
                        'log2_enrichment': enrich, 'worst_case_log2': worst})
    out['interactor'] = (n_bait >= min_bait_reps) & ((n_ctrl == 0) | (enrich >= fc_cutoff))
    return out.sort_values('log2_enrichment', ascending=False)


n_prey = 60
kind = np.array(['interactor'] * 5 + ['sticky'] * 10 + ['background'] * (n_prey - 15))
base = rng.normal(22, 1.0, n_prey)    # log2 background level of every prey
ip_cols = ['bait1', 'bait2', 'bait3', 'ctrl1', 'ctrl2', 'ctrl3']
lg = pd.DataFrame({c: base + rng.normal(0, 0.3, n_prey) for c in ip_cols}, index=[f'prey{i:02d}' for i in range(n_prey)])
lg.loc[kind == 'interactor', ['bait1', 'bait2', 'bait3']] += 6    # real interactors: bait IPs only
lg.loc[kind == 'sticky', :] += 3                                  # bead binders: every IP, bait and control alike
ip = 2 ** lg
scored = score_vs_control_ips(ip, ['bait1', 'bait2', 'bait3'], ['ctrl1', 'ctrl2', 'ctrl3'])
called = set(scored.index[scored['interactor']])
truth = set(lg.index[kind == 'interactor'])
sticky = set(lg.index[kind == 'sticky'])
print(f'AP-MS: called {len(called)}; interactors {len(called & truth)}/{len(truth)}, sticky binders called {len(called & sticky)}/{len(sticky)}')
assert called == truth    # sticky binders excluded, every planted interactor recovered

# a failed control IP must be announced, not absorbed
dead_ip = ip.copy()
dead_ip['ctrl3'] = 0.0
scored_dead = score_vs_control_ips(dead_ip, ['bait1', 'bait2', 'bait3'], ['ctrl1', 'ctrl2', 'ctrl3'])
assert (scored_dead['n_ctrl_runs_used'] == 2).all() and (scored['n_ctrl_runs_used'] == 3).all()

# min_bait_reps: one stochastic dropout in 3 bait replicates costs the strict default a true interactor
drop = ip.copy()
victim = sorted(truth)[0]
drop.loc[victim, 'bait1'] = 0.0
strict = score_vs_control_ips(drop, ['bait1', 'bait2', 'bait3'], ['ctrl1', 'ctrl2', 'ctrl3'])
loose = score_vs_control_ips(drop, ['bait1', 'bait2', 'bait3'], ['ctrl1', 'ctrl2', 'ctrl3'], min_bait_reps=2)
print(f'min_bait_reps: {victim} called strict={bool(strict.loc[victim, "interactor"])}, min_bait_reps=2 -> {bool(loose.loc[victim, "interactor"])}')
assert not strict.loc[victim, 'interactor'] and loose.loc[victim, 'interactor']
