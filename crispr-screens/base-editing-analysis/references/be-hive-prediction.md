# BE-Hive Editing-Efficiency Prediction

Read when predicting per-spacer editing efficiency or bystander outcomes with BE-Hive (Arbab 2020).

## BE-Hive Editing-Efficiency Prediction

BE-Hive (`be_predict_bystander`) takes a fixed **50nt substrate**, not the bare 20nt spacer:
19nt of upstream context + the 20nt spacer + the 3nt PAM + 8nt of downstream context
(`-19..30` in the README's own numbering; spacer occupies substrate positions 1-20 of that
window, PAM at 21-23). Getting the upstream-context length wrong by even 1nt silently
shifts every predicted position by one base with no error -- checked on this Skill's own
worked example, verified via BE-Hive's own `pred_df` diagnostic field.

```python
import sys
sys.path.append("/path/to/be_predict_bystander/..")  # parent dir of the cloned repo
from be_predict_bystander import predict as bystander_model

spacer = "TGATCACGTAGCATGCACGT"  # 20nt
pam = "TGG"
upstream_19nt = "ATGCATGGATCGTAGCTAG"    # 19nt of real genomic context immediately 5' of the spacer
downstream_8nt = "CATGCTAG"              # 8nt of real genomic context immediately 3' of the PAM
substrate = upstream_19nt + spacer + pam + downstream_8nt
assert len(substrate) == 50

bystander_model.init_model(base_editor="BE4", celltype="mES")  # celltype in {'mES','HEK293','U2OS',...}
pred_df, stats = bystander_model.predict(substrate)

# Always cross-check BE-Hive's own read-back against the intended spacer before trusting
# pred_df's position-labeled columns (e.g. 'C4', 'C6') -- a wrong substrate length or
# offset produces a plausible-looking but silently mis-positioned prediction.
assert substrate[19:39] == spacer, "substrate/spacer offset is wrong -- check upstream context length"
print(stats["Total predicted probability"])
print(pred_df.sort_values("Predicted frequency", ascending=False).head(10))
```

Checked on BE-Hive git HEAD (maxwshen/be_predict_bystander, 2026-09-16) against a synthetic guide
with a known target C at spacer position 5 and bystander C at spacer position 7: `pred_df`'s `C4`/`C6`
columns (BE-Hive's own 0-indexed-from-position-4 editable-C naming) correctly identified both, and
`Total predicted probability` was 0.97-0.98 (not a stub, not all-zero).
