## Factor-Based Analysis

For complex perturbation responses, decompose the per-cell perturbation effect into shared latent factors:

```python
import pertpy as pt

# FR-Perturb ("Factorize-Recover") decomposes perturbation effects into shared factors.
# It is NOT part of pertpy: it is a standalone CLI from douglasyao/FR-Perturb
# (Yao et al. 2023 Nat Biotechnol). Run it outside Python:
#   python run_FR_Perturb.py --input <expression> --perturbations <matrix> --out <prefix>
```
