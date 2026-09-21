# GNINA reference (moved from SKILL.md)

## GNINA with CNN Scoring (modern default)

```bash
gnina -r receptor.pdb -l ligand.sdf \
      --autobox_ligand reference_ligand.sdf \
      --cnn_scoring rescore \
      -o poses.sdf.gz \
      --num_modes 9 --exhaustiveness 8
```

`--cnn_scoring`:
- `none`: no CNN; use the selected empirical scoring function throughout
- `rescore` (default): use empirical scoring during the search, then CNN-rerank the final poses; least computationally expensive CNN option
- `refinement`: use the CNN to refine poses after Monte Carlo chains and to rank the final poses; approximately 10 times slower than `rescore` on a GPU in the official documentation
- `metrorescore`: use CNN scoring in the Metropolis search and rescore the resulting poses
- `metrorefine`: use CNN scoring in the Metropolis search and refine the resulting poses
- `all`: use the CNN scoring function throughout; the official documentation describes this as extremely computationally intensive and not recommended

The six choices above are from GNINA 1.3. Earlier releases expose a smaller set; check `gnina --help` for the installed executable rather than assuming every mode is available.

`--autobox_ligand`: define box from reference ligand SDF/PDB. Otherwise specify `--center_x/y/z` + `--size_x/y/z`.

**Critical:** GNINA distributions include multiple named CNN models/ensembles rather than one universally described "PDBbind 2019" model. Record the selected model or ensemble and validate it with known co-crystal redocking and, when relevant, cross-docking controls.
