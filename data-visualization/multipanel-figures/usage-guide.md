# Multi-Panel Figures - Usage Guide

## Overview

Use this Skill to compose two or more analytical plots as one publication figure. It supports patchwork, cowplot, and gridExtra in R plus matplotlib GridSpec, subfigures, and mosaics in Python.

The implementation details, version checks, journal-size table, export commands, and failure modes live in [SKILL.md](SKILL.md). In particular, a shared legend is only honest when the mappings and scales are identical, and a requested page size requires avoiding tight bounding-box export.

## Example prompts

- "Combine four identical-scale ggplots in a 2 by 2 patchwork grid, collect axes and one shared Condition guide, with lowercase 8 pt bold tags, at 183 by 140 mm."
- "Make a 183 mm Nature double-column matplotlib figure with a named mosaic: one large PCA panel and two smaller right-hand panels."
- "Use gridExtra with a layout matrix that makes the first panel span two cells; keep separate legends because the panels use different encodings."
- "Create two matplotlib subfigures: two stacked scatter plots on the left and a heatmap with its local colourbar on the right."

## Related Skills

- data-visualization/ggplot2-fundamentals
- data-visualization/matplotlib-fundamentals
- reporting/figure-export
