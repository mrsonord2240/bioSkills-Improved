---
name: bio-data-visualization-multipanel-figures
description: Compose multi-panel publication figures with patchwork, cowplot, gridExtra (R), or matplotlib GridSpec, subfigures, and subplot_mosaic (Python). Covers truthful shared-guide and shared-axis constraints, panel tags, exact journal dimensions, and vector-font PDF export.
tool_type: mixed
primary_tool: patchwork
license: MIT
author: GPTomics
---

# Multi-Panel Figures

Use this Skill to combine two or more finished panels into one submission-ready figure. Start by deciding the panel geometry, which encodings are genuinely shared, and the target journal dimensions. A visually similar layout is not enough: the exported PDF page must have the requested size and use vector text.

## Version and export checks

The examples were checked with patchwork 1.3.2, ggplot2 4.0.3, cowplot 1.2.0, gridExtra 2.3, Cairo 1.7.0, and matplotlib 3.11.2. `axes = "collect"` was introduced in patchwork 1.2.0 (CRAN release 2024-01-08); use patchwork >= 1.3 with ggplot2 4.x. Older patchwork versions reject `axes` with an unused-argument error rather than silently collecting axes.

Before adapting a recipe, check `packageVersion("patchwork")` / `packageVersion("ggplot2")` or `matplotlib.__version__`. For PDF inspection, use `pdfinfo figure.pdf` for page dimensions and `pdffonts figure.pdf` to check for embedded TrueType/CID TrueType rather than Type 3 fonts.

### Nature export constants

| item | Nature setting used here |
| --- | --- |
| single column | 89 mm wide |
| double column | 183 mm wide |
| example canvas | 183 x 140 mm |
| body text | 7 pt |
| panel tag | 8 pt, bold |
| R PDF device | `Cairo::CairoPDF` |
| Python PDF font | `matplotlib.rcParams["pdf.fonttype"] = 42` |

Use the target journal's own current specifications for other journals. Do not use `bbox_inches="tight"` when the page dimensions are contractual: it changes the saved page bounding box. Keep labels inside each axes instead.

## Choose a layout

| need | R | Python |
| --- | --- | --- |
| regular grid with identical scales | `wrap_plots()` plus collection | `plt.subplots(..., layout="constrained")` |
| unequal spans | patchwork design string | `GridSpec` or `subplot_mosaic` |
| nested independent regions | patchwork nesting | `fig.subfigures()` |
| precise base-grid placement | `gridExtra::arrangeGrob(layout_matrix=)` | `GridSpec` |
| one legend | only genuinely identical mappings/scales | `fig.legend()` from a representative axes |

## patchwork: flat grid, shared axes and a shared legend

`axes = "collect"` and `axis_titles = "collect"` work on a flat `wrap_plots()` composition when panels have compatible, identical coordinate scales. They do not transform different variables, limits, or scale types into a common axis. `guides = "collect"` only merges equivalent guides: same aesthetic, labels, scale name, limits/breaks, and palette. Do not suppress N-1 legends when their mappings differ; that would imply equivalence that is not present.

```r
library(ggplot2)
library(patchwork)

set.seed(1)
df <- data.frame(x = rnorm(160), y = rnorm(160),
                 group = rep(c("Control", "Treatment"), each = 80))
theme_publication <- theme_classic(base_size = 7) +
  theme(plot.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

palette <- c(Control = "#4DBBD5", Treatment = "#E64B35")
make_panel <- function(dat, title) {
  ggplot(dat, aes(x, y, colour = group)) +
    geom_point(size = 1.2, alpha = 0.75) +
    scale_colour_manual(values = palette, name = "Condition") +
    coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) +
    labs(x = "Standardized x", y = "Standardized y", title = title) +
    theme_publication
}

# All four panels have the same coordinate limits and the same colour scale.
fig <- wrap_plots(lapply(letters[1:4], \(tag) make_panel(df, paste("Panel", tag))), ncol = 2) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"),
        plot.tag.location = "panel", plot.tag.position = c(0.02, 0.98),
        legend.position = "bottom")

save_cairo_pdf <- function(path, plot, width_mm = 183, height_mm = 140) {
  Cairo::CairoPDF(path, width = width_mm / 25.4, height = height_mm / 25.4)
  on.exit(dev.off(), add = TRUE)
  print(plot)
}
save_cairo_pdf("figure.pdf", fig)
```

The `& theme(...)` is intentional: styling tags inside `plot_annotation(theme = ...)` does not reliably propagate to patchwork tags. Tags above are positioned inside each panel to preserve the requested page size.

For a nested composition, collect axes at the nesting level that owns the shared boundary, or prefer the flat `wrap_plots()` form. Confirm collection in the rendered figure; it is a layout convenience, not semantic proof that axes are comparable.

## Other R arrangements

Use cowplot when alignment rather than guide collection is the priority. Its `align` result depends on the requested `axis` sides and the grobs' margins, so do not promise that `align = "v"` fails for different axis-label widths; inspect the output.

```r
library(cowplot)
library(ggplot2)
set.seed(1)
d <- data.frame(x = rnorm(80), y = rnorm(80))
panels <- lapply(letters[1:4], \(tag) ggplot(d, aes(x, y)) + geom_point() +
                  labs(title = tag) + theme_classic(base_size = 7))
combined <- plot_grid(plotlist = panels, ncol = 2, labels = "auto",
                      label_size = 8, label_fontface = "bold", align = "hv")
save_cairo_pdf <- function(path, plot, width_mm = 183, height_mm = 140) {
  Cairo::CairoPDF(path, width = width_mm / 25.4, height = height_mm / 25.4)
  on.exit(dev.off(), add = TRUE)
  print(plot)
}
save_cairo_pdf("cowplot.pdf", combined)
```

Use `gridExtra` for an explicit layout matrix. It arranges grobs but does not collect guides; construct one shared legend separately only when the panel scales are identical.

```r
library(gridExtra)
library(ggplot2)
set.seed(1)
d <- data.frame(x = rnorm(80), y = rnorm(80))
panels <- lapply(letters[1:4], \(tag) ggplot(d, aes(x, y)) + geom_point() +
                  labs(title = tag) + theme_classic(base_size = 7))
layout <- rbind(c(1, 1, 2), c(3, 4, 4))
g <- do.call(arrangeGrob, c(panels, list(layout_matrix = layout)))
save_cairo_pdf <- function(path, plot, width_mm = 183, height_mm = 140) {
  Cairo::CairoPDF(path, width = width_mm / 25.4, height = height_mm / 25.4)
  on.exit(dev.off(), add = TRUE)
  grid::grid.draw(plot)
}
save_cairo_pdf("gridextra.pdf", g)
```

## matplotlib: exact-size grid and one qualified shared legend

Set the font type before saving, enable constrained layout explicitly, and place tags inside the axes with a fixed point offset. `fig.legend()` does not merge arbitrary legends: reuse it only when panels use the same labels and the same mapping/palette.

```python
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
import numpy as np

width_mm, height_mm, dpi = 183, 140, 300
fig, axs = plt.subplots(2, 2, figsize=(width_mm / 25.4, height_mm / 25.4),
                        dpi=dpi, layout="constrained", sharex=True, sharey=True)
palette = {"Control": "#4DBBD5", "Treatment": "#E64B35"}
rng = np.random.default_rng(1)
x = rng.normal(size=160)
y = rng.normal(size=160)
groups = np.repeat(np.array(["Control", "Treatment"]), 80)
for ax, label in zip(axs.flat, "abcd"):
    for group, colour in palette.items():
        mask = groups == group
        ax.scatter(x[mask], y[mask], s=10, color=colour, label=group, rasterized=True)
    ax.set(xlim=(-3, 3), ylim=(-3, 3), xlabel="Standardized x", ylabel="Standardized y")
    tag_transform = offset_copy(ax.transAxes, fig=fig, x=5, y=-5, units="points")
    ax.text(0, 1, label, transform=tag_transform, fontsize=8, fontweight="bold",
            va="top", ha="left")

handles, labels = axs.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, title="Condition", loc="outside lower center", ncols=2)
fig.savefig("figure.pdf")
fig.savefig("figure.png", dpi=dpi)
```

`sharex=True` and `sharey=True` share locator/limit state; use them only when the plotted quantities and units are comparable. Matplotlib's constrained layout is opt-in (via `layout="constrained"` or `constrained_layout=True`), not a global default.

### GridSpec, subfigures, and mosaic

```python
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
fig = plt.figure(figsize=(183 / 25.4, 140 / 25.4), layout="constrained")
gs = GridSpec(2, 3, figure=fig)
ax_left = fig.add_subplot(gs[:, :2])
ax_top = fig.add_subplot(gs[0, 2])
ax_bottom = fig.add_subplot(gs[1, 2], sharex=ax_top)

# A named, spanning layout; each label identifies a single rectangular axes.
fig2, axd = plt.subplot_mosaic([["A", "A", "B"], ["A", "A", "C"]],
                               figsize=(183 / 25.4, 140 / 25.4), layout="constrained")
fig.savefig("gridspec.pdf")
fig2.savefig("mosaic.pdf")
```

For a nested region with its own layout engine, use this self-contained pattern; the colourbar remains local to the heatmap subfigure.

```python
import matplotlib.pyplot as plt
import numpy as np
fig = plt.figure(figsize=(183 / 25.4, 140 / 25.4), layout="constrained")
left, right = fig.subfigures(1, 2, width_ratios=(2, 1))
axes = left.subplots(2, 1, sharex=True)
axes[0].plot(np.arange(10)); axes[1].scatter(np.arange(10), np.arange(10))
heatmap = right.subplots().imshow(np.arange(100).reshape(10, 10))
right.colorbar(heatmap, ax=right.axes[0])
fig.savefig("subfigures.pdf")
```

## Failure modes and checks

| symptom | cause | correction |
| --- | --- | --- |
| tag style remains regular/default sized | tag theme was passed inside `plot_annotation()` | apply `& theme(plot.tag = element_text(size=8, face="bold"))` |
| repeated axes remain | nested patchwork layout or incompatible scales | use flat `wrap_plots()` for compatible panels, or collect at the owning nesting level |
| two legends remain | aesthetics/scales are not equivalent | keep both, or make the mapping, palette, title, breaks, and limits identical |
| R page is unexpectedly enormous | `units` omitted; `ggsave()` normally errors above its size limit | set `units="mm"`; never disable the guard to export an accidental inch canvas |
| Python page is larger than requested | `bbox_inches="tight"` changed the output bounding box | omit it and keep text/tags inside the canvas |
| Type 3 fonts in Python PDF | default `pdf.fonttype` used | set `matplotlib.rcParams["pdf.fonttype"] = 42` before save |
| typography/font acceptance is uncertain | PDF backend and installed fonts differ | inspect `pdffonts`; `Cairo::CairoPDF` and Type 42 are tested paths, not a blanket guarantee |

After export, inspect page size and fonts, then render the PDF to a PNG to detect clipped labels. The bundled examples export 183 x 140 mm at 300 dpi as 2161 x 1653 pixels; calculate expected pixels from the configured mm and dpi for a different backend or canvas.

## Related Skills

- data-visualization/ggplot2-fundamentals - individual R panels
- data-visualization/matplotlib-fundamentals - individual Python panels
- reporting/figure-export - journal-format checks
- data-visualization/color-palettes - consistent encodings
