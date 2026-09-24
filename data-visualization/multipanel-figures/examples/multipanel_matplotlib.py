"""Create exact-size 183 x 140 mm multi-panel figures with matplotlib.

Run: python examples/multipanel_matplotlib.py
The four panels deliberately have identical x/y quantities and a shared palette;
only under that condition is one figure-level legend and shared axes truthful.
"""

from pathlib import Path

import matplotlib

matplotlib.rcParams["pdf.fonttype"] = 42

import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
import numpy as np
from PIL import Image

WIDTH_MM, HEIGHT_MM, DPI = 183, 140, 300
PALETTE = {"Control": "#4DBBD5", "Treatment": "#E64B35"}
OUT = Path(".")


def add_fixed_point_tag(ax, label, fig):
    """Place an 8 pt tag five points inside its axes' upper-left corner."""
    transform = offset_copy(ax.transAxes, fig=fig, x=5, y=-5, units="points")
    ax.text(0, 1, label, transform=transform, fontsize=8, fontweight="bold",
            va="top", ha="left")


def add_identical_mapping(ax, x, y, groups, label, fig):
    for group, colour in PALETTE.items():
        mask = groups == group
        ax.scatter(x[mask], y[mask], color=colour, label=group,
                   s=10, alpha=0.75, rasterized=True)
    ax.set(xlim=(-3.2, 3.2), ylim=(-3.2, 3.2), xlabel="Standardized x",
           ylabel="Standardized y", title=f"Panel {label}")
    add_fixed_point_tag(ax, label.lower(), fig)


rng = np.random.default_rng(42)
x = rng.normal(size=240)
groups = np.repeat(np.array(["Control", "Treatment"]), 120)
y = rng.normal(size=240) + np.where(groups == "Treatment", 0.35, -0.35)

# Regular 2 x 2 layout: the physical canvas is exact because no tight bbox is used.
fig, axs = plt.subplots(
    2, 2, figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), dpi=DPI,
    layout="constrained", sharex=True, sharey=True,
)
for index, (ax, label) in enumerate(zip(axs.flat, "ABCD")):
    add_identical_mapping(ax, x, y + rng.normal(scale=0.15, size=y.size), groups, label, fig)

handles, labels = axs.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, title="Condition", loc="outside lower center", ncols=2)
fig.savefig(OUT / "multipanel_2x2.pdf")
fig.savefig(OUT / "multipanel_2x2.png", dpi=DPI)
plt.close(fig)

# A named, spanning layout. `subplot_mosaic` is preferable to nested indexing
# when the intended geometry is easier to read as a matrix.
fig, axd = plt.subplot_mosaic(
    [["A", "A", "B"], ["A", "A", "C"]],
    figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), dpi=DPI,
    layout="constrained", sharex=True, sharey=True,
)
for label, ax in axd.items():
    add_identical_mapping(ax, x, y + rng.normal(scale=0.15, size=y.size), groups, label, fig)
handles, labels = axd["A"].get_legend_handles_labels()
fig.legend(handles, labels, title="Condition", loc="outside lower center", ncols=2)
fig.savefig(OUT / "multipanel_mosaic.pdf")
fig.savefig(OUT / "multipanel_mosaic.png", dpi=DPI)
plt.close(fig)

for path in (OUT / "multipanel_2x2.png", OUT / "multipanel_mosaic.png"):
    with Image.open(path) as image:
        width_px, height_px = image.size
    print(
        f"{path.name}: {width_px} x {height_px} px; "
        f"{width_px / DPI * 25.4:.2f} x {height_px / DPI * 25.4:.2f} mm"
    )
print("PDFs use pdf.fonttype=42; confirm page size/fonts with pdfinfo and pdffonts.")
