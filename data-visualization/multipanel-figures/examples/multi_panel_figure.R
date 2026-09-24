# Create a 183 x 140 mm multi-panel figure with one qualified shared guide.
# Run: Rscript examples/multi_panel_figure.R
# The four panels intentionally use identical x/y limits and the same colour scale,
# so collecting axes and the Condition guide is semantically valid.

library(ggplot2)
library(patchwork)

width_mm <- 183
height_mm <- 140
dpi <- 300
set.seed(42)

df <- data.frame(
  x = rnorm(240),
  y = rnorm(240),
  group = factor(rep(c("Control", "Treatment"), each = 120))
)
df$y <- df$y + ifelse(df$group == "Treatment", 0.35, -0.35)

theme_publication <- theme_classic(base_size = 7) +
  theme(
    plot.title = element_text(size = 7, face = "plain"),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25)
  )
palette <- c(Control = "#4DBBD5", Treatment = "#E64B35")

make_panel <- function(seed, title) {
  set.seed(seed)
  d <- transform(df, y = y + rnorm(nrow(df), sd = 0.18))
  ggplot(d, aes(x, y, colour = group)) +
    geom_point(size = 1.1, alpha = 0.72) +
    scale_colour_manual(values = palette, name = "Condition") +
    coord_cartesian(xlim = c(-3.2, 3.2), ylim = c(-3.2, 3.2)) +
    labs(x = "Standardized x", y = "Standardized y", title = title) +
    theme_publication
}

panels <- Map(make_panel, 1:4, c("Discovery", "Replication", "Validation", "Sensitivity"))
figure <- wrap_plots(panels, ncol = 2) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "a") &
  # `& theme()` is required: the theme in plot_annotation() does not style tags.
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.location = "panel",
    plot.tag.position = c(0.02, 0.98),
    legend.position = "bottom"
  )

save_cairo_pdf <- function(path, plot, width_mm, height_mm) {
  Cairo::CairoPDF(path, width = width_mm / 25.4, height = height_mm / 25.4)
  on.exit(dev.off(), add = TRUE)
  print(plot)
}
save_cairo_pdf("Figure1.pdf", figure, width_mm, height_mm)
ggsave("Figure1.png", figure, width = width_mm, height = height_mm,
       units = "mm", dpi = dpi, device = ragg::agg_png)

# PNG dimensions provide a portable measurement during the example run; use
# `pdfinfo Figure1.pdf` as the final PDF page-size check before submission.
png_dim <- dim(png::readPNG("Figure1.png"))
png_width <- png_dim[2]
png_height <- png_dim[1]
message(sprintf("Figure1.pdf target page: %d x %d mm (verify with pdfinfo)", width_mm, height_mm))
message(sprintf("Figure1.png measured: %d x %d px = %.2f x %.2f mm at %d dpi",
                png_width, png_height, png_width / dpi * 25.4, png_height / dpi * 25.4, dpi))
