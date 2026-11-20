# =============================================================================
# Rethinking stream fish assessments: environmental DNA broadens the
# multidimensional diversity of assemblages in Eastern Amazonia
#
# Authors: Fabricio dos Anjos Santa Rosa, Derlan José Ferreira Silva,
#          Cintia Oliveira Carvalho, Silvia Britto Barreto,
#          Silvana Melo Sviggum, Anna Karolina Oliveira de Queiroz,
#          Enki Gautier, Marcia da Silva Anjos, Jarl Andreas Anmarkrud,
#          João Bráullio de Luna Sales, Quentin Mauvisseau, Hugo de Boer,
#          Jonathan Stuart Ready
# Contact: Fabricio dos Anjos Santa Rosa (fabriciodosanjos14@gmail.com)
#
# Script:  05_complementarity_FigS1_S2_S3.R
# Purpose: Exclusive contribution of each method per stream-campaign unit
#          (eDNA = (combined - traditional) / combined; traditional =
#          (combined - eDNA) / combined), paired Wilcoxon and Friedman tests;
#          Figure S1 (Spearman correlation among the 33 traits), Figure S2
#          (pool tree coloured by detection category) and Figure S3
#          (functional space with the convex hull of each method).
# Input:   phylogeny/, functional/raw/ and integrated/ outputs (from 01-04)
#          silhouettes/silhouettes_polygons.csv (fish outlines for Figures
#          S2-S3; folder set by the environment variable SIL_DIR)
# Output:  integrated/{results,figures}
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(geometry)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(openxlsx)
})

out_res <- "integrated/results"
out_fig <- "integrated/figures"
for (d in c(out_res, out_fig)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# colours per detection category (checked for colour blindness)
detection_cols <- c("Traditional only" = "#d1495b", "eDNA only" = "#2a78b5", "Both methods" = "#b8860b")
index_cols <- c(Richness = "#8c6d46", PD = "#3e9e5e", FRic = "#6b4c9a")

# style shared by all figures (main and supplementary)
axis_grey <- "grey20"
title_grey <- "grey20"  # 80% black
# legend key with slightly rounded corners
draw_key_round <- function(data, params, size) {
  grid::roundrectGrob(width = unit(0.85, "npc"), height = unit(0.85, "npc"), r = unit(0.18, "snpc"),
                      gp = grid::gpar(col = if (is.null(data$colour) || is.na(data$colour)) NA else data$colour,
                                      fill = scales::alpha(data$fill, data$alpha), lwd = 0.8))
}
fig_theme <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          axis.line = element_line(colour = axis_grey, linewidth = 0.4),
          axis.ticks = element_line(colour = axis_grey, linewidth = 0.4),
          axis.text = element_text(colour = axis_grey), axis.title = element_text(colour = axis_grey),
          legend.text = element_text(colour = axis_grey), legend.title = element_text(colour = axis_grey),
          plot.title = element_text(face = "bold", size = base_size + 1.5, colour = title_grey, hjust = 0.5),
          plot.caption = element_text(colour = axis_grey, hjust = 0.5),
          legend.position = "bottom")
}

# -----------------------------------------------------------------------------
# 1. Detection category of each taxon in the pool
# -----------------------------------------------------------------------------
mat <- read.csv("phylogeny/results/community_matrix_pa.csv", row.names = 1, check.names = FALSE)
tip_names <- read.csv("phylogeny/results/taxon_tip_match.csv")
row_method <- sub("\\|.*", "", rownames(mat))
det_trad <- colSums(mat[row_method == "Traditional", , drop = FALSE]) > 0
det_edna <- colSums(mat[row_method == "eDNA", , drop = FALSE]) > 0

detection <- data.frame(Tree_tip = colnames(mat), Trad = det_trad, eDNA = det_edna) %>%
  mutate(Category = case_when(Trad & eDNA ~ "Both methods", Trad ~ "Traditional only", TRUE ~ "eDNA only")) %>%
  left_join(tip_names %>% select(Taxon, Tree_tip, Substitute), by = "Tree_tip") %>%
  mutate(Category = factor(Category, levels = names(detection_cols)))
stopifnot(nrow(detection) == 119, !anyNA(detection$Taxon))
write.csv(detection, file.path(out_res, "detection_by_taxon.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. Figure 3 - gain per sampling unit
# -----------------------------------------------------------------------------
idx <- read.csv("integrated/results/indices_per_unit.csv")
wide <- idx %>%
  select(Method, Locality, S, PD, FRic) %>%
  pivot_longer(c(S, PD, FRic), names_to = "Index", values_to = "value") %>%
  pivot_wider(names_from = Method, values_from = value) %>%
  filter(!is.na(Traditional), !is.na(eDNA), !is.na(Combined), Combined > 0) %>%
  mutate(gain_eDNA = (Combined - Traditional) / Combined,
         gain_Trad = (Combined - eDNA) / Combined,
         Index = factor(recode(Index, S = "Richness"), levels = names(index_cols)))

tests <- wide %>% group_by(Index) %>%
  summarise(n = n(),
            mean_gain_eDNA = mean(gain_eDNA), mean_gain_Trad = mean(gain_Trad),
            median_gain_eDNA = median(gain_eDNA), median_gain_Trad = median(gain_Trad),
            V = unname(suppressWarnings(wilcox.test(gain_eDNA, gain_Trad, paired = TRUE, exact = FALSE))$statistic),
            p = suppressWarnings(wilcox.test(gain_eDNA, gain_Trad, paired = TRUE, exact = FALSE))$p.value,
            .groups = "drop") %>%
  mutate(p_holm = p.adjust(p, "holm"))

# Friedman test of the eDNA gain among indices (same units in all three)
m_fr <- wide %>% select(Locality, Index, gain_eDNA) %>%
  pivot_wider(names_from = Index, values_from = gain_eDNA) %>% na.omit()
fr <- friedman.test(as.matrix(m_fr[, names(index_cols)]))
friedman_gain <- data.frame(n = nrow(m_fr), Chi2 = unname(fr$statistic), df = unname(fr$parameter),
                            p = fr$p.value, W_Kendall = unname(fr$statistic) / (nrow(m_fr) * 2))
pairs <- combn(names(index_cols), 2, simplify = FALSE)
posthoc_gain <- bind_rows(lapply(pairs, function(pr) {
  tt <- suppressWarnings(wilcox.test(m_fr[[pr[1]]], m_fr[[pr[2]]], paired = TRUE, exact = FALSE))
  data.frame(Index1 = pr[1], Index2 = pr[2], n = nrow(m_fr),
             Median_diff = median(m_fr[[pr[1]]] - m_fr[[pr[2]]]), V = unname(tt$statistic), p = tt$p.value)
})) %>% mutate(p_holm = p.adjust(p, "holm"))

write.csv(wide, file.path(out_res, "gain_per_unit.csv"), row.names = FALSE)
write.csv(tests, file.path(out_res, "gain_wilcoxon_eDNA_vs_Trad.csv"), row.names = FALSE)
write.csv(friedman_gain, file.path(out_res, "gain_friedman_indices.csv"), row.names = FALSE)
write.csv(posthoc_gain, file.path(out_res, "gain_posthoc_indices.csv"), row.names = FALSE)

gain_long <- wide %>%
  select(Locality, Index, gain_eDNA, gain_Trad) %>%
  pivot_longer(c(gain_eDNA, gain_Trad), names_to = "Source", values_to = "gain") %>%
  mutate(Source = factor(recode(Source, gain_eDNA = "Added by eDNA", gain_Trad = "Added by traditional"),
                         levels = c("Added by traditional", "Added by eDNA")))
source_cols <- c("Added by traditional" = "#d1495b", "Added by eDNA" = "#2a78b5")

gain_summary <- gain_long %>% group_by(Index, Source) %>%
  summarise(m = mean(gain), ci = qt(0.975, n() - 1) * sd(gain) / sqrt(n()), .groups = "drop")
# index + statistics as x-axis label
x_labels <- tests %>%
  mutate(lab = paste0(Index, "\nWilcoxon ", ifelse(p_holm < 0.001, "P < 0.001", sprintf("P = %.3f", p_holm)),
                      " (n = ", n, ")"))
x_labels <- setNames(x_labels$lab, x_labels$Index)

dodge <- position_dodge(width = 0.8)
fig3 <- ggplot(gain_long, aes(Index, gain, fill = Source, colour = Source)) +
  geom_violin(position = dodge, alpha = 0.25, trim = TRUE, linewidth = 0.4, scale = "width", width = 0.75,
              key_glyph = draw_key_round) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8, seed = 1),
             size = 1, alpha = 0.55, show.legend = FALSE) +
  geom_errorbar(data = gain_summary, aes(Index, ymin = m - ci, ymax = m + ci, group = Source), inherit.aes = FALSE,
                position = dodge, width = 0.12, linewidth = 0.6, colour = axis_grey) +
  geom_point(data = gain_summary, aes(Index, m, group = Source, fill = Source), inherit.aes = FALSE,
             position = dodge, shape = 21, size = 3, colour = axis_grey, show.legend = FALSE) +
  scale_fill_manual(values = source_cols, name = NULL) +
  scale_colour_manual(values = source_cols, name = NULL) +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = function(x) paste0(x * 100, "%"),
                     expand = expansion(mult = c(0.01, 0.03))) +
  labs(x = NULL, y = "Share of the combined diversity\ncontributed exclusively by each method",
       title = "Exclusive Contribution of Each Sampling Method to the Combined Diversity") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.line = element_line(colour = axis_grey, linewidth = 0.4),
        axis.ticks = element_line(colour = axis_grey, linewidth = 0.4),
        axis.text = element_text(colour = axis_grey), axis.title = element_text(colour = axis_grey),
        axis.text.x = element_text(size = 10, lineheight = 1.2),
        plot.title = element_text(face = "bold", size = 11.5, colour = title_grey, hjust = 0.5),
        legend.position = "bottom")
ggsave(file.path(out_fig, "fig3_gain_per_unit.png"), fig3, width = 7.5, height = 5.4, dpi = 300, bg = "white")
ggsave(file.path(out_fig, "fig3_gain_per_unit.pdf"), fig3, width = 7.5, height = 5.4, device = cairo_pdf)

# -----------------------------------------------------------------------------
# 3. Figure S2 - pool tree coloured by detection
# -----------------------------------------------------------------------------
tree <- read.tree("phylogeny/data/pool_tree_119taxa.nwk")
tip_category <- detection$Category[match(tree$tip.label, detection$Tree_tip)]
tip_label <- detection$Taxon[match(tree$tip.label, detection$Tree_tip)]
tip_label <- ifelse(detection$Substitute[match(tree$tip.label, detection$Tree_tip)], paste0(tip_label, "*"), tip_label)
stopifnot(!anyNA(tip_category))

# terminal branches coloured by category; internal branches in grey
branch_col <- rep("grey55", nrow(tree$edge))
terminal <- tree$edge[, 2] <= Ntip(tree)
branch_col[terminal] <- detection_cols[as.character(tip_category[tree$edge[terminal, 2]])]

legend_labels <- paste0(names(detection_cols), " (", as.vector(table(tip_category)[names(detection_cols)]), " taxa)")

# labels with the name used in the study (not the tip name); substitutes marked with an asterisk
plot_tree <- tree; plot_tree$tip.label <- gsub("_", " ", tip_label)

# silhouette at the tip of each taxon (silhouettes/silhouettes_polygons.csv, vector outlines generated
# from the PNGs in the same folder), coloured by detection category, in a ring beyond
# the longest label, head towards the centre, linked to the name by a dashed line of the same colour
sil_dir <- Sys.getenv("SIL_DIR", "silhouettes")
sil_polygons <- read.csv(file.path(sil_dir, "silhouettes_polygons.csv"))
tip_taxon <- detection$Taxon[match(tree$tip.label, detection$Tree_tip)]
draw_silhouettes <- function(r0, max_w, max_h, r_label) {
  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  for (i in seq_len(Ntip(plot_tree))) {
    pol <- sil_polygons[sil_polygons$Taxon == tip_taxon[i], ]
    if (nrow(pol) == 0) next
    fill_col <- detection_cols[as.character(tip_category[i])]
    x0 <- pol$body_x0[1]                       # start of the body: barbels ahead of the head stay inside the ring
    w <- max(pol$x) - x0; h <- max(pol$y); s <- min(max_w / w, max_h / h)
    th <- atan2(lp$yy[i], lp$xx[i]); u <- c(cos(th), sin(th)); nr <- c(-u[2], u[1])
    # outline in the radial frame: axis 1 along the radius (head inwards), axis 2 perpendicular
    ra <- r0 + (pol$x - x0) * s
    pe <- (pol$y - h / 2) * s
    if (u[1] < 0) pe <- -pe                       # left half: mirror so the silhouette is not upside down
    x <- ra * u[1] + pe * nr[1]; y <- ra * u[2] + pe * nr[2]
    rings <- split(seq_along(x), pol$ring)
    xs <- unlist(lapply(rings, function(k) c(x[k], NA))); ys <- unlist(lapply(rings, function(k) c(y[k], NA)))
    polypath(head(xs, -1), head(ys, -1), col = fill_col, border = NA, rule = "evenodd")
    # dashed line from the end of the label to the silhouette
    if (r0 - r_label[i] > 0.02 * r0)
      segments((r_label[i] + 0.012 * r0) * u[1], (r_label[i] + 0.012 * r0) * u[2],
               (r0 - 0.012 * r0) * u[1], (r0 - 0.012 * r0) * u[2], col = fill_col, lty = "22", lwd = 0.6)
  }
}

# rounded rectangle (base graphics) for the legend keys
rounded_rect <- function(xc, yc, w, h, r = 0.25, col, border) {
  rx <- r * min(w, h) * (diff(par("usr")[1:2]) / diff(par("usr")[3:4])) * (par("pin")[2] / par("pin")[1])
  ry <- r * min(w, h)
  t <- seq(0, pi / 2, length.out = 8)
  cx <- c(xc + w / 2 - rx, xc - w / 2 + rx, xc - w / 2 + rx, xc + w / 2 - rx)
  cy <- c(yc + h / 2 - ry, yc + h / 2 - ry, yc - h / 2 + ry, yc - h / 2 + ry)
  ang <- list(t, t + pi / 2, t + pi, t + 3 * pi / 2)
  x <- unlist(lapply(1:4, function(i) cx[i] + rx * cos(ang[[i]])))
  y <- unlist(lapply(1:4, function(i) cy[i] + ry * sin(ang[[i]])))
  polygon(x, y, col = col, border = border, lwd = 1.2)
}

draw_s2 <- function() {
  par(mar = c(4.2, 0.5, 2.2, 0.5), xpd = NA, family = "sans")
  h <- max(node.depth.edgelength(plot_tree))
  plot(plot_tree, type = "fan", tip.color = detection_cols[as.character(tip_category)], font = 3,
       edge.color = branch_col, edge.width = 1.1, cex = 0.42, label.offset = h * 0.01,
       no.margin = FALSE, x.lim = c(-1.62, 1.62) * h, y.lim = c(-1.62, 1.62) * h)
  r_label <- h * 1.01 + strwidth(plot_tree$tip.label, cex = 0.42, font = 3)   # radius at the end of each label
  draw_silhouettes(max(r_label) + h * 0.03, max_w = 0.17 * h, max_h = 0.055 * h, r_label = r_label)
  title(main = "Phylogenetic Distribution of the Taxa Detected by Each Sampling Method",
        col.main = title_grey, font.main = 2, cex.main = 1.15, line = 0.6)
  # legend at the bottom with rounded keys
  lg <- legend("bottom", legend = legend_labels, horiz = TRUE, bty = "n", text.col = axis_grey, cex = 0.8,
               fill = NA, border = NA, x.intersp = 0.6, text.width = strwidth(legend_labels, cex = 0.8) * 1.25,
               inset = c(0, -0.035))
  w <- strheight("M", cex = 0.8) * 1.5 * (par("pin")[2] / par("pin")[1]) * diff(par("usr")[1:2]) / diff(par("usr")[3:4])
  hh <- strheight("M", cex = 0.8) * 1.5
  for (i in seq_along(legend_labels))
    rounded_rect(lg$text$x[i] - w * 0.85, lg$text$y[i], w, hh,
                 col = adjustcolor(detection_cols[i], 0.25), border = detection_cols[i])
}
png(file.path(out_fig, "figS2_tree_detection.png"), width = 9, height = 9.4, units = "in", res = 300, bg = "white")
draw_s2(); invisible(dev.off())
cairo_pdf(file.path(out_fig, "figS2_tree_detection.pdf"), width = 9, height = 9.4)
draw_s2(); invisible(dev.off())

# -----------------------------------------------------------------------------
# 4. Figure S3 - functional space by detection category
# -----------------------------------------------------------------------------
pc <- read.csv("functional/raw/data/functional_space_pcoa.csv", row.names = 1)
quality <- read.csv("functional/raw/results/functional_space_quality.csv")
pts <- pc %>% mutate(Taxon = rownames(pc)) %>%
  left_join(detection %>% select(Taxon, Category), by = "Taxon")
stopifnot(!anyNA(pts$Category))

# % variance per axis: variance of the coordinates, scaled to the cumulative variance of the 5 axes
v <- sapply(pc, var)
pct_var <- 100 * v / sum(v) * quality$Cumulative_variance[quality$Axes == ncol(pc)]
axis_labels <- setNames(sprintf("%s (%.1f%%)", names(pc), pct_var), names(pc))

hull <- function(d, axes, group) {
  h <- chull(d[[axes[1]]], d[[axes[2]]])
  data.frame(x = d[[axes[1]]][h], y = d[[axes[2]]][h], Group = group)
}
# silhouettes in the functional space: each taxon near its point, pushed apart where they overlap
# and linked to the actual point by a thin line in the taxon colour when displaced (panel coordinates in 0-1)
place_silhouettes <- function(px, py, width, height, sx = px, sy = py) {
  # spiral search around each point (all directions) for the nearest free position;
  # among free positions at the same distance, prefer those not covering another taxon's point;
  # sx, sy: centre of the search (the point itself, except for silhouettes pulled inside the plot)
  n <- length(px); x <- px; y <- py; placed <- rep(FALSE, n)
  n_near <- rowSums(outer(px, px, function(a, b) abs(a - b) < 0.06) & outer(py, py, function(a, b) abs(a - b) < 0.06))
  ang <- seq(0, 2 * pi, length.out = 33)[-33]
  for (i in order(n_near)) {                    # isolated taxa first: they stay on their own point
    for (r in c(0, seq(0.006, 0.6, by = 0.006))) {
      cx <- sx[i] + r * cos(ang); cy <- sy[i] + r * sin(ang)
      if (r == 0) { cx <- sx[i]; cy <- sy[i] }
      ok <- cx >= width[i] / 2 & cx <= 1 - width[i] / 2 & cy >= height[i] / 2 & cy <= 1 - height[i] / 2
      j <- which(placed)
      if (length(j)) ok <- ok & sapply(seq_along(cx), function(k)
        !any(abs(cx[k] - x[j]) < (width[i] + width[j]) / 2 & abs(cy[k] - y[j]) < (height[i] + height[j]) / 2))
      if (any(ok)) {
        others <- setdiff(seq_len(n), i)
        covers <- sapply(seq_along(cx), function(k)
          sum(abs(px[others] - cx[k]) < width[i] / 2 & abs(py[others] - cy[k]) < height[i] / 2))
        k <- which(ok)[which.min(covers[ok])]
        x[i] <- cx[k]; y[i] <- cy[k]; placed[i] <- TRUE; break
      }
    }
  }
  data.frame(x = x, y = y)
}
sil_dims <- sil_polygons %>% group_by(Taxon) %>% summarise(w = max(x), h = max(y), .groups = "drop")

hull_cols <- c("Pool (119 taxa)" = "grey35", Traditional = "#d1495b", eDNA = "#2a78b5")
# single legend: the hull of each method uses the symbol of the taxa detected only by that method
key_map <- c("Traditional only" = "Traditional", "eDNA only" = "eDNA", "Both methods" = "Both methods")
key_cols <- c(Traditional = "#d1495b", eDNA = "#2a78b5", "Both methods" = "#b8860b", "Pool (119 taxa)" = "grey35")
key_shapes <- c(Traditional = 17, eDNA = 15, "Both methods" = 16, "Pool (119 taxa)" = 32)  # 32: empty symbol
panel_s3 <- function(axes) {
  e1 <- pts[[axes[1]]]; e2 <- pts[[axes[2]]]
  env <- bind_rows(hull(pts, axes, "Pool (119 taxa)"),
                   hull(pts %>% filter(Category != "eDNA only"), axes, "Traditional"),
                   hull(pts %>% filter(Category != "Traditional only"), axes, "eDNA")) %>%
    mutate(Group = factor(Group, levels = names(hull_cols)))
  xr <- range(e1) + c(-1, 1) * 0.09 * diff(range(e1)); yr <- range(e2) + c(-1, 1) * 0.09 * diff(range(e2))
  nx <- (e1 - xr[1]) / diff(xr); ny <- (e2 - yr[1]) / diff(yr)
  d <- sil_dims[match(pts$Taxon, sil_dims$Taxon), ]
  s <- pmin(0.066 / d$w, 0.03 / d$h); lw <- d$w * s; lh <- d$h * s
  # silhouettes of taxa at hull vertices, pulled inwards so they are not confused with the lines
  pull_in <- pts$Taxon %in% c("Potamotrygon orbignyi", "Synbranchus marmoratus")
  dc <- sqrt((0.5 - nx)^2 + (0.5 - ny)^2)
  sx <- ifelse(pull_in, nx + 0.13 * (0.5 - nx) / dc, nx); sy <- ifelse(pull_in, ny + 0.13 * (0.5 - ny) / dc, ny)
  pos <- place_silhouettes(nx, ny, lw * 1.06, lh * 1.12, sx, sy)
  sil <- bind_rows(lapply(seq_len(nrow(pts)), function(i) {
    p <- sil_polygons[sil_polygons$Taxon == pts$Taxon[i], ]
    data.frame(Taxon = pts$Taxon[i], K = key_map[as.character(pts$Category[i])], ring = p$ring,
               X = xr[1] + (pos$x[i] + (p$x - d$w[i] / 2) * s[i]) * diff(xr),
               Y = yr[1] + (pos$y[i] + (p$y - d$h[i] / 2) * s[i]) * diff(yr))
  }))
  # line from the point to the nearest edge of the silhouette, when it does not cover the point
  bx <- pmin(pmax(nx, pos$x - lw / 2), pos$x + lw / 2); by <- pmin(pmax(ny, pos$y - lh / 2), pos$y + lh / 2)
  links <- data.frame(x = e1, y = e2, xend = xr[1] + bx * diff(xr), yend = yr[1] + by * diff(yr),
                      K = key_map[as.character(pts$Category)])[sqrt((bx - nx)^2 + (by - ny)^2) > 0.006, ]
  ggplot() +
    geom_polygon(data = env, aes(x, y, colour = Group, fill = Group), alpha = 0.06, linewidth = 0.7,
                 key_glyph = draw_key_round) +
    geom_segment(data = links, aes(x, y, xend = xend, yend = yend, colour = K), linewidth = 0.12,
                 show.legend = FALSE) +
    geom_polygon(data = sil, aes(X, Y, group = Taxon, subgroup = ring, fill = K), colour = NA,
                 alpha = 0.95, show.legend = FALSE) +
    geom_point(data = pts %>% mutate(K = key_map[as.character(Category)]),
               aes(.data[[axes[1]]], .data[[axes[2]]], shape = K, colour = K), size = 1.3, alpha = 0.9) +
    scale_colour_manual(values = key_cols, breaks = names(key_cols), limits = names(key_cols), name = NULL) +
    scale_fill_manual(values = key_cols, breaks = names(key_cols), limits = names(key_cols), name = NULL) +
    scale_shape_manual(values = key_shapes, breaks = names(key_cols), limits = names(key_cols), name = NULL) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2.5)),
           fill = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
    coord_cartesian(xlim = xr, ylim = yr, expand = FALSE) +
    labs(x = axis_labels[axes[1]], y = axis_labels[axes[2]]) +
    fig_theme(10) +
    theme(aspect.ratio = 1)
}
fig_s3 <- (panel_s3(c("PC1", "PC2")) + panel_s3(c("PC3", "PC4"))) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",
                  title = "Functional Space Occupied by the Taxa Detected by Each Sampling Method",
                  theme = fig_theme(10)) &
  theme(legend.position = "bottom", legend.box = "vertical", plot.tag = element_text(colour = title_grey, face = "bold"))
ggsave(file.path(out_fig, "figS3_functional_space.png"), fig_s3, width = 14, height = 7.6, dpi = 300, bg = "white")
ggsave(file.path(out_fig, "figS3_functional_space.pdf"), fig_s3, width = 14, height = 7.6, device = cairo_pdf)

# -----------------------------------------------------------------------------
# 4b. Figure S1 - Spearman correlation among the 33 traits
# -----------------------------------------------------------------------------
tr_all <- read.csv("functional/raw/data/traits_standardized_all.csv", row.names = 1, check.names = FALSE)
trait_defs <- read.csv("functional/raw/results/trait_definitions_selection.csv")
trait_labels <- c(
  body_mass = "Body mass", mouth_gape = "Mouth gape", snout_length = "Snout length",
  gut_length = "Gut length", longest_barbel_length = "Longest barbel length",
  mouth_width = "Mouth width", body_width = "Body width",
  caudal_peduncle_width = "Caudal peduncle width", head_width = "Head width",
  mouth_protrusion = "Mouth protrusion", head_depth = "Head depth",
  anal_fin_height = "Anal fin height", caudal_fin_height = "Caudal fin height",
  dorsal_fin_height = "Dorsal fin height", pectoral_fin_height = "Pectoral fin height",
  max_body_depth = "Maximum body depth",
  body_depth_at_pectoral_insertion = "Body depth at pectoral fin insertion",
  caudal_peduncle_depth = "Caudal peduncle depth", caudal_fin_area = "Caudal fin area",
  pectoral_fin_area = "Pectoral fin area", body_area = "Body area",
  head_length = "Head length", anal_fin_length = "Anal fin length",
  caudal_fin_length = "Caudal fin length", dorsal_fin_length = "Dorsal fin length",
  pectoral_fin_length = "Pectoral fin length", pelvic_fin_length = "Pelvic fin length",
  standard_length = "Standard length", eye_area = "Eye area", eye_diameter = "Eye diameter",
  mouth_height = "Mouth height", eye_position = "Eye position", pectoral_fin_position = "Pectoral fin position")
stopifnot(all(colnames(tr_all) %in% names(trait_labels)))
retained <- trait_defs$Trait[trait_defs$Retained]

rho_all <- cor(tr_all, method = "spearman")
trait_order <- colnames(rho_all)[hclust(as.dist(1 - abs(rho_all)))$order]
d_s1 <- as.data.frame(as.table(rho_all)) %>%
  mutate(Var1 = factor(Var1, levels = trait_order), Var2 = factor(Var2, levels = trait_order))
s1_labels <- setNames(ifelse(trait_order %in% retained, paste0(trait_labels[trait_order], " *"),
                             trait_labels[trait_order]), trait_order)
s1_face <- ifelse(trait_order %in% retained, "bold", "plain")
mark <- d_s1 %>% filter(abs(Freq) >= 0.7, Var1 != Var2) %>% mutate(Mark = "|ρ| ≥ 0.7")

fig_s1 <- ggplot(d_s1, aes(Var1, Var2, fill = Freq)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_point(data = mark, aes(shape = Mark), size = 0.9, colour = axis_grey) +
  scale_fill_gradient2(low = "#2a78b5", mid = "white", high = "#d1495b", limits = c(-1, 1),
                       name = "Spearman's ρ",
                       guide = guide_colourbar(barwidth = unit(5, "cm"), barheight = unit(0.35, "cm"),
                                               title.vjust = 0.8, order = 1)) +
  scale_shape_manual(values = c("|ρ| ≥ 0.7" = 4), name = NULL, guide = guide_legend(order = 2)) +
  scale_x_discrete(labels = s1_labels) + scale_y_discrete(labels = s1_labels) +
  labs(x = NULL, y = NULL, title = "Correlation among the 33 Morphological Traits",
       caption = sprintf("* trait retained after the correlation filter (%d of %d)", length(retained), ncol(tr_all))) +
  coord_fixed() +
  fig_theme(9) +
  theme(panel.grid = element_blank(),
        axis.text.x = suppressWarnings(element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, face = s1_face)),
        axis.text.y = suppressWarnings(element_text(size = 7, face = s1_face)))
suppressWarnings({
  ggsave(file.path(out_fig, "figS1_trait_correlation.png"), fig_s1, width = 8.5, height = 9.2, dpi = 300, bg = "white")
  ggsave(file.path(out_fig, "figS1_trait_correlation.pdf"), fig_s1, width = 8.5, height = 9.2, device = cairo_pdf)
})

# -----------------------------------------------------------------------------
# 5. Summary table
# -----------------------------------------------------------------------------
wb <- createWorkbook()
for (nm in c("detection", "tests", "friedman_gain", "posthoc_gain")) {
  addWorksheet(wb, nm); writeData(wb, nm, get(nm))
}
saveWorkbook(wb, file.path(out_res, "complementarity_summary.xlsx"), overwrite = TRUE)

print(table(detection$Category)); print(tests); print(friedman_gain); print(posthoc_gain)
