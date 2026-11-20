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
# Script:  04_indices_per_unit.R
# Purpose: Richness (S), PD and FRic standardized by the regional pool
#          (S / 119; PD / PD of the pool; FRic is already relative to the
#          volume of the pool) and compared among indices within each
#          sampling category (Friedman + paired Wilcoxon).
# Input:   phylogeny/results and functional/raw/results (from 01-03)
# Output:  integrated/{results,figures}
# =============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

out_res <- "integrated/results"
out_fig <- "integrated/figures"
for (d in c(out_res, out_fig)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------
pd   <- read.csv("phylogeny/results/communities_alpha_pd_ses.csv")
fric <- read.csv("functional/raw/results/communities_fric_ses.csv")
pd_pool <- read.csv("phylogeny/results/pd_pool_summary.csv")$PD
S_pool <- 119

indices_unit <- pd %>%
  select(Method, Locality, S, PD) %>%
  inner_join(fric %>% select(Method, Locality, FRic), by = c("Method", "Locality")) %>%
  mutate(S_rel = S / S_pool, PD_rel = PD / pd_pool, FRic_rel = FRic)

methods <- c("Traditional", "eDNA", "Combined")
indices <- c(S_rel = "Richness", PD_rel = "PD", FRic_rel = "FRic")

long <- indices_unit %>%
  select(Method, Locality, all_of(names(indices))) %>%
  pivot_longer(all_of(names(indices)), names_to = "Index", values_to = "value") %>%
  # only samples with all three indices (FRic requires S > 5)
  group_by(Method, Locality) %>% filter(sum(!is.na(value)) == 3) %>% ungroup() %>%
  mutate(Method = factor(Method, levels = methods), Index = factor(Index, levels = names(indices)))

# -----------------------------------------------------------------------------
# 2. Descriptive statistics
# -----------------------------------------------------------------------------
descriptive <- long %>%
  group_by(Method, Index) %>%
  summarise(n = n(), Mean = mean(value), SD = sd(value), Median = median(value),
            Q1 = quantile(value, 0.25), Q3 = quantile(value, 0.75), .groups = "drop")

# -----------------------------------------------------------------------------
# 3. Friedman per category + paired Wilcoxon (Holm)
# -----------------------------------------------------------------------------
friedman <- bind_rows(lapply(methods, function(mt) {
  mm <- long %>% filter(Method == mt) %>%
    pivot_wider(names_from = Index, values_from = value) %>%
    select(all_of(names(indices))) %>% as.matrix()
  ft <- friedman.test(mm)
  data.frame(Method = mt, n_samples = nrow(mm), Chi2 = unname(ft$statistic), df = unname(ft$parameter),
             p = ft$p.value, W_Kendall = unname(ft$statistic) / (nrow(mm) * 2))
}))
friedman$p_holm <- p.adjust(friedman$p, "holm")

pairs <- list(c("S_rel", "PD_rel"), c("S_rel", "FRic_rel"), c("PD_rel", "FRic_rel"))
posthoc <- bind_rows(lapply(methods, function(mt) {
  w <- long %>% filter(Method == mt) %>% pivot_wider(names_from = Index, values_from = value)
  res <- bind_rows(lapply(pairs, function(pr) {
    x <- w[[pr[1]]]; y <- w[[pr[2]]]
    tt <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE))
    data.frame(Method = mt, Index1 = indices[pr[1]], Index2 = indices[pr[2]], n = length(x),
               Median_diff = median(x - y), V = unname(tt$statistic), p = tt$p.value)
  }))
  res$p_holm <- p.adjust(res$p, "holm")
  res
}))

# letters (same letter = no difference); only if the Friedman test is significant
letters_tab <- bind_rows(lapply(methods, function(mt) {
  ph <- posthoc %>% filter(Method == mt)
  differs <- function(a, b) {
    r0 <- ph %>% filter((Index1 == a & Index2 == b) | (Index1 == b & Index2 == a))
    r0$p_holm < 0.05
  }
  ord <- descriptive %>% filter(Method == mt) %>% arrange(desc(Median)) %>% pull(Index) %>% as.character()
  ord <- unname(indices[ord])
  groups <- list()
  for (g in ord) {
    fits <- vapply(groups, function(G) !any(vapply(G, differs, logical(1), b = g)), logical(1))
    if (any(fits)) {
      for (i in which(fits)) groups[[i]] <- c(groups[[i]], g)
    } else {
      prev <- ord[seq_len(match(g, ord) - 1)]
      groups[[length(groups) + 1]] <- c(prev[!vapply(prev, differs, logical(1), b = g)], g)
    }
  }
  if (friedman$p_holm[friedman$Method == mt] >= 0.05) groups <- list(ord)
  data.frame(Method = mt, Index = ord, row.names = NULL,
             Letter = sapply(ord, function(g) paste(letters[which(sapply(groups, function(G) g %in% G))], collapse = "")))
}))

# -----------------------------------------------------------------------------
# 4. Export
# -----------------------------------------------------------------------------
write.csv(indices_unit, file.path(out_res, "indices_per_unit.csv"), row.names = FALSE)
write.csv(descriptive, file.path(out_res, "descriptive.csv"), row.names = FALSE)
write.csv(friedman, file.path(out_res, "friedman_by_category.csv"), row.names = FALSE)
write.csv(posthoc, file.path(out_res, "posthoc_wilcoxon.csv"), row.names = FALSE)
write.csv(letters_tab, file.path(out_res, "letters.csv"), row.names = FALSE)

wb <- createWorkbook()
for (nm in c("indices_unit", "descriptive", "friedman", "posthoc", "letters_tab")) {
  addWorksheet(wb, nm)
  writeData(wb, nm, get(nm))
}
saveWorkbook(wb, file.path(out_res, "comparison_S_PD_FRic.xlsx"), overwrite = TRUE)

# -----------------------------------------------------------------------------
# 5. Figure: violins of the 3 indices side by side, one panel per category
# -----------------------------------------------------------------------------
index_cols <- c(Richness = "#8c6d46", PD = "#3e9e5e", FRic = "#6b4c9a")

d <- long %>% mutate(Index = factor(indices[as.character(Index)], levels = indices),
                     Panel = factor(as.character(Method), levels = methods))
means <- d %>% group_by(Panel, Index) %>%
  summarise(m = mean(value), n = n(), se = sd(value) / sqrt(n()), .groups = "drop") %>%
  mutate(lo = m - qt(0.975, n - 1) * se, hi = m + qt(0.975, n - 1) * se)
fr_lab <- friedman %>% mutate(lab = sprintf("%s\nFriedman χ² = %.1f, %s; W = %.2f (n = %d)", Method, Chi2,
                ifelse(p_holm < 0.001, "P < 0.001", sprintf("P = %.3f", p_holm)), W_Kendall, n_samples))
panel_labels <- setNames(fr_lab$lab, fr_lab$Method)
axis_grey <- "grey20"
title_grey <- "grey20"  # 80% black
# legend key with slightly rounded corners
draw_key_round <- function(data, params, size) {
  grid::roundrectGrob(width = unit(0.85, "npc"), height = unit(0.85, "npc"), r = unit(0.18, "snpc"),
                      gp = grid::gpar(col = if (is.null(data$colour) || is.na(data$colour)) NA else data$colour,
                                      fill = scales::alpha(data$fill, data$alpha), lwd = 0.8))
}

# bracket under the three indices showing that the category spans all three
y_bracket <- -0.075; h_bracket <- 0.012
bracket <- data.frame(Panel = factor(levels(d$Panel), levels = levels(d$Panel)),
                      x0 = 0.75, x1 = 3.25, y = y_bracket, yt = y_bracket + h_bracket)

# category and statistics below each panel; indices along the x axis
fig <- ggplot(d, aes(Index, value)) +
  geom_violin(aes(fill = Index, colour = Index), alpha = 0.25, linewidth = 0.4, trim = TRUE, scale = "width", width = 0.8,
              key_glyph = draw_key_round) +
  geom_line(aes(group = Locality), colour = "grey80", linewidth = 0.25) +
  geom_jitter(aes(colour = Index), width = 0.08, size = 1, alpha = 0.55, show.legend = FALSE) +
  geom_segment(data = bracket, aes(x = x0, xend = x1, y = y, yend = y), inherit.aes = FALSE,
               colour = axis_grey, linewidth = 0.4) +
  geom_segment(data = bracket, aes(x = x0, xend = x0, y = y, yend = yt),
               inherit.aes = FALSE, colour = axis_grey, linewidth = 0.4) +
  geom_segment(data = bracket, aes(x = x1, xend = x1, y = y, yend = yt),
               inherit.aes = FALSE, colour = axis_grey, linewidth = 0.4) +
  geom_line(data = means, aes(Index, m, group = 1), colour = axis_grey, linewidth = 0.5) +
  geom_errorbar(data = means, aes(Index, ymin = lo, ymax = hi), inherit.aes = FALSE, width = 0.12,
                linewidth = 0.7, colour = axis_grey) +
  geom_point(data = means, aes(Index, m, fill = Index), shape = 21, size = 3, stroke = 0.8,
             colour = axis_grey, show.legend = FALSE) +
  facet_wrap(~Panel, nrow = 1, strip.position = "bottom", labeller = as_labeller(panel_labels)) +
  coord_cartesian(ylim = c(0, NA), clip = "off") +
  scale_fill_manual(values = index_cols, name = NULL) +
  scale_colour_manual(values = index_cols, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.02, 0.05))) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.25))) +
  labs(x = NULL, y = "Proportion of the regional pool (119 taxa)",
       title = "Richness, PD and FRic Recovered from the Regional Pool by Sampling Method") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.line = element_line(colour = axis_grey, linewidth = 0.4),
        axis.ticks = element_line(colour = axis_grey, linewidth = 0.4),
        axis.text = element_text(colour = axis_grey), axis.title = element_text(colour = axis_grey),
        strip.placement = "outside",
        strip.text = element_text(size = 9.5, colour = axis_grey, lineheight = 1.2, margin = margin(t = 12, b = 4)),
        plot.title = element_text(face = "bold", size = 11.5, colour = title_grey, hjust = 0.5),
        legend.position = "bottom", panel.spacing = unit(1.2, "lines"))

ggsave(file.path(out_fig, "fig_violins_S_PD_FRic.png"), fig, width = 11, height = 5.4, dpi = 300, bg = "white")
ggsave(file.path(out_fig, "fig_violins_S_PD_FRic.pdf"), fig, width = 11, height = 5.4, device = cairo_pdf)

# -----------------------------------------------------------------------------
cat("\n== Descriptive ==\n"); print(as.data.frame(descriptive), digits = 3)
cat("\n== Friedman ==\n"); print(friedman, digits = 3)
cat("\n== Post hoc ==\n"); print(posthoc, digits = 3)
cat("\n== Letters ==\n"); print(letters_tab)
