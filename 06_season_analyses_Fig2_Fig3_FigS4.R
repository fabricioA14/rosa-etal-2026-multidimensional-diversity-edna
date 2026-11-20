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
# Script:  06_season_analyses_Fig2_Fig3_FigS4.R
# Purpose: Analyses within and between hydrological seasons (dry: August
#          2022 and 2023; rainy: February 2023 and 2024). Within each season
#          (stream-campaign units as blocks): Friedman + paired Wilcoxon among
#          sampling categories and among indices, and exclusive contribution
#          of each method. Between seasons (mean of the two campaigns per
#          stream, n = 12): paired Wilcoxon for indices, exclusive
#          contributions and Jaccard similarity between methods; regional
#          (gamma) diversity of each method. P-values are not adjusted.
#          Figures 2, 3 and S4; Table S3.
# Input:   SuppInfo_S1.xlsx, sheet 3_analysis_input; integrated/, phylogeny/
#          and functional/raw/ outputs (from 01-05)
# Output:  integrated/{results,figures}
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(picante)
  library(geometry)
  library(readxl)
  library(openxlsx)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

occ_file <- "SuppInfo_S1.xlsx"
out_res  <- "integrated/results"
out_fig  <- "integrated/figures"
for (d in c(out_res, out_fig)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# style shared with the other figures
axis_grey <- "grey20"
title_grey <- "grey20"  # 80% black
season_cols <- c(Dry = "#c2703d", Rainy = "#3a7d8c")
method_cols <- c(Traditional = "#d1495b", eDNA = "#2a78b5", Combined = "#b8860b")
source_cols <- c("Added by traditional" = "#d1495b", "Added by eDNA" = "#2a78b5")
draw_key_round <- function(data, params, size) {
  grid::roundrectGrob(width = unit(0.85, "npc"), height = unit(0.85, "npc"), r = unit(0.18, "snpc"),
                      gp = grid::gpar(col = if (is.null(data$colour) || is.na(data$colour)) NA else data$colour,
                                      fill = scales::alpha(data$fill, data$alpha), lwd = 0.8))
}
fig_theme <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.line = element_line(colour = axis_grey, linewidth = 0.4),
          axis.ticks = element_line(colour = axis_grey, linewidth = 0.4),
          axis.text = element_text(colour = axis_grey), axis.title = element_text(colour = axis_grey),
          legend.text = element_text(colour = axis_grey), legend.title = element_text(colour = axis_grey),
          strip.text = element_text(colour = axis_grey, size = base_size),
          plot.title = element_text(face = "bold", size = base_size + 1.5, colour = title_grey, hjust = 0.5),
          plot.caption = element_text(colour = axis_grey, hjust = 0.5),
          legend.position = "bottom")
}
p_txt <- function(p) ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p))
season_of <- function(loc) ifelse(grepl("_08_", loc), "Dry", "Rainy")
stream_of <- function(loc) sub("_.*$", "", loc)
wil <- function(x, y) suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE))
index_labels <- c(S_rel = "Richness", PD_rel = "PD", FRic_rel = "FRic")

# Friedman + paired Wilcoxon among the levels of 'group', with complete blocks
friedman_wilcoxon <- function(data, group, levels_) {
  w <- data %>% select(block, {{ group }}, value) %>%
    pivot_wider(names_from = {{ group }}, values_from = value) %>% na.omit()
  mm <- as.matrix(w[, levels_]); ft <- friedman.test(mm)
  global <- data.frame(n = nrow(mm), Chi2 = unname(ft$statistic), df = unname(ft$parameter), p = ft$p.value,
                       W_Kendall = unname(ft$statistic) / (nrow(mm) * (length(levels_) - 1)))
  pairwise <- bind_rows(lapply(combn(levels_, 2, simplify = FALSE), function(pr) {
    tt <- wil(w[[pr[1]]], w[[pr[2]]])
    data.frame(Level1 = pr[1], Level2 = pr[2], n = nrow(w), Median_diff = median(w[[pr[1]]] - w[[pr[2]]]),
               V = unname(tt$statistic), p = tt$p.value)
  }))
  list(global = global, pairwise = pairwise)
}

# -----------------------------------------------------------------------------
# 1. Data per unit
# -----------------------------------------------------------------------------
idx <- read.csv("integrated/results/indices_per_unit.csv") %>%
  mutate(Season = season_of(Locality), Stream = stream_of(Locality))
long <- idx %>%
  select(Method, Locality, Season, Stream, all_of(names(index_labels))) %>%
  pivot_longer(all_of(names(index_labels)), names_to = "Index", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(Index = factor(index_labels[Index], levels = index_labels),
         Method = factor(Method, levels = names(method_cols)),
         Season = factor(Season, levels = names(season_cols)), block = Locality)

descr_season <- long %>% group_by(Season, Index, Method) %>%
  summarise(n = n(), Mean = mean(value), SD = sd(value), Median = median(value), .groups = "drop")

gain <- read.csv("integrated/results/gain_per_unit.csv") %>%
  mutate(Season = factor(season_of(Locality), levels = names(season_cols)), Stream = stream_of(Locality),
         Index = factor(Index, levels = index_labels), block = Locality)

# -----------------------------------------------------------------------------
# 2. A and B: Friedman within each season
# -----------------------------------------------------------------------------
A <- list(); B <- list(); C2 <- list()
for (pe in names(season_cols)) {
  for (ix in index_labels) {
    r <- friedman_wilcoxon(long %>% filter(Season == pe, Index == ix), Method, names(method_cols))
    A[[paste(pe, ix)]] <- list(global = cbind(Season = pe, Index = ix, r$global),
                               pairwise = cbind(Season = pe, Index = ix, r$pairwise))
  }
  for (mt in names(method_cols)) {
    r <- friedman_wilcoxon(long %>% filter(Season == pe, Method == mt), Index, unname(index_labels))
    B[[paste(pe, mt)]] <- list(global = cbind(Season = pe, Method = mt, r$global),
                               pairwise = cbind(Season = pe, Method = mt, r$pairwise))
  }
  r <- friedman_wilcoxon(gain %>% filter(Season == pe) %>% transmute(block, Index, value = gain_eDNA),
                         Index, unname(index_labels))
  C2[[pe]] <- list(global = cbind(Season = pe, r$global), pairwise = cbind(Season = pe, r$pairwise))
}
friedman_categories <- bind_rows(lapply(A, `[[`, "global")); pairwise_categories <- bind_rows(lapply(A, `[[`, "pairwise"))
friedman_indices <- bind_rows(lapply(B, `[[`, "global")); pairwise_indices <- bind_rows(lapply(B, `[[`, "pairwise"))
friedman_gain_eDNA <- bind_rows(lapply(C2, `[[`, "global")); pairwise_gain_eDNA <- bind_rows(lapply(C2, `[[`, "pairwise"))

# C: eDNA vs traditional within each season
gain_methods <- gain %>% group_by(Season, Index) %>%
  summarise(n = n(), Mean_eDNA = mean(gain_eDNA), Mean_Trad = mean(gain_Trad),
            V = unname(wil(gain_eDNA, gain_Trad)$statistic), p = wil(gain_eDNA, gain_Trad)$p.value,
            .groups = "drop")

# -----------------------------------------------------------------------------
# 3. D: dry vs rainy (mean per stream, n = 12)
# -----------------------------------------------------------------------------
dry_vs_rainy <- function(data, ...) {
  data %>% group_by(..., Stream, Season) %>% summarise(value = mean(value), .groups = "drop") %>%
    pivot_wider(names_from = Season, values_from = value) %>% filter(!is.na(Dry), !is.na(Rainy)) %>%
    group_by(...) %>%
    summarise(n = n(), Mean_dry = mean(Dry), Mean_rainy = mean(Rainy), Median_diff = median(Dry - Rainy),
              V = unname(wil(Dry, Rainy)$statistic), p = wil(Dry, Rainy)$p.value, .groups = "drop")
}
seasons_indices <- dry_vs_rainy(long, Index, Method)
seasons_gain <- dry_vs_rainy(gain %>% pivot_longer(c(gain_eDNA, gain_Trad), names_to = "Source", values_to = "value"),
                             Index, Source)

occ <- as.data.frame(read_excel(occ_file, sheet = "3_analysis_input"))
occ$Identification <- trimws(occ$Identification)
pa <- occ %>% distinct(Sampling_Method, Locality, sp = Identification)
jac <- pa %>% group_by(Locality) %>%
  summarise(a = length(intersect(sp[Sampling_Method == "Traditional"], sp[Sampling_Method == "eDNA"])),
            n_union = length(unique(sp)), .groups = "drop") %>%
  mutate(value = a / n_union, Season = factor(season_of(Locality), levels = names(season_cols)),
         Stream = stream_of(Locality), Index = "Jaccard")
seasons_jaccard <- dry_vs_rainy(jac, Index)

# -----------------------------------------------------------------------------
# 4. E: regional (gamma) diversity per method and season
# -----------------------------------------------------------------------------
tip_names <- read.csv("phylogeny/results/taxon_tip_match.csv")
tip <- setNames(tip_names$Tree_tip, trimws(tip_names$Taxon))
tree <- read.tree("phylogeny/data/pool_tree_119taxa.nwk")
coords <- as.matrix(read.csv("functional/raw/data/functional_space_pcoa.csv", row.names = 1, check.names = FALSE))
vol_pool <- convhulln(coords, options = "FA")$vol
pd_pool <- sum(tree$edge.length)
pd_set <- function(sp) pd(matrix(1, 1, length(sp), dimnames = list("a", tip[sp])), tree, include.root = TRUE)$PD
fric_set <- function(sp) convhulln(coords[sp, , drop = FALSE], options = "FA")$vol / vol_pool

pa_season <- pa %>% mutate(Season = season_of(Locality))
sets <- bind_rows(pa_season, pa_season %>% mutate(Sampling_Method = "Combined")) %>% distinct(Sampling_Method, Season, sp)
gamma_season <- sets %>% group_by(Method = Sampling_Method, Season) %>%
  summarise(S = n(), PD = pd_set(sp), FRic = fric_set(sp), .groups = "drop") %>%
  mutate(S_rel = S / nrow(tip_names), PD_rel = PD / pd_pool, FRic_rel = FRic)
turnover_season <- sets %>% group_by(Method = Sampling_Method, sp) %>%
  summarise(Per = paste(sort(unique(Season)), collapse = "+"), .groups = "drop") %>%
  mutate(Per = recode(Per, "Dry+Rainy" = "Both seasons", Dry = "Dry only", Rainy = "Rainy only")) %>%
  count(Method, Per)

# -----------------------------------------------------------------------------
# 5. Save results
# -----------------------------------------------------------------------------
tabs <- c("descr_season", "friedman_categories", "pairwise_categories", "friedman_indices", "pairwise_indices",
          "gain_methods", "friedman_gain_eDNA", "pairwise_gain_eDNA", "seasons_indices", "seasons_gain",
          "seasons_jaccard", "gamma_season", "turnover_season")
wb <- createWorkbook()
for (nm in tabs) {
  write.csv(get(nm), file.path(out_res, paste0("season_", nm, ".csv")), row.names = FALSE)
  addWorksheet(wb, nm); writeData(wb, nm, get(nm))
}
saveWorkbook(wb, file.path(out_res, "season_analyses.xlsx"), overwrite = TRUE)
options(width = 200)
for (nm in tabs) { cat("\n==", nm, "\n"); print(as.data.frame(get(nm)), digits = 4) }

# -----------------------------------------------------------------------------
# 6. Figure 2 - standardized indices by category and season
# -----------------------------------------------------------------------------
summary2 <- long %>% group_by(Index, Method, Season) %>%
  summarise(m = mean(value), ci = qt(0.975, n() - 1) * sd(value) / sqrt(n()), .groups = "drop")
lab2 <- seasons_indices %>% mutate(lab = p_txt(p))
dodge <- position_dodge(width = 0.8)
fig2 <- ggplot(long, aes(Method, value, fill = Season, colour = Season)) +
  geom_violin(position = dodge, alpha = 0.25, trim = TRUE, linewidth = 0.4, scale = "width", width = 0.75,
              key_glyph = draw_key_round) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8, seed = 1),
             size = 0.9, alpha = 0.5, show.legend = FALSE) +
  geom_errorbar(data = summary2, aes(Method, ymin = m - ci, ymax = m + ci, group = Season), inherit.aes = FALSE,
                position = dodge, width = 0.12, linewidth = 0.6, colour = axis_grey) +
  geom_point(data = summary2, aes(Method, m, group = Season, fill = Season), inherit.aes = FALSE,
             position = dodge, shape = 21, size = 2.6, colour = axis_grey, show.legend = FALSE) +
  geom_text(data = lab2, aes(Method, Inf, label = lab), inherit.aes = FALSE, vjust = 1.2, size = 2.8,
            colour = axis_grey) +
  facet_wrap(~Index, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = season_cols, name = "Season") +
  scale_colour_manual(values = season_cols, name = "Season") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"), expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = NULL, y = "Proportion of the regional pool",
       title = "Diversity Recovered by Each Sampling Category in the Dry and Rainy Seasons") +
  fig_theme(10)
ggsave(file.path(out_fig, "fig2_season_indices.png"), fig2, width = 10, height = 4.8, dpi = 300, bg = "white")
ggsave(file.path(out_fig, "fig2_season_indices.pdf"), fig2, width = 10, height = 4.8, device = cairo_pdf)

# -----------------------------------------------------------------------------
# 7. Figure 3 - exclusive contribution of each method by season
# -----------------------------------------------------------------------------
g3 <- gain %>%
  pivot_longer(c(gain_eDNA, gain_Trad), names_to = "Source", values_to = "gain") %>%
  mutate(Source = factor(recode(Source, gain_eDNA = "Added by eDNA", gain_Trad = "Added by traditional"),
                         levels = names(source_cols)))
summary3 <- g3 %>% group_by(Index, Season, Source) %>%
  summarise(m = mean(gain), ci = qt(0.975, n() - 1) * sd(gain) / sqrt(n()), .groups = "drop")
lab3 <- gain_methods %>% mutate(lab = p_txt(p))
fig3 <- ggplot(g3, aes(Season, gain, fill = Source, colour = Source)) +
  geom_violin(position = dodge, alpha = 0.25, trim = TRUE, linewidth = 0.4, scale = "width", width = 0.75,
              key_glyph = draw_key_round) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8, seed = 1),
             size = 0.9, alpha = 0.5, show.legend = FALSE) +
  geom_errorbar(data = summary3, aes(Season, ymin = m - ci, ymax = m + ci, group = Source), inherit.aes = FALSE,
                position = dodge, width = 0.12, linewidth = 0.6, colour = axis_grey) +
  geom_point(data = summary3, aes(Season, m, group = Source, fill = Source), inherit.aes = FALSE,
             position = dodge, shape = 21, size = 2.6, colour = axis_grey, show.legend = FALSE) +
  geom_text(data = lab3, aes(Season, 1.06, label = lab), inherit.aes = FALSE, size = 2.8, colour = axis_grey) +
  facet_wrap(~Index, nrow = 1) +
  scale_fill_manual(values = source_cols, name = NULL) +
  scale_colour_manual(values = source_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2), labels = function(x) paste0(x * 100, "%")) +
  labs(x = NULL, y = "Share of the combined diversity\ncontributed exclusively by each method",
       title = "Exclusive Contribution of Each Sampling Method in the Dry and Rainy Seasons") +
  fig_theme(10)
ggsave(file.path(out_fig, "fig3_season_complementarity.png"), fig3, width = 10, height = 4.8, dpi = 300, bg = "white")
ggsave(file.path(out_fig, "fig3_season_complementarity.pdf"), fig3, width = 10, height = 4.8, device = cairo_pdf)

# -----------------------------------------------------------------------------
# 8. Figure S4 - taxa by season, similarity between methods and regional diversity
# -----------------------------------------------------------------------------
turnover_cols <- c("Dry only" = season_cols[["Dry"]], "Both seasons" = "grey60", "Rainy only" = season_cols[["Rainy"]])
sA <- turnover_season %>% mutate(Method = factor(Method, levels = names(method_cols)),
                                 Per = factor(Per, levels = names(turnover_cols)))
pA <- ggplot(sA, aes(Method, n, fill = Per)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3, key_glyph = draw_key_round) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3, colour = "white") +
  scale_fill_manual(values = turnover_cols, name = NULL) +
  labs(x = NULL, y = "Number of taxa", title = "Taxa Recorded in Each Season") +
  fig_theme(10)

jac_stream <- jac %>% group_by(Stream, Season) %>% summarise(value = mean(value), .groups = "drop")
mj <- jac %>% group_by(Season) %>%
  summarise(m = mean(value), ci = qt(0.975, n() - 1) * sd(value) / sqrt(n()), .groups = "drop")
pB <- ggplot(jac, aes(Season, value, fill = Season, colour = Season)) +
  geom_line(data = jac_stream, mapping = aes(Season, value, group = Stream), colour = "grey75", linewidth = 0.3,
            inherit.aes = FALSE) +
  geom_violin(alpha = 0.25, linewidth = 0.4, width = 0.6, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.08, seed = 1), size = 0.9, alpha = 0.5, show.legend = FALSE) +
  geom_errorbar(data = mj, aes(Season, ymin = m - ci, ymax = m + ci), inherit.aes = FALSE, width = 0.1,
                linewidth = 0.6, colour = axis_grey) +
  geom_point(data = mj, aes(Season, m, fill = Season), inherit.aes = FALSE, shape = 21, size = 2.6,
             colour = axis_grey, show.legend = FALSE) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.3, size = 2.8, colour = axis_grey,
           label = paste0("Wilcoxon ", p_txt(seasons_jaccard$p), " (n = ", seasons_jaccard$n, " streams)")) +
  scale_fill_manual(values = season_cols) + scale_colour_manual(values = season_cols) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Jaccard similarity between methods", title = "Similarity Between Methods") +
  fig_theme(10)

gC <- gamma_season %>% select(Method, Season, S_rel, PD_rel, FRic_rel) %>%
  pivot_longer(c(S_rel, PD_rel, FRic_rel), names_to = "Index", values_to = "value") %>%
  mutate(Index = factor(index_labels[Index], levels = index_labels), Method = factor(Method, levels = names(method_cols)),
         Season = factor(Season, levels = names(season_cols)))
pC <- ggplot(gC, aes(Method, value, fill = Season)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "white", linewidth = 0.3,
           key_glyph = draw_key_round) +
  facet_wrap(~Index, nrow = 1) +
  scale_fill_manual(values = season_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(x * 100, "%")) +
  labs(x = NULL, y = "Proportion of the regional pool", title = "Regional Diversity Recovered in Each Season") +
  fig_theme(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

figS4 <- (pA | pB) / pC + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", colour = title_grey))
ggsave(file.path(out_fig, "figS4_season_regional.png"), figS4, width = 10, height = 9, dpi = 300, bg = "white")
ggsave(file.path(out_fig, "figS4_season_regional.pdf"), figS4, width = 10, height = 9, device = cairo_pdf)
