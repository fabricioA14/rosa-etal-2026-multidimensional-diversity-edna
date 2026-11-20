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
# Script:  03_functional_richness_FRic.R
# Purpose: Functional richness (FRic, incidence-based) per stream-campaign unit
#          and sampling category. Run with the argument "raw" (the 33
#          morphological measurements; main analysis) or "ratios"
#          (ecomorphological ratios + size; Villéger et al. 2010; Brosse et
#          al. 2021). Steps: log(x + 1) + z-score of all traits (set
#          NO_LOG=1 for untransformed traits); Spearman correlation filter
#          |rho| >= 0.7 (Table S2); Gower distance -> PCoA, number of axes by
#          mAD (Maire et al. 2015; Magneville et al. 2022); FRic = convex hull
#          volume / volume of the pool (Villéger et al. 2008); SES.FRic with
#          999 randomizations of the pool.
# Input:   pool_taxa_traits.xlsx (morphological traits)
#          SuppInfo_S1.xlsx, sheet 3_analysis_input
# Output:  functional/<version>/{data,results,figures}
# =============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(readxl)
  library(ape)
  library(FD)
  library(geometry)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

set.seed(2024)
n_runs  <- 999
cor_threshold <- 0.7

pool_file <- "pool_taxa_traits.xlsx"
occ_file  <- "SuppInfo_S1.xlsx"
version <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(version)) version <- "ratios"
stopifnot(version %in% c("ratios", "raw"))

out_data <- file.path("functional", version, "data")
out_res  <- file.path("functional", version, "results")
out_fig  <- file.path("functional", version, "figures")
for (d in c(out_data, out_res, out_fig)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
version_label <- c(ratios = "ecomorphological ratios", raw = "raw measurements")[[version]]

# -----------------------------------------------------------------------------
# 1. Traits of the selected version
# -----------------------------------------------------------------------------
pool <- read.xlsx(pool_file, sheet = 1)
names(pool) <- trimws(names(pool))
pool$Taxon <- trimws(pool$Taxon)
m <- pool
rownames(m) <- m$Taxon

if (version == "raw") {
  traits <- m[, setdiff(names(m), c("Taxon", "Tree_tip"))]
  rownames(traits) <- m$Taxon
  definitions <- data.frame(Trait = names(traits), Formula = "raw measurement", Function = NA)
  # size trait protected in the correlation filter
  protected <- "standard_length"
  no_log <- character(0)
} else {

# safe ratio: zero denominator -> 0 (absent structure)
r <- function(a, b) ifelse(b == 0, 0, a / b)

traits <- with(m, data.frame(
  # size
  Size_logSL              = log10(standard_length),
  Mass_logBM              = log10(body_mass),
  # feeding / mouth
  Mouth_surface           = r(mouth_width * mouth_gape, body_width * max_body_depth),
  Mouth_shape             = r(mouth_gape, mouth_width),
  Mouth_position          = r(mouth_height, head_depth),
  Protrusion_rel          = r(mouth_protrusion, head_length),
  Snout_rel               = r(snout_length, head_length),
  Barbel_rel              = r(longest_barbel_length, head_length),
  Gut_rel                 = r(gut_length, standard_length),
  Head_length_rel         = r(head_length, standard_length),
  # sensory
  Eye_size                = r(eye_diameter, head_depth),
  Eye_position            = r(eye_position, head_depth),
  # body shape / habitat use
  Body_depth_rel          = r(max_body_depth, standard_length),
  Body_compression        = r(max_body_depth, body_width),
  Head_depth_rel          = r(head_depth, max_body_depth),
  Peduncle_depth_rel      = r(caudal_peduncle_depth, max_body_depth),
  # locomotion / fins
  Pectoral_position       = r(pectoral_fin_position, body_depth_at_pectoral_insertion),
  Pectoral_aspect         = r(pectoral_fin_length^2, pectoral_fin_area),
  Pectoral_length_rel     = r(pectoral_fin_length, standard_length),
  Caudal_aspect           = r(caudal_fin_height^2, caudal_fin_area),
  Caudal_throttling       = r(caudal_fin_height, caudal_peduncle_depth),
  Fins_surface            = r(2 * pectoral_fin_area + caudal_fin_area, body_area),
  Dorsal_length_rel       = r(dorsal_fin_length, standard_length),
  Anal_length_rel         = r(anal_fin_length, standard_length),
  Pelvic_length_rel       = r(pelvic_fin_length, standard_length),
  row.names = m$Taxon
))
stopifnot(all(is.finite(as.matrix(traits))))

definitions <- data.frame(
  Trait = names(traits),
  Formula = c(
    "log10(standard_length)", "log10(body_mass)",
    "(mouth_width x mouth_gape) / (body_width x max_body_depth)",
    "mouth_gape / mouth_width", "mouth_height / head_depth",
    "mouth_protrusion / head_length", "snout_length / head_length",
    "longest_barbel_length / head_length",
    "gut_length / standard_length", "head_length / standard_length",
    "eye_diameter / head_depth", "eye_position / head_depth",
    "max_body_depth / standard_length", "max_body_depth / body_width",
    "head_depth / max_body_depth", "caudal_peduncle_depth / max_body_depth",
    "pectoral_fin_position / body_depth_at_pectoral_insertion",
    "pectoral_fin_length² / pectoral_fin_area",
    "pectoral_fin_length / standard_length",
    "caudal_fin_height² / caudal_fin_area",
    "caudal_fin_height / caudal_peduncle_depth",
    "(2 x pectoral_fin_area + caudal_fin_area) / body_area",
    "dorsal_fin_length / standard_length",
    "anal_fin_length / standard_length",
    "pelvic_fin_length / standard_length"),
  Function = c("Size", "Size", rep("Feeding", 8), rep("Sensory", 2),
               rep("Body shape / habitat", 4), rep("Locomotion", 9))
)
protected <- "Size_logSL"
# size and mass are already log10-transformed
no_log <- c("Size_logSL", "Mass_logBM")
}

# -----------------------------------------------------------------------------
# 2. Standardization: log(x+1) + z-score
#    In the 'raw' version (used in the manuscript) all 33 traits have
#    |skewness| > 1 and are therefore all log(x+1)-transformed; the skewness
#    criterion below only changes the selection in the 'ratios' version.
# -----------------------------------------------------------------------------
skewness <- function(x) mean((x - mean(x))^3) / sd(x)^3
skew_before <- sapply(traits, skewness)
tr_log <- traits
to_log <- names(skew_before)[abs(skew_before) > 1 & !names(skew_before) %in% no_log]
# sensitivity analysis (Table S3, section G): NO_LOG=1 turns off the log transformation
if (Sys.getenv("NO_LOG") == "1") to_log <- character(0)
tr_log[to_log] <- lapply(tr_log[to_log], log1p)
skew_after <- sapply(tr_log, skewness)
tr_std <- as.data.frame(scale(tr_log))

definitions$Skewness_raw <- round(skew_before, 2)
definitions$Transformation <- ifelse(definitions$Trait %in% to_log, "log(x+1) + z-score", "z-score")
definitions$Skewness_final <- round(skew_after, 2)

# -----------------------------------------------------------------------------
# 3. Correlation filter (Spearman |rho| >= 0.7)
#    In each pair above the threshold, the trait with the higher mean |rho| with the others is removed.
#    The size trait (standard length) is always kept.
# -----------------------------------------------------------------------------
rho_all <- cor(tr_std, method = "spearman")
keep <- colnames(tr_std)
removed <- data.frame()
repeat {
  rr <- abs(rho_all[keep, keep]); diag(rr) <- 0
  if (max(rr) < cor_threshold) break
  ij <- which(rr == max(rr), arr.ind = TRUE)[1, ]
  pair <- keep[ij]
  mean_rho <- rowMeans(rr[pair, , drop = FALSE])
  out <- if (protected %in% pair) setdiff(pair, protected) else pair[which.max(mean_rho)]
  stays <- setdiff(pair, out)
  removed <- rbind(removed, data.frame(Removed = out, Correlated_with = stays,
                                       rho = rho_all[out, stays]))
  keep <- setdiff(keep, out)
}
definitions$Retained <- definitions$Trait %in% keep
tr_final <- tr_std[, keep]

# -----------------------------------------------------------------------------
# 4. Gower -> PCoA -> quality (mAD) to choose the number of axes
# -----------------------------------------------------------------------------
gow <- gowdis(tr_final)
pc <- pcoa(gow)
eig_pos <- pc$values$Eigenvalues[pc$values$Eigenvalues > 0]
var_expl <- eig_pos / sum(eig_pos)
coords_all <- pc$vectors[, seq_along(eig_pos)]

d_tr <- as.vector(gow) / max(gow)
quality <- bind_rows(lapply(2:8, function(k) {
  d_sp <- as.vector(dist(coords_all[, 1:k]))
  d_sp <- d_sp / max(d_sp)
  data.frame(Axes = k, mAD = mean(abs(d_sp - d_tr)), mSD = mean((d_sp - d_tr)^2),
             Cumulative_variance = sum(var_expl[1:k]))
}))
# number of axes: lowest mAD, limited to 5 so that small communities have FRic
k <- quality$Axes[which.min(ifelse(quality$Axes <= 5, quality$mAD, Inf))]
coords <- coords_all[, 1:k]
colnames(coords) <- paste0("PC", 1:k)

cor_axes <- cor(tr_final, coords, method = "spearman")

# -----------------------------------------------------------------------------
# 5. Communities and FRic
# -----------------------------------------------------------------------------
occ <- as.data.frame(read_excel(occ_file, sheet = "3_analysis_input"))
occ$Identification <- trimws(occ$Identification)
stopifnot(all(occ$Identification %in% rownames(coords)))

pa <- occ %>% distinct(Sampling_Method, Locality, sp = Identification)
pa_all <- bind_rows(pa, pa %>% distinct(Locality, sp) %>% mutate(Sampling_Method = "Combined"))

vol_pool <- convhulln(coords, options = "FA")$vol
fric <- function(sp) {
  if (length(sp) <= k) return(NA_real_)
  v <- tryCatch(convhulln(coords[sp, , drop = FALSE], options = "FA")$vol, error = function(e) NA_real_)
  v / vol_pool
}

comms <- pa_all %>%
  group_by(Method = Sampling_Method, Locality) %>%
  summarise(sp = list(sp), S = n(), .groups = "drop")
comms$FRic <- sapply(comms$sp, fric)

# -----------------------------------------------------------------------------
# 6. SES.FRic: null distribution per richness (random draw of S taxa from the 119-taxa pool)
# -----------------------------------------------------------------------------
null_by_S <- function(S) replicate(n_runs, fric(sample(rownames(coords), S)))
S_unique <- sort(unique(c(comms$S[comms$S > k], 74, 82)))
null_dists <- setNames(lapply(S_unique, null_by_S), S_unique)

ses <- function(obs, S) {
  if (is.na(obs) || !(as.character(S) %in% names(null_dists))) return(c(NA, NA, NA, NA))
  v <- null_dists[[as.character(S)]]; v <- v[!is.na(v)]
  c(mean(v), sd(v), (obs - mean(v)) / sd(v), (sum(v < obs) + 0.5 * sum(v == obs)) / length(v))
}
s_mat <- t(mapply(ses, comms$FRic, comms$S))
comms$FRic_null_mean <- s_mat[, 1]
comms$FRic_null_sd <- s_mat[, 2]
comms$SES_FRic <- s_mat[, 3]
comms$p_FRic <- s_mat[, 4]

alpha <- comms %>%
  select(-sp) %>%
  mutate(Stream = sub("_.*$", "", Locality),
         Campaign = sub("^[^_]+_", "", Locality),
         Method = factor(Method, levels = c("Traditional", "eDNA", "Combined"))) %>%
  arrange(Method, Stream, Campaign)

# Gamma per method
sp_trad <- unique(pa$sp[pa$Sampling_Method == "Traditional"])
sp_edna <- unique(pa$sp[pa$Sampling_Method == "eDNA"])
sp_comb <- union(sp_trad, sp_edna)
gamma <- data.frame(Set = c("Traditional", "eDNA", "Combined"),
                    S = c(length(sp_trad), length(sp_edna), length(sp_comb)))
gamma$FRic <- c(fric(sp_trad), fric(sp_edna), fric(sp_comb))
g_ses <- t(mapply(ses, gamma$FRic[1:2], gamma$S[1:2]))
gamma$SES_FRic <- c(g_ses[, 3], NA)
gamma$p_FRic <- c(g_ses[, 4], NA)

# -----------------------------------------------------------------------------
# 7. Traditional vs eDNA vs Combined (paired by sample)
# -----------------------------------------------------------------------------
methods <- c("Traditional", "eDNA", "Combined")
vars <- c(S = "Richness (S)", FRic = "FRic (proportion of the pool)", SES_FRic = "SES.FRic")
long <- alpha %>%
  select(Method, Locality, all_of(names(vars))) %>%
  pivot_longer(all_of(names(vars)), names_to = "Metric") %>%
  group_by(Metric, Locality) %>% filter(sum(!is.na(value)) == 3) %>% ungroup()

friedman <- bind_rows(lapply(names(vars), function(v) {
  mm <- long %>% filter(Metric == v) %>% pivot_wider(names_from = Method, values_from = value) %>%
    select(all_of(methods)) %>% as.matrix()
  ft <- friedman.test(mm)
  data.frame(Metric = v, n_blocks = nrow(mm), Chi2 = unname(ft$statistic), p = ft$p.value,
             W_Kendall = unname(ft$statistic) / (nrow(mm) * 2))
}))
friedman$p_holm <- p.adjust(friedman$p, "holm")

pairs <- list(c("Traditional", "eDNA"), c("Traditional", "Combined"), c("eDNA", "Combined"))
posthoc <- bind_rows(lapply(names(vars), function(v) {
  w <- long %>% filter(Metric == v) %>% pivot_wider(names_from = Method, values_from = value)
  res <- bind_rows(lapply(pairs, function(pr) {
    x <- w[[pr[1]]]; y <- w[[pr[2]]]
    tt <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE))
    data.frame(Metric = v, Group1 = pr[1], Group2 = pr[2], n = length(x),
               Median_diff = median(x - y), V = unname(tt$statistic), p = tt$p.value)
  }))
  res$p_holm <- p.adjust(res$p, "holm")
  res
}))

differs <- function(v, a, b) {
  r0 <- posthoc %>% filter(Metric == v, (Group1 == a & Group2 == b) | (Group1 == b & Group2 == a))
  r0$p_holm < 0.05
}
letters_tab <- bind_rows(lapply(names(vars), function(v) {
  ord <- long %>% filter(Metric == v) %>% group_by(Method) %>% summarise(md = median(value)) %>%
    arrange(desc(md)) %>% pull(Method) %>% as.character()
  groups <- list()
  for (g in ord) {
    fits <- vapply(groups, function(G) !any(vapply(G, function(h) differs(v, h, g), logical(1))), logical(1))
    if (any(fits)) {
      for (i in which(fits)) groups[[i]] <- c(groups[[i]], g)
    } else {
      prev <- ord[seq_len(match(g, ord) - 1)]
      groups[[length(groups) + 1]] <- c(prev[!vapply(prev, function(h) differs(v, h, g), logical(1))], g)
    }
  }
  if (friedman$p_holm[friedman$Metric == v] >= 0.05) groups <- list(ord)
  data.frame(Metric = v, Method = ord, row.names = NULL,
             Letter = sapply(ord, function(g) paste(letters[which(sapply(groups, function(G) g %in% G))], collapse = "")))
}))

descriptive <- long %>% group_by(Metric, Method) %>%
  summarise(n = n(), Mean = mean(value), SD = sd(value), Median = median(value), .groups = "drop")

signif_structure <- alpha %>% group_by(Method) %>%
  summarise(n_with_FRic = sum(!is.na(SES_FRic)),
            FRic_lower = sum(p_FRic <= 0.025, na.rm = TRUE),
            FRic_higher = sum(p_FRic >= 0.975, na.rm = TRUE), .groups = "drop")

# -----------------------------------------------------------------------------
# 8. Export
# -----------------------------------------------------------------------------
write.csv(traits, file.path(out_data, "traits_unstandardized.csv"))
write.csv(tr_std, file.path(out_data, "traits_standardized_all.csv"))
write.csv(tr_final, file.path(out_data, "traits_standardized_retained.csv"))
write.csv(coords, file.path(out_data, "functional_space_pcoa.csv"))

write.csv(definitions, file.path(out_res, "trait_definitions_selection.csv"), row.names = FALSE)
write.csv(removed, file.path(out_res, "traits_removed_correlation.csv"), row.names = FALSE)
write.csv(quality, file.path(out_res, "functional_space_quality.csv"), row.names = FALSE)
write.csv(cor_axes, file.path(out_res, "trait_axis_correlation.csv"))
write.csv(alpha, file.path(out_res, "communities_fric_ses.csv"), row.names = FALSE)
write.csv(gamma, file.path(out_res, "gamma_fric_by_method.csv"), row.names = FALSE)
write.csv(friedman, file.path(out_res, "comparison3_friedman.csv"), row.names = FALSE)
write.csv(posthoc, file.path(out_res, "comparison3_posthoc_wilcoxon.csv"), row.names = FALSE)

wb <- createWorkbook()
for (nm in c("definitions", "removed", "quality", "alpha", "gamma", "descriptive", "friedman", "posthoc",
             "letters_tab", "signif_structure")) {
  addWorksheet(wb, nm)
  writeData(wb, nm, get(nm))
}
addWorksheet(wb, "cor_axes"); writeData(wb, "cor_axes", as.data.frame(cor_axes), rowNames = TRUE)
saveWorkbook(wb, file.path(out_res, "functional_diversity_fric.xlsx"), overwrite = TRUE)

# -----------------------------------------------------------------------------
# 9. Figures
# -----------------------------------------------------------------------------
method_cols <- c(Traditional = "#d1495b", eDNA = "#2a78b5", Combined = "grey45")
fig_theme <- theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank())
save_fig <- function(p, name, w, h) {
  ggsave(file.path(out_fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(out_fig, paste0(name, ".pdf")), p, width = w, height = h)
}

# Fig 1 - correlation matrix of the standardized traits
trait_order <- colnames(rho_all)[hclust(as.dist(1 - abs(rho_all)))$order]
d1 <- as.data.frame(as.table(rho_all)) %>%
  mutate(Var1 = factor(Var1, levels = trait_order), Var2 = factor(Var2, levels = trait_order))
trait_labels <- setNames(ifelse(trait_order %in% keep, trait_order, paste0(trait_order, " (removed)")), trait_order)
f1 <- ggplot(d1, aes(Var1, Var2, fill = Freq)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_point(data = d1 %>% filter(abs(Freq) >= cor_threshold, Var1 != Var2), shape = 4, size = 1, colour = "grey20") +
  scale_fill_gradient2(low = "#2a78b5", mid = "white", high = "#d1495b", limits = c(-1, 1), name = "rho") +
  scale_x_discrete(labels = trait_labels) + scale_y_discrete(labels = trait_labels) +
  labs(x = NULL, y = NULL,
       title = sprintf("Spearman correlation among the %d standardized traits [%s]", ncol(tr_std), version_label),
       subtitle = sprintf("x in cells: |rho| >= %.1f; %d traits retained", cor_threshold, length(keep))) +
  coord_fixed() + fig_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 7))
save_fig(f1, "fig1_trait_correlation", 8.5, 8)

# Fig 2 - quality of the functional space
f2 <- ggplot(quality, aes(Axes, mAD)) +
  geom_line(colour = "grey40") + geom_point(size = 2.5, colour = "grey40") +
  geom_point(data = quality %>% filter(Axes == k), size = 4, colour = "#d1495b") +
  geom_text(aes(label = sprintf("%.0f%%", 100 * Cumulative_variance)), vjust = -1, size = 2.8) +
  scale_x_continuous(breaks = 2:8) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  labs(x = "Number of PCoA axes", y = "mAD (mean absolute deviation)",
       title = sprintf("Quality of the functional space [%s]", version_label),
       subtitle = sprintf("Selected: %d axes (red). Label = cumulative variance.", k)) + fig_theme
save_fig(f2, "fig2_space_quality", 5.5, 4)

# Fig 3 - functional space with the gamma convex hull of each method
cd <- as.data.frame(coords) %>% mutate(Taxon = rownames(coords),
  Group = case_when(Taxon %in% intersect(sp_trad, sp_edna) ~ "Both",
                    Taxon %in% sp_trad ~ "Traditional only", TRUE ~ "eDNA only"))
hull_df <- function(sp, ax, lab) {
  xy <- coords[sp, ax]; h <- chull(xy)
  data.frame(x = xy[h, 1], y = xy[h, 2], Set = lab)
}
panel <- function(ax) {
  hh <- bind_rows(hull_df(rownames(coords), ax, "Pool (119)"),
                  hull_df(sp_trad, ax, "Traditional"), hull_df(sp_edna, ax, "eDNA"))
  ggplot() +
    geom_polygon(data = hh, aes(x, y, colour = Set, fill = Set), alpha = 0.08, linewidth = 0.7) +
    geom_point(data = cd, aes(.data[[colnames(coords)[ax[1]]]], .data[[colnames(coords)[ax[2]]]], shape = Group),
               size = 1.6, alpha = 0.8) +
    scale_colour_manual(values = c("Pool (119)" = "grey30", Traditional = "#d1495b", eDNA = "#2a78b5")) +
    scale_fill_manual(values = c("Pool (119)" = "grey30", Traditional = "#d1495b", eDNA = "#2a78b5")) +
    scale_shape_manual(values = c(Both = 16, `Traditional only` = 17, `eDNA only` = 15)) +
    labs(x = sprintf("%s (%.1f%%)", colnames(coords)[ax[1]], 100 * var_expl[ax[1]]),
         y = sprintf("%s (%.1f%%)", colnames(coords)[ax[2]], 100 * var_expl[ax[2]]),
         colour = NULL, fill = NULL, shape = NULL) +
    fig_theme + coord_equal()
}
f3 <- (panel(c(1, 2)) + panel(c(3, 4))) + plot_layout(guides = "collect") +
  plot_annotation(title = sprintf("Functional space (PCoA on Gower) and convex hulls by method [%s]", version_label),
                  subtitle = sprintf("Gamma FRic: Traditional = %.2f, eDNA = %.2f, Combined = %.2f of the pool volume (%dD)",
                                     gamma$FRic[1], gamma$FRic[2], gamma$FRic[3], k)) &
  theme(legend.position = "bottom")
save_fig(f3, "fig3_functional_space", 11, 6)

# Fig 4 - correlation between traits and axes
d4 <- as.data.frame(as.table(cor_axes)) %>%
  mutate(Var1 = factor(Var1, levels = rev(keep)))
f4 <- ggplot(d4, aes(Var2, Var1, fill = Freq)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", Freq)), size = 2.6) +
  scale_fill_gradient2(low = "#2a78b5", mid = "white", high = "#d1495b", limits = c(-1, 1), name = "rho") +
  labs(x = NULL, y = NULL, title = sprintf("Spearman correlation between traits and PCoA axes [%s]", version_label)) +
  fig_theme
save_fig(f4, "fig4_trait_axes", 6, 6)

# Fig 5 - Traditional vs eDNA vs Combined (violin + mean ± 95% CI)
d5 <- long %>% mutate(Metric = factor(vars[Metric], levels = vars),
                      Method = factor(Method, levels = methods))
means <- d5 %>% group_by(Metric, Method) %>%
  summarise(m = mean(value), n = n(), se = sd(value) / sqrt(n()), .groups = "drop") %>%
  mutate(lo = m - qt(0.975, n - 1) * se, hi = m + qt(0.975, n - 1) * se)
lab_let <- letters_tab %>%
  left_join(long %>% group_by(Metric) %>% summarise(top = max(value), amp = diff(range(value))), by = "Metric") %>%
  mutate(Metric = factor(vars[Metric], levels = vars), Method = factor(Method, levels = methods),
         y = top + 0.08 * amp)
lab_fr <- friedman %>% mutate(Metric = factor(vars[Metric], levels = vars),
  lab = sprintf("Friedman X² = %.1f, %s; W = %.2f (n = %d)", Chi2,
                ifelse(p_holm < 0.001, "p < 0.001", sprintf("p = %.3f", p_holm)), W_Kendall, n_blocks))
f5 <- ggplot(d5, aes(Method, value)) +
  geom_hline(data = data.frame(Metric = factor(vars["SES_FRic"], levels = vars)),
             aes(yintercept = 0), linetype = 2, colour = "grey60") +
  geom_violin(aes(fill = Method), colour = NA, alpha = 0.35, trim = TRUE, scale = "width", width = 0.8) +
  geom_line(aes(group = Locality), colour = "grey80", linewidth = 0.25) +
  geom_jitter(aes(colour = Method), width = 0.08, size = 0.7, alpha = 0.5) +
  geom_line(data = means, aes(Method, m, group = 1), colour = "black", linewidth = 0.5) +
  geom_errorbar(data = means, aes(Method, ymin = lo, ymax = hi), inherit.aes = FALSE, width = 0.12, linewidth = 0.7) +
  geom_point(data = means, aes(Method, m, fill = Method), shape = 21, size = 3, stroke = 0.8) +
  geom_text(data = lab_let, aes(Method, y, label = Letter), fontface = "bold", size = 3.6) +
  geom_text(data = lab_fr, aes(x = 2, y = -Inf, label = lab), vjust = -0.6, size = 2.6, colour = "grey30") +
  facet_wrap(~Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = method_cols, guide = "none") +
  scale_colour_manual(values = method_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.1))) +
  labs(x = NULL, y = NULL,
       title = sprintf("Functional richness per sample: Traditional, eDNA and Combined [%s]", version_label),
       caption = paste0("Samples with S > ", k, " in all three categories. Point = mean ± 95% CI; grey lines link the same sample.\n",
                        "Letters: paired Wilcoxon (Holm) after a significant Friedman test; same letter = no difference.")) +
  fig_theme + theme(panel.grid.major.x = element_blank())
save_fig(f5, "fig5_trad_edna_combined", 10, 4.5)

# -----------------------------------------------------------------------------
cat("\nLog-transformed traits:", paste(to_log, collapse = ", "), "\n")
cat("\n== Removed ==\n"); print(removed)
cat("\nRetained (", length(keep), "):", paste(keep, collapse = ", "), "\n")
cat("\n== Quality ==\n"); print(quality); cat("Axes used:", k, "\n")
cat("\n== Gamma ==\n"); print(gamma)
cat("\n== Friedman ==\n"); print(friedman)
cat("\n== Post hoc ==\n"); print(posthoc)
cat("\n== Descriptive ==\n"); print(descriptive)
cat("\n== Significant structure ==\n"); print(signif_structure)
