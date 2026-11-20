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
# Script:  02_phylogenetic_diversity_PD.R
# Purpose: Faith's PD per stream-campaign unit and sampling category
#          (traditional, eDNA, combined); 12 streams x 4 campaigns.
# Input:   SuppInfo_S1.xlsx, sheet 3_analysis_input (occurrences)
#          pool_taxa_traits.xlsx (column Tree_tip links taxa to tree tips)
#          phylogeny/data/3167_Neotrop_Freshwater_Fish_spp_time_tree.tre (from 01)
# Output:  phylogeny/{results,figures}
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(picante)
  library(readxl)
  library(openxlsx)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

set.seed(2024)
n_runs <- 999

pool_file <- "pool_taxa_traits.xlsx"
occ_file  <- "SuppInfo_S1.xlsx"
f_tree    <- "phylogeny/data/3167_Neotrop_Freshwater_Fish_spp_time_tree.tre"
out_res   <- "phylogeny/results"
out_fig   <- "phylogeny/figures"

# -----------------------------------------------------------------------------
# 1. Data and name -> tree tip link
# -----------------------------------------------------------------------------
pool <- read.xlsx(pool_file, sheet = 1)
pool$Taxon <- trimws(pool$Taxon)
pool$Tree_tip <- trimws(pool$Tree_tip)
taxon2tip <- setNames(pool$Tree_tip, pool$Taxon)

tree <- read.nexus(f_tree)
stopifnot(all(pool$Tree_tip %in% tree$tip.label))
# The local pool (119 taxa) is the regional pool of the null model
pool_tree <- keep.tip(tree, pool$Tree_tip)

occ <- as.data.frame(read_excel(occ_file, sheet = "3_analysis_input"))
occ$Identification <- trimws(occ$Identification)
stopifnot(all(occ$Identification %in% pool$Taxon))
occ$tip <- unname(taxon2tip[occ$Identification])

# Presence/absence per sample (Locality = stream_campaign) and method
pa <- occ %>%
  distinct(Sampling_Method, Locality, tip) %>%
  mutate(pres = 1L)
pa_comb <- pa %>% distinct(Locality, tip) %>% mutate(Sampling_Method = "Combined", pres = 1L)
pa_all <- bind_rows(pa, pa_comb) %>%
  mutate(sample = paste(Sampling_Method, Locality, sep = "|"))

comm <- pa_all %>%
  select(sample, tip, pres) %>%
  pivot_wider(names_from = tip, values_from = pres, values_fill = 0L) %>%
  as.data.frame()
rownames(comm) <- comm$sample
comm <- as.matrix(comm[, -1])
# include pool species absent from all samples (zero columns)
missing_tips <- setdiff(pool_tree$tip.label, colnames(comm))
if (length(missing_tips)) {
  comm <- cbind(comm, matrix(0L, nrow(comm), length(missing_tips), dimnames = list(NULL, missing_tips)))
}
comm <- comm[, pool_tree$tip.label]

dist_pool <- cophenetic(pool_tree)

# -----------------------------------------------------------------------------
# 2. Alpha metrics + SES (taxa.labels null model on the 119-taxa pool)
# -----------------------------------------------------------------------------
pd_obs   <- pd(comm, pool_tree, include.root = TRUE)
ses_pd   <- ses.pd(comm, pool_tree, null.model = "taxa.labels", runs = n_runs, include.root = TRUE)
ses_mpd  <- ses.mpd(comm, dist_pool, null.model = "taxa.labels", runs = n_runs)
ses_mntd <- ses.mntd(comm, dist_pool, null.model = "taxa.labels", runs = n_runs)

alpha <- data.frame(
  sample = rownames(comm),
  S = pd_obs$SR,
  PD = pd_obs$PD,
  PD_null_mean = ses_pd$pd.rand.mean,
  SES_PD = ses_pd$pd.obs.z,
  p_PD = ses_pd$pd.obs.p,
  MPD = ses_mpd$mpd.obs,
  SES_MPD = ses_mpd$mpd.obs.z,
  NRI = -ses_mpd$mpd.obs.z,
  p_MPD = ses_mpd$mpd.obs.p,
  MNTD = ses_mntd$mntd.obs,
  SES_MNTD = ses_mntd$mntd.obs.z,
  NTI = -ses_mntd$mntd.obs.z,
  p_MNTD = ses_mntd$mntd.obs.p
) %>%
  separate(sample, into = c("Method", "Locality"), sep = "\\|", remove = FALSE) %>%
  mutate(
    Stream = sub("_.*$", "", Locality),
    Group = substr(Stream, 1, 2),
    Campaign = sub("^[^_]+_", "", Locality),
    Month = ifelse(substr(Campaign, 1, 2) == "08", "August", "February"),
    Year = paste0("20", substr(Campaign, 4, 5)),
    Method = factor(Method, levels = c("Traditional", "eDNA", "Combined"))
  ) %>%
  arrange(Method, Stream, Year, Month)

# -----------------------------------------------------------------------------
# 3. Traditional vs eDNA (paired by sample)
# -----------------------------------------------------------------------------
vars <- c("S", "PD", "SES_PD", "MPD", "NRI", "MNTD", "NTI")
paired_long <- alpha %>%
  filter(Method != "Combined") %>%
  select(Method, Locality, all_of(vars)) %>%
  pivot_longer(all_of(vars), names_to = "Metric")

tests <- bind_rows(lapply(vars, function(v) {
  w <- paired_long %>% filter(Metric == v) %>%
    pivot_wider(names_from = Method, values_from = value) %>%
    filter(!is.na(Traditional), !is.na(eDNA))
  tt <- suppressWarnings(wilcox.test(w$Traditional, w$eDNA, paired = TRUE, exact = FALSE))
  data.frame(Metric = v, n_pairs = nrow(w),
             Median_Trad = median(w$Traditional), Median_eDNA = median(w$eDNA),
             Mean_Trad = mean(w$Traditional), Mean_eDNA = mean(w$eDNA),
             V = unname(tt$statistic), p = tt$p.value)
}))
tests$p_holm <- p.adjust(tests$p, "holm")

# Proportion of communities with significant structure (two-tailed, alpha = 0.05)
signif_structure <- alpha %>%
  group_by(Method) %>%
  summarise(
    n = n(),
    PD_lower = sum(p_PD <= 0.025), PD_higher = sum(p_PD >= 0.975),
    MPD_clustered = sum(p_MPD <= 0.025, na.rm = TRUE), MPD_overdispersed = sum(p_MPD >= 0.975, na.rm = TRUE),
    MNTD_clustered = sum(p_MNTD <= 0.025, na.rm = TRUE), MNTD_overdispersed = sum(p_MNTD >= 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# 4. Gamma: total, shared and exclusive PD per method
# -----------------------------------------------------------------------------
# PD including the path to the root of the pool tree, as for the assemblages (Faith, 1992)
pd_set <- function(tips) pd(matrix(1L, 1, length(tips), dimnames = list("a", tips)), pool_tree, include.root = TRUE)$PD
sp_trad <- unique(pa$tip[pa$Sampling_Method == "Traditional"])
sp_edna <- unique(pa$tip[pa$Sampling_Method == "eDNA"])
sp_all  <- union(sp_trad, sp_edna)
pd_all  <- pd_set(sp_all)
gamma <- data.frame(
  Set = c("Traditional", "eDNA", "Combined"),
  S = c(length(sp_trad), length(sp_edna), length(sp_all)),
  PD = c(pd_set(sp_trad), pd_set(sp_edna), pd_all)
)
gamma$PD_pct_of_total <- 100 * gamma$PD / pd_all
# Exclusive PD = total PD - PD of the other method
partition <- data.frame(
  Component = c("Shared", "Traditional only", "eDNA only"),
  PD = c(pd_set(sp_trad) + pd_set(sp_edna) - pd_all, pd_all - pd_set(sp_edna), pd_all - pd_set(sp_trad)),
  S = c(length(intersect(sp_trad, sp_edna)), length(setdiff(sp_trad, sp_edna)), length(setdiff(sp_edna, sp_trad)))
)
partition$PD_pct <- 100 * partition$PD / pd_all

# Gamma per stream (4 campaigns pooled)
gamma_stream <- pa_all %>%
  mutate(Stream = sub("_.*$", "", Locality)) %>%
  distinct(Sampling_Method, Stream, tip) %>%
  group_by(Sampling_Method, Stream) %>%
  summarise(S = n(), PD = pd_set(tip), .groups = "drop")

# -----------------------------------------------------------------------------
# 5. Export
# -----------------------------------------------------------------------------
write.csv(alpha, file.path(out_res, "communities_alpha_pd_ses.csv"), row.names = FALSE)
write.csv(tests, file.path(out_res, "communities_trad_vs_edna_wilcoxon.csv"), row.names = FALSE)
write.csv(signif_structure, file.path(out_res, "communities_significant_structure.csv"), row.names = FALSE)
write.csv(gamma, file.path(out_res, "gamma_pd_by_method.csv"), row.names = FALSE)
write.csv(partition, file.path(out_res, "gamma_pd_partition.csv"), row.names = FALSE)
write.csv(gamma_stream, file.path(out_res, "gamma_pd_by_stream.csv"), row.names = FALSE)
write.csv(comm, file.path(out_res, "community_matrix_pa.csv"))

wb <- createWorkbook()
for (nm in c("alpha", "tests", "signif_structure", "gamma", "partition", "gamma_stream")) {
  addWorksheet(wb, nm)
  writeData(wb, nm, get(nm))
}
saveWorkbook(wb, file.path(out_res, "community_phylogenetic_diversity.xlsx"), overwrite = TRUE)

# -----------------------------------------------------------------------------
# 6. Figures
# -----------------------------------------------------------------------------
method_cols <- c(Traditional = "#d1495b", eDNA = "#2a78b5", Combined = "grey45")
fig_theme <- theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank())

save_fig <- function(p, name, w, h) {
  ggsave(file.path(out_fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(out_fig, paste0(name, ".pdf")), p, width = w, height = h)
}

# Fig 6 - Traditional vs eDNA, paired
metric_labels <- c(S = "Richness (S)", PD = "PD (Ma)", SES_PD = "SES.PD",
                   NRI = "NRI (-SES.MPD)", NTI = "NTI (-SES.MNTD)")
lab_p <- tests %>% filter(Metric %in% names(metric_labels)) %>%
  mutate(Metric = factor(metric_labels[Metric], levels = metric_labels),
         lab = ifelse(p_holm < 0.001, "p < 0.001", sprintf("p = %.3f", p_holm)))
d6 <- paired_long %>% filter(Metric %in% names(metric_labels)) %>%
  mutate(Metric = factor(metric_labels[Metric], levels = metric_labels))
f6 <- ggplot(d6, aes(Method, value)) +
  geom_hline(data = data.frame(Metric = factor(metric_labels[c("SES_PD", "NRI", "NTI")], levels = metric_labels)),
             aes(yintercept = 0), linetype = 2, colour = "grey60") +
  geom_line(aes(group = Locality), colour = "grey80", linewidth = 0.3) +
  geom_boxplot(aes(fill = Method), width = 0.45, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.06, size = 0.8, alpha = 0.6) +
  geom_text(data = lab_p, aes(x = 1.5, y = Inf, label = lab), vjust = 1.3, size = 3) +
  facet_wrap(~Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = method_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = NULL,
       title = "Phylogenetic diversity per sample: Traditional vs eDNA (48 stream-campaign pairs)",
       caption = "Paired Wilcoxon, Holm-adjusted p. Grey lines link the same sample.") +
  fig_theme
save_fig(f6, "fig6_trad_vs_edna_paired", 12, 4)

# Fig 7 - PD vs richness with null expectation
f7 <- ggplot(alpha %>% filter(Method != "Combined"), aes(S, PD)) +
  geom_line(data = alpha %>% group_by(S) %>% summarise(PD_null_mean = mean(PD_null_mean)),
            aes(S, PD_null_mean), colour = "grey40", linetype = 2) +
  geom_point(aes(colour = Method), size = 2, alpha = 0.8) +
  scale_colour_manual(values = method_cols) +
  labs(x = "Richness (S)", y = "PD (Ma)", colour = NULL,
       title = "Observed PD vs richness",
       subtitle = "Dashed line: null mean (random draw of S taxa from the 119-taxa pool)") +
  fig_theme + theme(legend.position = "top")

# Fig 8 - gamma: PD partition between methods
f8 <- ggplot(partition %>% mutate(Component = factor(Component, levels = rev(Component))),
             aes(PD_pct, Component, fill = Component)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%% (%.0f Ma; %d taxa)", PD_pct, PD, S)), hjust = -0.05, size = 3) +
  scale_fill_manual(values = c("Shared" = "grey55", "Traditional only" = "#d1495b",
                               "eDNA only" = "#2a78b5"), guide = "none") +
  scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(x = sprintf("%% of the total PD detected (%.0f Ma, %d taxa)", pd_all, length(sp_all)), y = NULL,
       title = "Gamma PD: shared and exclusive by method") +
  fig_theme + theme(panel.grid.major.y = element_blank())
save_fig(f7 + f8 + plot_layout(widths = c(1, 1.2)), "fig7_pd_richness_gamma", 12, 4.2)

# Fig 8 (file) - PD per stream and campaign
d9 <- alpha %>% filter(Method != "Combined") %>%
  mutate(Campaign = factor(paste(substr(Month, 1, 3), Year), levels = c("Aug 2022", "Feb 2023", "Aug 2023", "Feb 2024")))
f9 <- ggplot(d9, aes(Campaign, PD, colour = Method, group = Method)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  facet_wrap(~Stream, ncol = 6) +
  scale_colour_manual(values = method_cols) +
  labs(x = NULL, y = "PD (Ma)", colour = NULL, title = "PD per stream across campaigns") +
  fig_theme + theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
save_fig(f9, "fig8_pd_stream_campaign", 11, 5)

# -----------------------------------------------------------------------------
cat("\n== Paired tests ==\n"); print(tests)
cat("\n== Significant structure ==\n"); print(signif_structure)
cat("\n== Gamma ==\n"); print(gamma); print(partition)
cat("\n== Alpha summary by method ==\n")
print(alpha %>% group_by(Method) %>% summarise(across(c(S, PD, SES_PD, NRI, NTI), ~round(mean(.x, na.rm = TRUE), 2))))
