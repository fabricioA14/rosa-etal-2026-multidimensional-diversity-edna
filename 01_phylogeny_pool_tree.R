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
# Script:  01_phylogeny_pool_tree.R
# Purpose: Phylogenetic diversity of the regional species pool (119 taxa):
#          pruned time-calibrated tree of Cassemiro et al. (2023) and
#          correspondence between taxa and tips (Table S1).
# Input:   pool_taxa_traits.xlsx (taxon list and tree tips)
#          South-American-freshwater-fish-Fishes.zip (tree
#          3167_Neotrop_Freshwater_Fish_spp_time_tree.tre; Cassemiro et al.
#          2023, PNAS 120: e2211974120)
# Output:  phylogeny/{data,results,figures}
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(picante)
  library(openxlsx)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

set.seed(2024)
n_runs <- 999

pool_file <- "pool_taxa_traits.xlsx"
tree_zip  <- "South-American-freshwater-fish-Fishes.zip"
out_data  <- "phylogeny/data"
out_res   <- "phylogeny/results"
out_fig   <- "phylogeny/figures"
for (d in c(out_data, out_res, out_fig)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------
tmp <- file.path(tempdir(), "repo")
unzip(tree_zip, exdir = tmp)
f_tree <- list.files(tmp, pattern = "time_tree\\.tre$", recursive = TRUE, full.names = TRUE)

tree <- read.nexus(f_tree)
file.copy(f_tree, file.path(out_data, basename(f_tree)), overwrite = TRUE)

pool <- read.xlsx(pool_file, sheet = 1)
pool$Tree_tip <- trimws(pool$Tree_tip)
pool$Family <- sub("_.*$", "", pool$Tree_tip)

orders <- c(
  Potamotrygonidae = "Myliobatiformes", Belonidae = "Beloniformes",
  Rivulidae = "Cyprinodontiformes", Cichlidae = "Cichliformes",
  Synbranchidae = "Synbranchiformes", Tetraodontidae = "Tetraodontiformes",
  Acestrorhynchidae = "Characiformes", Anostomidae = "Characiformes",
  Characidae = "Characiformes", Chilodontidae = "Characiformes",
  Crenuchidae = "Characiformes", Curimatidae = "Characiformes",
  Erythrinidae = "Characiformes", Hemiodontidae = "Characiformes",
  Iguanodectidae = "Characiformes", Lebiasinidae = "Characiformes",
  Serrasalmidae = "Characiformes",
  Apteronotidae = "Gymnotiformes", Gymnotidae = "Gymnotiformes",
  Hypopomidae = "Gymnotiformes", Rhamphichthyidae = "Gymnotiformes",
  Sternopygidae = "Gymnotiformes",
  Aspredinidae = "Siluriformes", Auchenipteridae = "Siluriformes",
  Callichthyidae = "Siluriformes", Cetopsidae = "Siluriformes",
  Doradidae = "Siluriformes", Heptapteridae = "Siluriformes",
  Loricariidae = "Siluriformes", Pimelodidae = "Siluriformes",
  Pseudopimelodidae = "Siluriformes", Trichomycteridae = "Siluriformes"
)
stopifnot(all(pool$Family %in% names(orders)))
pool$Order <- unname(orders[pool$Family])

# -----------------------------------------------------------------------------
# 2. Name-tip matching and pruning
# -----------------------------------------------------------------------------
match_tab <- data.frame(
  Taxon = pool$Taxon,
  Tree_tip = pool$Tree_tip,
  In_tree = pool$Tree_tip %in% tree$tip.label,
  Substitute = sub("^[^_]+_", "", pool$Tree_tip) != gsub(" ", "_", trimws(pool$Taxon))
)
stopifnot(all(match_tab$In_tree))

pool_tree <- keep.tip(tree, pool$Tree_tip)
write.tree(pool_tree, file.path(out_data, "pool_tree_119taxa.nwk"))
write.nexus(pool_tree, file = file.path(out_data, "pool_tree_119taxa.nex"))

# -----------------------------------------------------------------------------
# 3. Reference tree for the null model: full original tree (3,167 tips)
#    The pool is linked to the tree by the names in column Tree_tip.
# -----------------------------------------------------------------------------
ref_tree <- tree
# picante fails with a one-row matrix; duplicate the row and keep only the first
comm <- matrix(as.integer(ref_tree$tip.label %in% pool$Tree_tip), nrow = 2,
               ncol = Ntip(ref_tree), byrow = TRUE,
               dimnames = list(c("Paragominas", "dup"), ref_tree$tip.label))
dist_ref <- cophenetic(ref_tree)

# -----------------------------------------------------------------------------
# 4. Pool metrics
# -----------------------------------------------------------------------------
pool_metrics <- function(tips, tr = tree) {
  sub_t <- keep.tip(tr, tips)
  d <- cophenetic(sub_t)
  diag(d) <- NA
  data.frame(
    S = length(tips),
    PD = sum(sub_t$edge.length),
    PD_with_root = sum(sub_t$edge.length) +
      (max(node.depth.edgelength(tr)) - max(node.depth.edgelength(sub_t))),
    MRCA_age = max(node.depth.edgelength(sub_t)),
    MPD = mean(d, na.rm = TRUE),
    MNTD = mean(apply(d, 1, min, na.rm = TRUE))
  )
}

stingray <- "Potamotrygonidae_Potamotrygon_orbignyi"
tilapia <- "Cichlidae_Retroculus_lapidifer"   # Oreochromis niloticus (non-native)
scenarios <- list(
  "Full pool" = pool$Tree_tip,
  "Without Potamotrygon" = setdiff(pool$Tree_tip, stingray),
  "Without Oreochromis" = setdiff(pool$Tree_tip, tilapia),
  "Without both" = setdiff(pool$Tree_tip, c(stingray, tilapia))
)
sensitivity <- bind_rows(lapply(scenarios, pool_metrics), .id = "Scenario")

# -----------------------------------------------------------------------------
# 5. Null model (taxa.labels: random draw of S tips from the full tree)
# -----------------------------------------------------------------------------
ses_pd   <- ses.pd(comm, ref_tree, null.model = "taxa.labels", runs = n_runs, include.root = FALSE)
ses_mpd  <- ses.mpd(comm, dist_ref, null.model = "taxa.labels", runs = n_runs)
ses_mntd <- ses.mntd(comm, dist_ref, null.model = "taxa.labels", runs = n_runs)
ses_pd <- ses_pd[1, ]; ses_mpd <- ses_mpd[1, ]; ses_mntd <- ses_mntd[1, ]

# Keep the null distributions for the figure (same model, explicit draws)
null_draws <- t(replicate(n_runs, {
  s <- sample(ref_tree$tip.label, nrow(pool))
  dd <- dist_ref[s, s]
  diag(dd) <- NA
  c(PD = sum(keep.tip(ref_tree, s)$edge.length),
    MPD = mean(dd, na.rm = TRUE),
    MNTD = mean(apply(dd, 1, min, na.rm = TRUE)))
}))

ses_tab <- data.frame(
  Metric = c("PD", "MPD", "MNTD"),
  Observed = c(ses_pd$pd.obs, ses_mpd$mpd.obs, ses_mntd$mntd.obs),
  Null_mean = c(ses_pd$pd.rand.mean, ses_mpd$mpd.rand.mean, ses_mntd$mntd.rand.mean),
  Null_SD = c(ses_pd$pd.rand.sd, ses_mpd$mpd.rand.sd, ses_mntd$mntd.rand.sd),
  SES = c(ses_pd$pd.obs.z, ses_mpd$mpd.obs.z, ses_mntd$mntd.obs.z),
  p_rank = c(ses_pd$pd.obs.p, ses_mpd$mpd.obs.p, ses_mntd$mntd.obs.p),
  Runs = n_runs,
  Ref_tree_S = Ntip(ref_tree)
)
# NRI and NTI = -SES
ses_tab$NRI_NTI <- c(NA, -ses_tab$SES[2], -ses_tab$SES[3])

# -----------------------------------------------------------------------------
# 6. Evolutionary distinctiveness (ED, fair proportion) and contribution per taxon
# -----------------------------------------------------------------------------
ed <- evol.distinct(pool_tree, type = "fair.proportion")
names(ed) <- c("Tree_tip", "ED")
pd_total <- sum(pool_tree$edge.length)

# PD lost when each species is removed (length of its unique terminal branch)
terminal <- setNames(pool_tree$edge.length[match(seq_len(Ntip(pool_tree)), pool_tree$edge[, 2])],
                     pool_tree$tip.label)

by_species <- pool %>%
  select(Taxon, Tree_tip, Family, Order) %>%
  left_join(ed, by = "Tree_tip") %>%
  mutate(Terminal_branch = terminal[Tree_tip],
         ED_pct_PD = 100 * ED / pd_total) %>%
  arrange(desc(ED))

by_order <- by_species %>%
  group_by(Order) %>%
  summarise(S = n(), ED_sum = sum(ED), .groups = "drop") %>%
  mutate(pct_S = 100 * S / sum(S), pct_PD = 100 * ED_sum / sum(ED_sum)) %>%
  arrange(desc(ED_sum))

by_family <- by_species %>%
  group_by(Order, Family) %>%
  summarise(S = n(), ED_sum = sum(ED), ED_mean = mean(ED), .groups = "drop") %>%
  mutate(pct_PD = 100 * ED_sum / sum(ED_sum)) %>%
  arrange(desc(ED_sum))

# PD accumulation with richness (random subsamples of the local pool)
accumulation <- bind_rows(lapply(c(2, seq(5, nrow(pool), by = 5), nrow(pool)), function(k) {
  v <- replicate(199, sum(keep.tip(pool_tree, sample(pool_tree$tip.label, k))$edge.length))
  data.frame(S = k, PD_mean = mean(v), lo = quantile(v, 0.025), hi = quantile(v, 0.975))
}))

# -----------------------------------------------------------------------------
# 7. Export tables
# -----------------------------------------------------------------------------
pool_summary <- pool_metrics(pool$Tree_tip)
write.csv(pool_summary, file.path(out_res, "pd_pool_summary.csv"), row.names = FALSE)
write.csv(sensitivity, file.path(out_res, "pd_sensitivity.csv"), row.names = FALSE)
write.csv(ses_tab, file.path(out_res, "ses_null_model.csv"), row.names = FALSE)
write.csv(by_species, file.path(out_res, "ed_by_species.csv"), row.names = FALSE)
write.csv(by_order, file.path(out_res, "pd_by_order.csv"), row.names = FALSE)
write.csv(by_family, file.path(out_res, "pd_by_family.csv"), row.names = FALSE)
write.csv(match_tab, file.path(out_res, "taxon_tip_match.csv"), row.names = FALSE)

wb <- createWorkbook()
for (nm in c("pool_summary", "sensitivity", "ses_tab", "by_species", "by_order", "by_family", "match_tab")) {
  addWorksheet(wb, nm)
  writeData(wb, nm, get(nm))
}
saveWorkbook(wb, file.path(out_res, "pool_phylogenetic_diversity.xlsx"), overwrite = TRUE)

# -----------------------------------------------------------------------------
# 8. Figures
# -----------------------------------------------------------------------------
order_cols <- c(
  Characiformes = "#2a78b5", Siluriformes = "#d1495b", Gymnotiformes = "#edae49",
  Cichliformes = "#3e9e5e", Myliobatiformes = "#6b4c9a", Beloniformes = "#00a3a3",
  Cyprinodontiformes = "#e07b39", Synbranchiformes = "#8c6d46",
  Tetraodontiformes = "#c45ca2"
)
order_levels <- by_order$Order

save_fig <- function(p, name, w, h) {
  ggsave(file.path(out_fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(out_fig, paste0(name, ".pdf")), p, width = w, height = h)
}

# Fig 1 - phylogeny of the pool
plot_tree <- ladderize(pool_tree)
tip_order <- setNames(pool$Order, pool$Tree_tip)[plot_tree$tip.label]
tip_taxon <- setNames(pool$Taxon, pool$Tree_tip)[plot_tree$tip.label]
tip_label <- ifelse(match_tab$Substitute[match(plot_tree$tip.label, match_tab$Tree_tip)],
                    paste0(tip_taxon, "*"), tip_taxon)

# branch colour: a branch takes the order colour when all descendant tips belong to that order
n_tip <- Ntip(plot_tree)
desc_orders <- vector("list", n_tip + plot_tree$Nnode)
for (i in seq_len(n_tip)) desc_orders[[i]] <- tip_order[i]
po <- postorder(plot_tree)
for (e in po) {
  p <- plot_tree$edge[e, 1]; ch <- plot_tree$edge[e, 2]
  desc_orders[[p]] <- unique(c(desc_orders[[p]], desc_orders[[ch]]))
}
edge_col <- sapply(plot_tree$edge[, 2], function(n) {
  o <- desc_orders[[n]]
  if (length(o) == 1) order_cols[o] else "grey55"
})

labelled_tree <- plot_tree
labelled_tree$tip.label <- tip_label

draw_tree <- function() {
  par(mar = c(1, 1, 1, 1))
  plot(labelled_tree, type = "fan", show.tip.label = TRUE, label.offset = 4,
       tip.color = order_cols[tip_order], edge.color = edge_col, edge.width = 1.4,
       cex = 0.42, no.margin = TRUE, open.angle = 8)
  legend("bottomleft", legend = order_levels, col = order_cols[order_levels], lwd = 3,
         bty = "n", cex = 0.65, title = "Order", title.adj = 0)
  legend("topleft", bty = "n", cex = 0.62,
         legend = c(sprintf("S = %d taxa", nrow(pool)),
                    sprintf("PD = %.1f Ma", pd_total),
                    "* substitute tip in the tree"))
}
png(file.path(out_fig, "fig1_pool_phylogeny.png"), width = 9, height = 9, units = "in", res = 300)
draw_tree(); dev.off()
pdf(file.path(out_fig, "fig1_pool_phylogeny.pdf"), width = 9, height = 9)
draw_tree(); dev.off()

# Fig 2 - ED per species (top 30)
top_ed <- by_species %>% slice_max(ED, n = 30) %>%
  mutate(Taxon = factor(Taxon, levels = rev(Taxon)), Order = factor(Order, levels = order_levels))
f2 <- ggplot(top_ed, aes(ED, Taxon, colour = Order)) +
  geom_segment(aes(x = 0, xend = ED, yend = Taxon), linewidth = 0.8) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = order_cols, drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Evolutionary distinctiveness - ED (Ma, fair proportion)", y = NULL, colour = "Order",
       title = "30 most evolutionarily distinct taxa in the pool") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(face = "italic"), panel.grid.major.y = element_blank(),
        legend.position = c(0.8, 0.3), legend.background = element_rect(fill = "white", colour = NA))
save_fig(f2, "fig2_ed_top30", 7, 7)

# Fig 3 - contribution of each order to richness and PD
order_long <- by_order %>%
  select(Order, `Richness (S)` = pct_S, `PD (sum of ED)` = pct_PD) %>%
  pivot_longer(-Order, names_to = "Component", values_to = "pct") %>%
  mutate(Order = factor(Order, levels = rev(order_levels)),
         Component = factor(Component, levels = c("Richness (S)", "PD (sum of ED)")))
f3 <- ggplot(order_long, aes(pct, Order, fill = Component)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_dodge(width = 0.75),
            hjust = -0.15, size = 2.8) +
  scale_fill_manual(values = c("Richness (S)" = "grey65", "PD (sum of ED)" = "#2a78b5")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "% of the pool total", y = NULL, fill = NULL,
       title = "Contribution of each order to richness and phylogenetic diversity") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), legend.position = "top")
save_fig(f3, "fig3_order_contribution", 7, 5)

# Fig 4 - null model
null_long <- as.data.frame(null_draws) %>% pivot_longer(everything(), names_to = "Metric")
obs <- ses_tab %>% mutate(lab = sprintf("obs = %.1f\nSES = %.2f\np = %.3f", Observed, SES, p_rank))
null_long$Metric <- factor(null_long$Metric, levels = c("PD", "MPD", "MNTD"))
obs$Metric <- factor(obs$Metric, levels = c("PD", "MPD", "MNTD"))
f4 <- ggplot(null_long, aes(value)) +
  geom_histogram(bins = 40, fill = "grey75", colour = "white") +
  geom_vline(data = obs, aes(xintercept = Observed), colour = "#d1495b", linewidth = 1) +
  geom_label(data = obs, aes(x = Inf, y = Inf, label = lab), hjust = 1.05, vjust = 1.1,
             size = 2.8, label.size = 0) +
  facet_wrap(~Metric, scales = "free") +
  labs(x = "Ma", y = "Frequency",
       title = sprintf("Observed pool vs. %d random draws of %d taxa from the full tree (%d taxa)",
                       n_runs, nrow(pool), Ntip(ref_tree))) +
  theme_minimal(base_size = 10)
save_fig(f4, "fig4_null_model", 10, 3.8)

# Fig 5 - PD accumulation and sensitivity
f5a <- ggplot(accumulation, aes(S, PD_mean)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#2a78b5", alpha = 0.2) +
  geom_line(colour = "#2a78b5", linewidth = 1) +
  geom_point(data = pool_summary, aes(S, PD), colour = "#d1495b", size = 2.5) +
  labs(x = "Number of taxa", y = "PD (Ma)", title = "PD accumulation in the pool") +
  theme_minimal(base_size = 10)
f5b <- ggplot(sensitivity %>% mutate(Scenario = factor(Scenario, levels = rev(names(scenarios)))),
              aes(PD, Scenario)) +
  geom_col(fill = "#2a78b5", width = 0.6) +
  geom_text(aes(label = sprintf("%.0f Ma (S = %d)", PD, S)), hjust = 1.05, colour = "white", size = 3) +
  labs(x = "PD (Ma)", y = NULL, title = "PD sensitivity") +
  theme_minimal(base_size = 10) + theme(panel.grid.major.y = element_blank())
save_fig(f5a + f5b + plot_layout(widths = c(1.2, 1)), "fig5_accumulation_sensitivity", 10, 4)

# -----------------------------------------------------------------------------
cat("\n== Summary ==\n"); print(pool_summary)
cat("\n== Sensitivity ==\n"); print(sensitivity)
cat("\n== Null model ==\n"); print(ses_tab)
cat("\n== By order ==\n"); print(by_order)
cat("\n== Top 10 ED ==\n"); print(head(by_species, 10))
