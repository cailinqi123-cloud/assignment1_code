##***************************
## Assignment 1, Trichinella nativa comparing with Trichinella Spiralis
##
## Jordyn Anene  
## Luis Humberto Villasana Avila  
## Shariq Mumtaz Hashmi
## Sumayyah Khan
##
## 2026-11-02
## Assignment #2 secondary contributor: Linqi Cai (2026-11-23)
## Edit 1: locate input files reliably, validate inputs, and create output directory.
## Edit 2: consistently classify TE annotations in both species and flag unknown classes.
## Edit 3: make plot captions data-driven and PCA robust to constant features.
## Original authors' numerical interpretations below should be rechecked after rerunning.
##
##***************************

## Cleaning the environment to make sure that the code can run on its own ----
rm(list = ls())

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
library(Biostrings)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())


## Preparing CSV files ----
# EDIT 1 (reproducibility): work from the repository root, assignment1/, or assignment1/R/.
# Paths are checked explicitly instead of assuming a particular working directory.
find_data_file <- function(filename) {
  candidates <- file.path(c("assignment1/data", "data", "../data"), filename)
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop("Cannot find ", filename, ". Working directory: ", getwd(),
         ". Open the complete GitHub repository (including assignment1/data).")
  }
  normalizePath(existing[1], mustWork = TRUE)
}

nat_file <- find_data_file("GCA_001447565.2.csv")
spiralis_file <- find_data_file("GCA_001447595.1.csv")

# Save figures beside the data folder, regardless of where the R session started.
output_dir <- file.path(dirname(dirname(nat_file)), "Output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read inputs using the same directory capitalization on all operating systems.
df_te_nat <- read.csv(nat_file, stringsAsFactors = FALSE)
df_trich_te <- read.csv(spiralis_file, stringsAsFactors = FALSE)

required_columns <- c("qry_id", "matching_class", "qry_left", "in_repeat_end",
                      "perc_div", "perc_del", "qry_width")
for (species_input in list(T_nativa = df_te_nat, T_spiralis = df_trich_te)) {
  missing_columns <- base::setdiff(required_columns, names(species_input))
  if (length(missing_columns) > 0L) {
    stop("Input CSV is missing required column(s): ",
         paste(missing_columns, collapse = ", "))
  }
}


#fixing the issue with qry_left and in_repeat_end columns for Nativa:
df_te_nat = df_te_nat |>
  mutate(qry_left2 = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")), .after = qry_left) |>
  mutate(in_repeat_end2 = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")), .after = in_repeat_end)

#fixing the issue with qry_left and in_repeat_end columns for Spiralis:
df_trich_te = df_trich_te |>
  mutate(qry_left2 = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")), .after = qry_left) |>
  mutate(in_repeat_end2 = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")), .after = in_repeat_end)


## Removing non-TE elements and classifying TE annotations ----
# EDIT 2 (correctness): the original used %notin% (undefined in this script)
# and applied different prefix rules across species. Classification should be
# based on each annotated TE row, rather than assuming all remaining rows are retrotransposons.
non_te_classes <- c("Simple_repeat", "simple_repeat", "Low_complexity",
                    "Satellite", "rRNA", "snRNA", "ARTEFACT", "Unknown")

# Keep original copies for the repeat-density calculation below.
df_te_nat_org <- df_te_nat
df_trich_te_org <- df_trich_te

classify_te <- function(data, species_name) {
  classified <- data |>
    filter(!is.na(matching_class), !(matching_class %in% non_te_classes)) |>
    mutate(TE_type = case_when(
      str_detect(matching_class, "^(DNA/|RC)") ~ "DNA Transposons",
      str_detect(matching_class, "^(PLE|LTR|LINE|SINE)") ~ "Retrotransposons",
      TRUE ~ NA_character_
    ))

  # Report unexpected labels instead of silently treating them as retrotransposons.
  unknown <- classified |>
    filter(is.na(TE_type)) |>
    count(matching_class, sort = TRUE)
  if (nrow(unknown) > 0L) {
    warning(species_name, ": excluded TE classes that need manual review: ",
            paste(unknown$matching_class, collapse = ", "))
  }
  classified |>
    filter(!is.na(TE_type))
}

df_te_nat <- classify_te(df_te_nat_org, "T. nativa")
df_trich_te <- classify_te(df_trich_te_org, "T. spiralis")

# Quality checks: the two cleaned datasets contain only the intended TE groups.
stopifnot(!any(df_te_nat$matching_class %in% non_te_classes),
          !any(df_trich_te$matching_class %in% non_te_classes),
          all(df_te_nat$TE_type %in% c("DNA Transposons", "Retrotransposons")),
          all(df_trich_te$TE_type %in% c("DNA Transposons", "Retrotransposons")))

# Reproducible summaries of the classifications in each species.
df_te_nat |> count(TE_type)
df_trich_te |> count(TE_type)

# To observe what are the TE elements present in Nativa
df_te_nat$matching_class |> unique()

# To observe what are the TE elements present in Spiralis
df_trich_te$matching_class |> unique()

# To identify how many unique qry_ids we have in Nativa
df_te_nat$qry_id |> unique() |> length()
#=============> we have 438 unique contigs

# To identify how many unique qry_ids we have in Spiralis
df_trich_te$qry_id |> unique() |> length()
#=============> we have 315 unique contigs 

# To obtain a list of the dna transposons for T. nativa 
ls_dna_qry = df_te_nat |> 
  filter(TE_type == "DNA Transposons") |>
  pull(qry_id) |> unique()
ls_dna_qry

# To obtain a list of the dna transposons for T. spiralis 
ls_ts_dna_qry = df_trich_te |>
  filter(TE_type == "DNA Transposons") |>
  pull(qry_id) |> unique()
ls_ts_dna_qry

# To obtain a list of a all the retrotransposons in Nativa
ls_retro_qry = df_te_nat |>
  filter(TE_type == "Retrotransposons") |>
  pull(qry_id) |> unique()
ls_retro_qry

# To obtain a list of a all the retrotransposons in Spiralis
ls_ts_retro_qry = df_trich_te |>
  filter(TE_type == "Retrotransposons") |>
  pull(qry_id) |> unique()
ls_ts_retro_qry

# Checking that the number of DNA transposons and retrotransposons match in Nativa
Biostrings::setdiff(ls_dna_qry, ls_retro_qry) |> length()
Biostrings::setdiff(ls_retro_qry, ls_dna_qry) |> length()

# Checking that the number of DNA transposons and retrotransposons match in Spiralis
Biostrings::setdiff(ls_ts_dna_qry, ls_ts_retro_qry) |> length()
Biostrings::setdiff(ls_ts_retro_qry, ls_ts_dna_qry) |> length()


## Creating different data frames according to the type of TE ----
#Data frame with only DNA transponsos in Nativa
df_dna_only = df_te_nat |>
  filter(qry_id %in% Biostrings::setdiff(ls_dna_qry, ls_retro_qry))
df_dna_only

#Data frame with only DNA transponsos in Spiralis
df_ts_dnaonly = df_trich_te |>
  filter(qry_id %in% Biostrings::setdiff(ls_ts_dna_qry, ls_ts_retro_qry))
df_ts_dnaonly

#Data frame with only retrotransponsos in Nativa
df_retro_only = df_te_nat |>
  filter(qry_id %in% Biostrings::setdiff(ls_retro_qry, ls_dna_qry))
df_retro_only

#Data frame with only retrotransponsos in spiralis
df_ts_retroonly = df_trich_te |>
  filter(qry_id %in% Biostrings::setdiff(ls_ts_retro_qry, ls_ts_dna_qry))
df_ts_retroonly

# Data frame with both DNA transponsos and retrotransponsos in Nativa
df_bothTE_nat = df_te_nat |>
  filter(qry_id %in% Biostrings::intersect(ls_dna_qry, ls_retro_qry))
df_bothTE_nat

# Data frame with both DNA transponsos and retrotransponsos in Spiralis
df_ts_both = df_trich_te |>
  filter(qry_id %in% Biostrings::intersect(ls_ts_dna_qry, ls_ts_retro_qry))
df_ts_both


## Finding the inflection point for repeat-dense contigs ----
# EDIT 1 (reproducibility): save diagnostic plots as image files rather than
# drawing in RStudio's interactive Plots pane, which can be too small on screen.
# Each species is plotted on linear and log10 scales using the original rows.
save_contig_rank_plots <- function(df, species_label, file_suffix) {
  ranked_repeats <- df |>
    count(qry_id, name = "repeat_count") |>
    arrange(desc(repeat_count)) |>
    mutate(contig_rank = row_number())

  base_plot <- ggplot(ranked_repeats, aes(x = contig_rank, y = repeat_count)) +
    geom_point(size = 0.7, alpha = 0.65) +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "grey40") +
    labs(x = "Contig rank (ordered by repeat count)",
         y = "Repeat annotations per contig",
         title = paste(species_label, "repeat count by contig"),
         caption = "Dashed line: 50-repeat exploratory cutoff") +
    theme_light()

  ggsave(file.path(output_dir, paste0("repeat_count_linear_", file_suffix, ".png")),
         plot = base_plot, width = 9, height = 5)
  ggsave(file.path(output_dir, paste0("repeat_count_log10_", file_suffix, ".png")),
         plot = base_plot + scale_y_log10() +
           labs(y = "Repeat annotations per contig (log10 scale)"),
         width = 9, height = 5)
  invisible(ranked_repeats)
}

save_contig_rank_plots(df_trich_te_org, "T. spiralis", "spiralis")
save_contig_rank_plots(df_te_nat_org, "T. nativa", "nativa")
# The original script uses >50 repeat annotations (not sequence length).


#================> Cutoff for repeat-dense contigs: >50 repeats per contig (both species).
#================> the curve bends from a gentle to a steep decline at ~50
#================> Repeats per contig (all repeat classes >>>> org tables) were sorted and plotted on a log10 scale;

# To obtain from the original cvs file the contigs larger than 50 in nativa
ls_dense_nat = df_te_nat_org |>
  group_by(qry_id) |>
  summarise(count = n()) |>
  arrange(desc(count)) |>
  filter(count > 50) |>
  pull(qry_id)
ls_dense_nat

# To obtain number of large contigs present in nativa
ls_dense_nat |> length()
#=================> In T.nativa there are 268 contigs

# To obtain from the original cvs file the contigs larger than 50 in Spiralis 
ls_dense_ts = df_trich_te_org |>
  group_by(qry_id) |>
  summarise(count = n()) |>
  arrange(desc(count)) |>
  filter(count > 50) |>
  pull(qry_id)

# To obtain number of large contigs present in Spiralis
ls_dense_ts |> length()
#=================> In T.spiralis there are 220 contigs 

#=================>  The same cutoff is used for both species so the comparison is consistent.
#=================>  Note: we haven't used FASTA files, contigs are ranked by repeat count, not by length, so "repeat dense" rather than "largest".


# Creating a new table with TE rows that are in the ls_dense for which species, so we can plot them against perc_del and perc_div

df_dense_nat = df_te_nat |>
  filter(qry_id %in% ls_dense_nat)

df_dense_ts = df_trich_te |>
  filter(qry_id %in% ls_dense_ts)

nrow(df_te_nat)
nrow(df_dense_nat)

nrow(df_trich_te)
nrow(df_dense_ts)

#====> Dense contigs hold 823 of 1319 TEs in T.Nativa and 838 of 1037 TEs in T.Spiralis.
  
## Creating new columns classifying contigs as either DNA Transposons or Retrotransposons and Nativa or spiralis to later create a single table to make comparisons and plotting easier ----

# EDIT 2: TE_type is inherited from the validated classification above.
# Do not reclassify each selected subset or silently assign unknown classes to Retrotransposons.
stopifnot(all(df_dense_nat$TE_type %in% c("DNA Transposons", "Retrotransposons")),
          all(df_dense_ts$TE_type %in% c("DNA Transposons", "Retrotransposons")))
  
# Identifying the amount of DNA Transposons and Retrotransposons in both Nativa and Spiralis to check values
df_dense_nat |>
  count(TE_type)

df_dense_ts |>
  count(TE_type)

df_dense_nat$qry_id |> unique()
df_dense_ts$qry_id |> unique()
  
df_dense_nat |>
  filter(TE_type == "Retrotransposons") |>
  count(matching_class)  

df_dense_ts |>
  filter(TE_type == "Retrotransposons") |>
  count(matching_class)

#creating a new table with both species labelled in the new column
df_dense_nat = df_dense_nat |>
  mutate(Species = "T.nativa")
  
df_dense_ts = df_dense_ts |>
  mutate(Species = "T.spiralis")

df_both_species = bind_rows(df_dense_nat, df_dense_ts)

nrow(df_both_species)
df_both_species |> count(Species, TE_type)


## Determine the timings of the TEs ----
# Do DNA transposons and retrotransposons differ in age, and is that pattern the same in T.nativa and T.spiralis?
# Do older TEs appear more decayed? Use the relationship between perc_div and perc_del

# Exploratory: Checking for TE divergence
plot_divergence_hist <- df_both_species |> 
  ggplot(aes(x = perc_div, fill = TE_type)) +
  geom_histogram(binwidth = 1, color = "white") +
  labs(x = "Divergence from reference sequence %", y = "Number of TE copies", title = "Distribution of TE divergence by type and species") +
  facet_grid(TE_type ~ Species, scales = 'free_y') +
  theme_light() +
  scale_fill_brewer(palette = "Pastel2") +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "Fig4_TE_Species_dist.png"), plot = plot_divergence_hist, width = 10, height = 6)


#==> T. nativa: no impossible values (all between ~2% and ~37% divergence); data passes the check.
#==> Retrotransposons: one main peak at ~20%; a few very old copies (~37%) separated from the main distribution,
#    and a thin tail of younger copies (2–8%) too small to indicate a recent burst. Outliers kept as genuine values.
#==> DNA transposons: patchy distribution (n = 85), a few younger copies at ~9–10%; gaps reflect small sample size.

# Plotting perc_div (substitution, as a proxy for age) against perc_del (deletion) for both species.
# To examine whether older TE copies carry more deletions, and whether this relationship differs between DNA transposons and Retrotransposon, and T.spiralis species and T.nativa species.
  
plot_divergence_deletion <- df_both_species |>
  ggplot(aes(x = perc_div, y = perc_del, colour = TE_type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm) +
  facet_wrap(~ Species) +
  labs(x = "Divergence from the reference sequence in %", y = "Deletion in %", color = "TE type")

ggsave(file.path(output_dir, "fig_div_vs_del.png"), plot = plot_divergence_deletion, width = 10, height = 6)
  
# In both species, older retrotransposon copies show slightly more deletion, whereas DNA transposons show no clear relationship; the overlapping confidence bands indicate the differences between TE types and species are weak.

# Part 2
# Are DNA transposons and retrotransposons different ages, and is that difference the same in T.native and T.spiralis?

#Statistic analysis: Using wilcoxon ranked sum tested, based on the histogram we generated we noticed that the distibution of DNA TEs is regular. With wilcoxon ranked sum we will check whether the divergence differ between the two classes of TE (DNA vs Retrotransposon).
  
wilcox.test(perc_div ~ TE_type, data = df_both_species |> filter(Species == "T.nativa"))
#In T.nativa, divergence did not differ between DNA transposons and retrotransposons (Wilcoxon rank-sum test, W = 29899, p = 0.48; medians ≈ 21% for both).

wilcox.test(perc_div ~ TE_type, data = df_both_species |> filter(Species == "T.spiralis"))
    
# DNA transposons and retrotransposons did not differ in divergence in T.nativa, but DNA transposons were slightly more diverged in T.spiralis (W = 41606, p = 0.038). The T. spiralis difference is small (~1.4 percentage points), so the TE-type age pattern is broadly similar between species.      #=======> The small T.spiralis difference (p = 0.038) does not survive correction for running two tests, so it should be treated as a tentative trend, not a firm result.
      
#  Comparing divergence between species within each TE type, to test directly whether TE timing differs between T.nativa and T.spiralis.
      
wilcox.test(perc_div ~ Species, data = df_both_species |> filter(TE_type == "Retrotransposons"))

wilcox.test(perc_div ~ Species, data = df_both_species |> filter(TE_type == "DNA Transposons"))
      
df_median_div = df_both_species |>
  group_by(Species, TE_type) |>
  summarise(med = median(perc_div, na.rm = TRUE),
            n = sum(!is.na(perc_div)), .groups = "drop")

# EDIT 3: captions should reflect the data being analyzed, not hard-coded results.
median_caption <- paste(
  sprintf("%s | %s: median %.1f%% (n = %d)",
          df_median_div$Species, df_median_div$TE_type,
          df_median_div$med, df_median_div$n),
  collapse = "\n"
)
      
#=====> Medians were calculated per species and TE type to report the size of the age difference alongside the Wilcoxon p-values, since a p-value alone shows whether groups differ, not by how much. ANd to add captions to the Boxplot.
      
df_both_species |> filter(Species == "T.nativa") |>
  ggplot(aes(x = TE_type, y = perc_div, fill = TE_type)) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.5, colour = "grey30") +
  theme_light()
  
df_both_species |> filter(Species == "T.spiralis") |>
  ggplot(aes(x = TE_type, y = perc_div, fill = TE_type)) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.5, colour = "grey30") +
  theme_light()
  
plot_violin <- df_both_species |>
  ggplot(aes(x = TE_type, y = perc_div, fill = TE_type)) +
  geom_violin(alpha = 0.4) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.2, size = 0.8, colour = "grey30") +
  facet_wrap(~ Species) +
  labs(x = "TE type", y = "Divergence from reference sequence (%)", title = "TE divergence by type and species") +
  labs(caption = median_caption) +
    theme(legend.position = "none")
    
ggsave(file.path(output_dir, "fig_div_TEtype.png"), plot = plot_violin, width = 10, height = 6)
    
#Histograms showed irregular distributions for DNA transposons and unequal group sizes, so a non-parametric Wilcoxon test was used instead of a t-test.

#==> Median divergence is around 21% for both TE types in T.nativa, so DNA transposons and retrotransposons are of similar age. In T.spiralis, DNA transposons sit slightly higher (~23% vs ~21%), suggesting slightly older DNA TE copies.
#==> Retrotransposons look almost identical in both species (median ~21%, range ~3–38%), and span a wider range of ages than DNA TEs, i.e. activity over a longer period.
#==> All violins have a single bulge at ~20–25% and no peak at low divergence: no evidence of a recent burst of TE activity in either species; most copies are old.
#==> DNA groups are small (n = 85, 100) vs retros (n = 738 each).
  
# Part 3

#Do older TE copies tend to be shorter, and does that differ between DNA transposons and retrotransposons, and between T. nativa and T. spiralis?

plot_divergence_length <- df_both_species |>
  ggplot(aes(x = perc_div, y = qry_width, colour = TE_type)) +
  scale_y_log10() +
  geom_point(size = 0.8, alpha = 0.3) +
  geom_smooth(method = lm) +
  facet_grid(TE_type ~ Species, scales = "free_y") +
  labs(x = "Divergence from the reference sequence (%)",
       y = "TE copy length (bp, log scale)",
       title = "Older TE copies appear longer",
       caption = "Each point is one TE copy; lines are linear fits with 95% confidence bands.\nn: T. nativa DNA = 85, Retro = 738; T. spiralis DNA = 100, Retro = 738") +
  theme_light() +
  theme(legend.position = "none")
    
ggsave(file.path(output_dir, "fig_div_TEwidth_4.png"), plot = plot_divergence_length, width = 10, height = 8)
    
#==> Most TE copies are short fragments (<200 bp), far below full length TE size; consistent with old, decayed TEs. Contrary to expectation, copy length increases with divergence, clearly for retrotransposons (both species, n = 738), weakly for DNA transposons (wide confidence bands, n = 85, 100).
#==> Likely bias: highly diverged copies must be longer to pass RepeatMasker's score threshold, so old short fragments are under-detected. Length vs divergence therefore reflects detectability as well as decay.
    
# Part 4: Do contigs in T.nativa and T.spiralis differ in their composition of TE families?
# PCA on TE family copy counts per contig (one row per contig, one column per family),
# coloured by species, to see whether contigs from the two species carry similar or different TE family mixes.

# Creating a wide table to run PCA
    
df_pca_wide = df_both_species |>
  group_by(qry_id, Species, matching_class) |>
  summarise(count = n(), .groups = "drop") |>
  pivot_wider(names_from = matching_class,
              values_from = count,
              values_fill = 0)

dim(df_pca_wide)

# EDIT 3 (PCA): remove zero-variance families before standardization.
# This prevents prcomp(scale. = TRUE) from failing on constant columns.
pca_features <- df_pca_wide |>
  select(-qry_id, -Species)
if (anyNA(pca_features)) {
  stop("PCA feature matrix contains missing values; inspect the TE count table.")
}
feature_sd <- vapply(pca_features, sd, numeric(1))
keep_features <- is.finite(feature_sd) & feature_sd > 0
if (sum(keep_features) < 2L) {
  stop("At least two non-constant TE families are needed for a two-axis PCA.")
}
pca_input <- pca_features[, keep_features, drop = FALSE]
if (any(!keep_features)) {
  message("PCA removed constant TE families: ",
          paste(names(feature_sd)[!keep_features], collapse = ", "))
}

pca_te <- prcomp(pca_input, center = TRUE, scale. = TRUE)

# The explained variance must be calculated from the CURRENT PCA result.
pca_variance_percent <- 100 * pca_te$sdev^2 / sum(pca_te$sdev^2)
summary(pca_te)

# # Table of PC scores for each contig (from pca_te$x), with species added, for plotting PC1 vs PC2.

df_pca_score = pca_te$x |>
  as_tibble() |>
  mutate(Species = df_pca_wide$Species)

df_pca_score
    
#Plotting
plot_pca <- df_pca_score |>
  ggplot(aes(x = PC1, y = PC2, colour = Species)) +
  geom_point(alpha = 0.5) +
stat_ellipse(level = 0.95, segments = 51, type = "t") +
  labs(x = sprintf("PC1 (%.1f%% of variance)", pca_variance_percent[1]),
       y = sprintf("PC2 (%.1f%% of variance)", pca_variance_percent[2]),
       title = "PCA of TE-family composition by species") +
  labs(caption = sprintf("Each point is one contig; 95%% confidence ellipses by species.\nPCA on scaled TE family counts (%d non-constant families).", ncol(pca_input))) +
  theme_light()
  
pca_te$rotation[ ,1:2]
  
ggsave(file.path(output_dir, "PCA_species.png"), plot = plot_pca, width = 8, height = 6)