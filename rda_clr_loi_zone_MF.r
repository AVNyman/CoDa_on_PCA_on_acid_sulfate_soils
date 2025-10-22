# RDA Analysis: CLR Chemical Data vs LOI + Zone
# Redundancy Analysis with both continuous (LOI) and categorical (Zone) predictors

# =============================================================================
# PACKAGE LOADING
# =============================================================================

# Required packages
required_packages <- c("vegan", "ggplot2", "dplyr", "RColorBrewer", "gridExtra")

# Check and load packages
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# =============================================================================
# DATA LOADING
# =============================================================================

# Load CLR-transformed chemical data
clr_chem <- read.csv("clrChem.csv")
cat("✓ Loaded clrChem.csv:", nrow(clr_chem), "rows x", ncol(clr_chem), "columns\n")

# Load site data with LOI and Zone
site_data <- read.csv("SiteData.csv")
cat("✓ Loaded SiteData.csv:", nrow(site_data), "rows x", ncol(site_data), "columns\n")

# =============================================================================
# DATA PREPARATION
# =============================================================================

# Merge datasets by ID
merged_data <- merge(clr_chem, site_data[, c("ID", "LOI", "Zone")], by = "ID")
cat("✓ Merged datasets:", nrow(merged_data), "rows\n")

# Extract CLR chemical data (exclude ID column)
clr_matrix <- merged_data[, grepl("^clr_", names(merged_data))]
cat("✓ CLR chemical matrix:", nrow(clr_matrix), "samples x", ncol(clr_matrix), "elements\n")

# Prepare predictor variables
env_data <- data.frame(
  LOI = merged_data$LOI,
  Zone = factor(merged_data$Zone)
)

# Check data structure
cat("✓ Environmental data:\n")
cat("  LOI range:", round(min(env_data$LOI, na.rm = TRUE), 2), "-", 
    round(max(env_data$LOI, na.rm = TRUE), 2), "\n")
cat("  Zone levels:", paste(levels(env_data$Zone), collapse = ", "), "\n")
cat("  Zone counts:", paste(table(env_data$Zone), collapse = ", "), "\n")

# Check for missing values
complete_cases <- complete.cases(env_data)
cat("Missing values in environmental data:", sum(!complete_cases), "\n")

# Remove samples with missing data if any
if(any(!complete_cases)) {
  clr_matrix <- clr_matrix[complete_cases, ]
  env_data <- env_data[complete_cases, ]
  merged_data <- merged_data[complete_cases, ]
  cat("✓ Removed", sum(!complete_cases), "samples with missing data\n")
}

cat("Final dataset:", nrow(clr_matrix), "samples\n")

# =============================================================================
# REDUNDANCY ANALYSIS (RDA) - FULL MODEL
# =============================================================================

cat("\n=== PERFORMING RDA ANALYSIS ===\n")

# Full model with both predictors
rda_full <- rda(clr_matrix ~ LOI + Zone, data = env_data)
cat("✓ Full RDA model completed\n")

# Individual models for comparison
rda_loi_only <- rda(clr_matrix ~ LOI, data = env_data)
rda_zone_only <- rda(clr_matrix ~ Zone, data = env_data)
cat("✓ Individual models completed\n")

# Print summary of full model
cat("\nFull Model Summary:\n")
print(summary(rda_full))

# =============================================================================
# STATISTICAL TESTING
# =============================================================================

cat("\n=== STATISTICAL TESTING ===\n")

# Test overall model significance
set.seed(123)
anova_full <- anova(rda_full, permutations = 999)
cat("✓ Overall model test:\n")
print(anova_full)

# Test each predictor individually (marginal effects)
set.seed(123)
anova_marginal <- anova(rda_full, by = "margin", permutations = 999)
cat("\n✓ Marginal effects test:\n")
print(anova_marginal)

# Test each predictor sequentially (conditional effects)
set.seed(123)
anova_terms <- anova(rda_full, by = "terms", permutations = 999)
cat("\n✓ Sequential (conditional) effects test:\n")
print(anova_terms)

# Test RDA axes
set.seed(123)
anova_axes <- anova(rda_full, by = "axis", permutations = 999)
cat("\n✓ Axis-wise testing:\n")
print(anova_axes)

# =============================================================================
# VARIANCE PARTITIONING
# =============================================================================

cat("\n=== VARIANCE PARTITIONING ===\n")

# Partition variance between LOI and Zone
varpart_result <- varpart(clr_matrix, 
                         env_data[, "LOI", drop = FALSE],
                         env_data[, "Zone", drop = FALSE])

cat("✓ Variance partitioning completed\n")
print(varpart_result)

# Plot variance partitioning
plot(varpart_result, digits = 2, 
     Xnames = c("LOI", "Zone"),
     bg = c("red", "blue"))
title("Variance Partitioning: LOI vs Zone")

# =============================================================================
# VARIANCE EXPLAINED CALCULATIONS
# =============================================================================

cat("\n=== VARIANCE EXPLAINED ===\n")

# Extract variance components
total_var <- sum(rda_full$CA$eig) + sum(rda_full$CCA$eig)
constrained_var <- sum(rda_full$CCA$eig)
rda1_var <- rda_full$CCA$eig[1]
rda2_var <- rda_full$CCA$eig[2]

# Calculate percentages
total_explained <- (constrained_var / total_var) * 100
rda1_explained <- (rda1_var / total_var) * 100
rda2_explained <- (rda2_var / total_var) * 100

cat("Total variance explained by LOI + Zone:", round(total_explained, 2), "%\n")
cat("Variance explained by RDA1:", round(rda1_explained, 2), "%\n")
cat("Variance explained by RDA2:", round(rda2_explained, 2), "%\n")

# Individual model comparisons
loi_only_var <- (sum(rda_loi_only$CCA$eig) / sum(rda_loi_only$CA$eig + rda_loi_only$CCA$eig)) * 100
zone_only_var <- (sum(rda_zone_only$CCA$eig) / sum(rda_zone_only$CA$eig + rda_zone_only$CCA$eig)) * 100

cat("LOI alone explains:", round(loi_only_var, 2), "%\n")
cat("Zone alone explains:", round(zone_only_var, 2), "%\n")

# Adjusted R-squared
cat("Adjusted R-squared (full model):", round(RsquareAdj(rda_full)$adj.r.squared, 4), "\n")

# =============================================================================
# ORDINATION PLOTS
# =============================================================================

cat("\n=== CREATING ORDINATION PLOTS ===\n")

# Get scores
site_scores <- scores(rda_full, display = "sites", choices = 1:2)
species_scores <- scores(rda_full, display = "species", choices = 1:2)
biplot_scores <- scores(rda_full, display = "bp", choices = 1:2)

# Create data frame for sites
site_df <- data.frame(
  ID = merged_data$ID,
  RDA1 = site_scores[, 1],
  RDA2 = site_scores[, 2],
  LOI = env_data$LOI,
  Zone = env_data$Zone
)

# Create data frame for species
species_df <- data.frame(
  Element = gsub("clr_", "", rownames(species_scores)),
  RDA1 = species_scores[, 1],
  RDA2 = species_scores[, 2],
  RDA1_abs = abs(species_scores[, 1]),
  RDA2_abs = abs(species_scores[, 2])
)

# Sort by RDA1 importance
species_df <- species_df[order(species_df$RDA1_abs, decreasing = TRUE), ]

cat("Top 10 elements most strongly associated with RDA1:\n")
print(species_df[1:10, c("Element", "RDA1", "RDA2")])

# =============================================================================
# ENHANCED VISUALIZATIONS
# =============================================================================

# Set up color palette for zones
zone_colors <- RColorBrewer::brewer.pal(length(levels(env_data$Zone)), "Set1")
names(zone_colors) <- levels(env_data$Zone)

# Plot 1: Sites colored by Zone, sized by LOI
p1 <- ggplot(site_df, aes(x = RDA1, y = RDA2)) +
  geom_point(aes(color = Zone, size = LOI), alpha = 0.7) +
  scale_color_manual(values = zone_colors) +
  scale_size_continuous(name = "LOI", range = c(2, 6)) +
  labs(title = "RDA: Sites by Zone and LOI",
       x = paste0("RDA1 (", round(rda1_explained, 1), "%)"),
       y = paste0("RDA2 (", round(rda2_explained, 1), "%)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Plot 2: Sites colored by LOI, shaped by Zone
p2 <- ggplot(site_df, aes(x = RDA1, y = RDA2)) +
  geom_point(aes(color = LOI, shape = Zone), size = 3, alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red", name = "LOI") +
  scale_shape_manual(values = c(16, 17, 18)[1:length(levels(env_data$Zone))]) +
  labs(title = "RDA: Sites by LOI and Zone",
       x = paste0("RDA1 (", round(rda1_explained, 1), "%)"),
       y = paste0("RDA2 (", round(rda2_explained, 1), "%)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Plot 3: Species loadings
p3 <- ggplot(species_df, aes(x = RDA1, y = RDA2)) +
  geom_point(color = "darkgreen", size = 2) +
  geom_text(aes(label = Element), vjust = -0.5, hjust = 0.5, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = "Chemical Element Loadings",
       x = paste0("RDA1 (", round(rda1_explained, 1), "%)"),
       y = paste0("RDA2 (", round(rda2_explained, 1), "%)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Plot 4: Biplot with environmental vectors
p4 <- ggplot(site_df, aes(x = RDA1, y = RDA2)) +
  geom_point(aes(color = Zone), size = 3, alpha = 0.7) +
  scale_color_manual(values = zone_colors) +
  # Add environmental vectors
  geom_segment(data = data.frame(biplot_scores), 
               aes(x = 0, y = 0, xend = RDA1*2, yend = RDA2*2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "red", size = 1) +
  geom_text(data = data.frame(biplot_scores), 
            aes(x = RDA1*2.2, y = RDA2*2.2, label = rownames(biplot_scores)),
            color = "red", size = 4, fontface = "bold") +
  labs(title = "RDA Biplot with Environmental Vectors",
       x = paste0("RDA1 (", round(rda1_explained, 1), "%)"),
       y = paste0("RDA2 (", round(rda2_explained, 1), "%)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Display plots
print(p1)
print(p2)
print(p3)
print(p4)

# =============================================================================
# ZONE-SPECIFIC ANALYSIS
# =============================================================================

cat("\n=== ZONE-SPECIFIC ANALYSIS ===\n")

# LOI distribution by zone
loi_by_zone <- site_df %>%
  group_by(Zone) %>%
  summarise(
    n = n(),
    LOI_mean = mean(LOI),
    LOI_sd = sd(LOI),
    LOI_min = min(LOI),
    LOI_max = max(LOI),
    .groups = "drop"
  )

cat("LOI distribution by Zone:\n")
print(loi_by_zone)

# Boxplot of LOI by Zone
p5 <- ggplot(site_df, aes(x = Zone, y = LOI, fill = Zone)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = zone_colors) +
  labs(title = "LOI Distribution by Zone",
       x = "Zone", y = "Loss on Ignition (LOI)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

print(p5)

# Correlation between RDA axes and predictors
correlations <- data.frame(
  RDA1_LOI = cor(site_df$RDA1, site_df$LOI),
  RDA2_LOI = cor(site_df$RDA2, site_df$LOI)
)

cat("\nCorrelations between RDA axes and LOI:\n")
print(correlations)

# =============================================================================
# INTERACTION ANALYSIS
# =============================================================================

cat("\n=== TESTING LOI × ZONE INTERACTION ===\n")

# Model with interaction
rda_interaction <- rda(clr_matrix ~ LOI * Zone, data = env_data)

# Test interaction significance
set.seed(123)
anova_interaction <- anova(rda_interaction, rda_full, permutations = 999)
cat("Interaction test (comparing LOI*Zone vs LOI+Zone):\n")
print(anova_interaction)

# =============================================================================
# SUMMARY RESULTS
# =============================================================================

cat("\n=== COMPREHENSIVE ANALYSIS SUMMARY ===\n")
cat("✓ RDA Analysis Complete\n")
cat("✓ Total variance explained:", round(total_explained, 2), "%\n")
cat("✓ Adjusted R²:", round(RsquareAdj(rda_full)$adj.r.squared, 4), "\n")

cat("\nStatistical Significance:\n")
cat("- Overall model p-value:", round(anova_full$`Pr(>F)`[1], 4), "\n")
cat("- LOI marginal effect p-value:", round(anova_marginal$`Pr(>F)`[1], 4), "\n")
cat("- Zone marginal effect p-value:", round(anova_marginal$`Pr(>F)`[2], 4), "\n")

cat("\nVariance Partitioning:\n")
cat("- Pure LOI effect:", round(varpart_result$part$indfract$Adj.R.squared[1] * 100, 2), "%\n")
cat("- Pure Zone effect:", round(varpart_result$part$indfract$Adj.R.squared[2] * 100, 2), "%\n")
cat("- Shared LOI-Zone effect:", round(varpart_result$part$indfract$Adj.R.squared[3] * 100, 2), "%\n")

if(anova_interaction$`Pr(>F)`[2] < 0.05) {
  cat("✓ LOI × Zone interaction is SIGNIFICANT\n")
} else {
  cat("✗ LOI × Zone interaction is NOT significant\n")
}

# Save key results
cat("\n=== SAVED OBJECTS ===\n")
cat("- rda_full: Full RDA model (LOI + Zone)\n")
cat("- rda_loi_only: LOI-only model\n")
cat("- rda_zone_only: Zone-only model\n")
cat("- varpart_result: Variance partitioning results\n")
cat("- site_df: Site scores with environmental data\n")
cat("- species_df: Element loadings ranked by importance\n")
cat("- anova_marginal: Marginal effects test results\n")