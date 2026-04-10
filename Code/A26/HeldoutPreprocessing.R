# 1. remove outliers (only Metabolitics data, because Proteomics outliers looks like another batch)
# 2. comparing the mapped samples and heldout samples and swap them if necessary
# 3. compute k = (x-med)/mad, for Proteomics and metabolitics data, each of them would have one or two k values 
# (from original samples and heldout samples), and we unselect the one with maximum k values for the final mapping samples

### Use both PC1 and PC2, and compute the robust distance of each sample from the center:
# d_i = sqrt{(PC1_i - med(PC1))^2 + (PC2_i - med(PC2))^2}
# Then normalize this by the MAD of the distances:
# k_i = \frac{d_i}{MAD(d)}
# This way, k represents how far a sample is from the robust center of the PCA space (PC1, PC2).


######################################## 
############ Metabolomics Data #########
######################################## 
a26.met.data <- read_excel("../data/PostCombat_A26_Metabolite_beforeHeldout.xlsx")


#Function to detect outliers based on median ± k * MAD （detect outliers based on the median ± k × MAD (Median Absolute Deviation)）
mad_outliers <- function(x, k = 10) {
  med <- median(x, na.rm = TRUE)
  m <- mad(x, constant = 1, na.rm = TRUE)   # default constant=1
  out <- x[abs(x - med) > k * m]
  return(out)
}


pca_met_corrected <- prcomp(a26.met.data[,-c(1:3)], scale. = TRUE)
pca_met <- data.frame(pca_met_corrected$x)  # PCA scores
pca_met$plate <- a26.met.data$Plate
pca_met$sample_type <- a26.met.data$sampleType
pca_met$color_indicator <- ifelse(pca_met$sample_type == "SPQC", "black", pca_met$plate)
pca_met$sample <- as.numeric(a26.met.data$Sample)



# Detect outliers for PC1 and PC2
pc1_outliers <- mad_outliers(pca_met$PC1, k = 7)
pc2_outliers <- mad_outliers(pca_met$PC2, k = 7)
# Subset outliers
outlier_met <- pca_met[pca_met$PC1 %in% pc1_outliers | pca_met$PC2 %in% pc2_outliers, ]




ggplot(pca_met, aes(x=PC1, y=PC2, color=color_indicator)) +
  geom_point(size=2) +
  geom_text_repel(data=outlier_met, aes(label=sample), size=3, color="black") +  # Add outlier labels
  labs(title="PCA Plot of A26 proabolome Data using [ComBat]",
       x="PC 1",
       y="PC 2") +
  theme_minimal() +
  scale_color_manual(values=c("black", "red", "blue", "green", "purple", "orange", "pink")) +  # Set black for SPQC
  scale_shape_manual(values=c(16, 17))  

# remove 4 outliers: 2, 8, 246, 413
a26.met.data <- a26.met.data[!a26.met.data$Sample %in% outlier_met$sample,  ]
pca_met <- pca_met[!pca_met$sample %in% outlier_met$sample, ]

######################################## 
############ Proteomics Data ########### 
######################################## 

a26.pro.data <- read_excel("../data/PostCombat_A26_Protein.xlsx", sheet = "combat_corrected_unscaled")

# PCA + PLOT
pca_pro <- prcomp(a26.pro.data[,-c(1:3)], scale. = TRUE)
pca_pro <- data.frame(pca_pro$x)  # PCA scores
pca_pro$plate <- a26.pro.data$pro.plate
pca_pro$sample_type <- a26.pro.data$pro.sampletype
pca_pro$color_indicator <- ifelse(pca_pro$sample_type == "SPQC", "black", pca_pro$plate)
pca_pro$sample <- a26.pro.data$pro.sample

# convert Duke ID to Analysis ID
dukeid.df <- data.frame(read_excel("../data/A26_Bonnie_Flatmap_Grid_20260315.xlsx", sheet = "Proteome_DukeIDs"))
colnames(pca_pro)[431] <- "DukeID"
pca_pro$DukeID <- as.numeric(substr(pca_pro$DukeID, 3,8))
pca_df_mapped <- pca_pro %>%
  left_join(dukeid.df, by = "DukeID")

# Remove SPQC samples
pca_df_mapped <- pca_df_mapped[pca_df_mapped$sample_type == "Non-SPQC",]

# Remove NA in AnalysisID
pca_pro <- pca_df_mapped
pca_pro <- pca_pro[!is.na(pca_pro$AnalysisID), ]


# Detect outliers for PC1 and PC2
pc1_outliers <- mad_outliers(pca_pro$PC1, k = 7)
pc2_outliers <- mad_outliers(pca_pro$PC2, k = 7)
# Subset outliers
outlier_pro <- pca_pro[pca_pro$PC1 %in% pc1_outliers | pca_pro$PC2 %in% pc2_outliers, ]



ggplot(pca_pro, aes(x=PC1, y=PC2, color=color_indicator)) +
  geom_point(size=2) +
  geom_text_repel(data=outlier_pro, aes(label=AnalysisID), size=3, color="black") +  # Add outlier labels
  labs(title="PCA Plot of A26 proabolome Data using [ComBat]",
       x="PC 1",
       y="PC 2") +
  theme_minimal() +
  scale_color_manual(values=c("black", "red", "blue", "green", "purple", "orange", "pink")) +  # Set black for SPQC
  scale_shape_manual(values=c(16, 17))  


######################################## 
# Compare the heldout and mapped samples and swap them if necessary
######################################## 

# 397 common samples in both datasets
common_samples <- intersect(pca_pro$AnalysisID, pca_met$sample)
# update both pca data
pca_pro <- pca_pro[pca_pro$AnalysisID %in% common_samples,]
pca_met <- pca_met[pca_met$sample %in% common_samples,]


# HELDOUT: sample 413 is removed and replaced by the heldout1 sample 247
heldout.df <- data.frame(read_excel("../data/A26_Bonnie_Flatmap_Grid_20260315.xlsx", sheet = "Mapped_Samples"))
heldout.df <- heldout.df[-nrow(heldout.df),]

heldout.df[heldout.df$Mapped_AnalysisID %in% outlier_met$sample, ]

idx <- heldout.df$Mapped_AnalysisID %in% outlier_met$sample
heldout.df$Mapped_AnalysisID[idx] <- heldout.df$HeldOut_AnalysisID1[idx]
heldout.df$HeldOut_AnalysisID1[idx] <- NA


# keep only heldout rows where at least one candidate is in common_samples
heldout.df <- heldout.df %>%
  filter(
    Mapped_AnalysisID %in% common_samples |
      HeldOut_AnalysisID1 %in% common_samples |
      HeldOut_AnalysisID2 %in% common_samples
  )

# Function: robust distance in 2D
compute_k_2d <- function(df, pc1_col = "PC1", pc2_col = "PC2", id_col = "AnalysisID") {
  med1 <- median(df[[pc1_col]], na.rm = TRUE)
  med2 <- median(df[[pc2_col]], na.rm = TRUE)
  
  d <- sqrt((df[[pc1_col]] - med1)^2 + (df[[pc2_col]] - med2)^2)
  mad_d <- mad(d, na.rm = TRUE, constant = 1)
  k <- d / mad_d
  
  data.frame(AnalysisID = df[[id_col]], k = k)
}

# Apply to proteomics and metabolomics
k_pro <- compute_k_2d(pca_pro, "PC1", "PC2", "AnalysisID")
k_met <- compute_k_2d(pca_met, "PC1", "PC2", "sample")

# Merge k-values: max across modalities
k_df <- merge(k_pro, k_met, by = "AnalysisID", all = TRUE)
k_df$k_max <- pmax(k_df$k.x, k_df$k.y, na.rm = TRUE)


# Add a column for final selection
heldout.df$finalID <- heldout.df$Mapped_AnalysisID

for (i in seq_len(nrow(heldout.df))) {
  candidates <- c(
    heldout.df$Mapped_AnalysisID[i],
    heldout.df$HeldOut_AnalysisID1[i],
    heldout.df$HeldOut_AnalysisID2[i]
  )
  # Drop NAs
  candidates <- na.omit(candidates)
  
  # Get k-values
  k_vals <- k_df$k_max[match(candidates, k_df$AnalysisID)]
  
  # Keep the candidate with the smallest k
  if (length(k_vals) > 0) {
    best_idx <- which.min(k_vals)
    heldout.df$finalID[i] <- candidates[best_idx]
  } else {
    heldout.df$finalID[i] <- NA
  }
}

# Function to lookup k values
lookup_k <- function(ids) {
  k_df$k_max[match(as.character(ids), k_df$AnalysisID)]
}

# Add columns for each candidate’s k
k_table_wide <- heldout.df %>%
  mutate(
    k_mapped  = lookup_k(Mapped_AnalysisID),
    k_held1   = lookup_k(HeldOut_AnalysisID1),
    k_held2   = lookup_k(HeldOut_AnalysisID2)
  )


#write.xlsx(k_table_wide,"../data/Heldout_Updated_a26.xlsx")

