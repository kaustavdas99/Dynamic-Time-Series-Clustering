####################################################################################
# LOAD LIBRARIES
####################################################################################
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(dtwclust)
library(ggplot2)
library(plotly)
library(zoo)

####################################################################################
# READ DATA & PARSE MONTH
####################################################################################
immun <- read_excel("Child_Immunisation_FINAL.xlsx") %>%
  mutate(
    Month_clean = str_squish(Month),
    Month_start = str_extract(Month_clean, "^[A-Za-z]+\\s\\d{2}"),  # e.g., "Apr 20"
    Month_date  = parse_date_time(Month_start, orders = "b y")
  )

####################################################################################
# AGGREGATE DUPLICATES PER DISTRICT-MONTH
####################################################################################
ts_data <- immun %>%
  select(District, Month_date,
         `Percentage of BCG`,
         `Percentage of OPV 0 (Birth Dose)`,
         `Percentage of Hepatitis-B0 (Birth Dose)`) %>%
  group_by(District, Month_date) %>%
  summarise(
    BCG   = mean(`Percentage of BCG`, na.rm = TRUE),
    OPV0  = mean(`Percentage of OPV 0 (Birth Dose)`, na.rm = TRUE),
    HepB0 = mean(`Percentage of Hepatitis-B0 (Birth Dose)`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Month_date >= as.Date("2020-04-01"),
         Month_date <= as.Date("2022-03-01"))

####################################################################################
# PIVOT TO WIDE FORMAT FOR CLUSTERING & IMPUTE MISSING
####################################################################################
ts_wide <- ts_data %>%
  pivot_wider(names_from = Month_date, values_from = c(BCG, OPV0, HepB0))

ts_wide_imputed <- ts_wide
for(i in 2:ncol(ts_wide_imputed)) {
  x <- ts_wide_imputed[[i]]
  na_idx <- is.na(x)
  if(any(na_idx)) {
    ts_wide_imputed[[i]][na_idx] <- approx(which(!na_idx), x[!na_idx], which(na_idx), method="linear", rule=2)$y
  }
}

####################################################################################
# PREPARE MATRIX & DTW CLUSTERING
####################################################################################
ts_matrix <- as.matrix(ts_wide_imputed[,-1])
rownames(ts_matrix) <- ts_wide_imputed$District

set.seed(123)
ts_clusters <- tsclust(
  ts_matrix,
  type = "partitional",
  k = 4,
  distance = "dtw_basic",
  centroid = "pam",
  seed = 123,
  trace = TRUE
)

cluster_assignment <- data.frame(
  District = rownames(ts_matrix),
  Cluster  = ts_clusters@cluster
)

####################################################################################
# COMPUTE AVERAGE COVERAGE PER CLUSTER
####################################################################################
avg_ts <- ts_data %>%
  left_join(cluster_assignment, by = "District") %>%
  group_by(Cluster, Month_date) %>%
  summarise(
    BCG   = mean(BCG, na.rm = TRUE),
    OPV0  = mean(OPV0, na.rm = TRUE),
    HepB0 = mean(HepB0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(BCG, OPV0, HepB0),
               names_to = "Vaccine",
               values_to = "Percent")

####################################################################################
# FULL MONTH SEQUENCE FOR ALL CLUSTERS & VACCINES
####################################################################################
all_months <- seq.Date(from = as.Date("2020-04-01"), to = as.Date("2022-03-01"), by = "month")

avg_ts_full_complete <- expand.grid(
  Cluster    = unique(avg_ts$Cluster),
  Vaccine    = unique(avg_ts$Vaccine),
  Month_date = all_months
) %>%
  left_join(avg_ts, by = c("Cluster","Vaccine","Month_date")) %>%
  arrange(Cluster, Vaccine, Month_date) %>%
  group_by(Cluster, Vaccine) %>%
  mutate(Percent = na.locf(Percent, na.rm = FALSE, fromLast = FALSE)) %>%
  ungroup()

###################################################################################
# PLOT INTERACTIVE COVERAGE PER CLUSTER (APR 2020 - MAR 2022)
###################################################################################
plots_list <- lapply(sort(unique(avg_ts_full_complete$Cluster)), function(cl) {
  df <- avg_ts_full_complete %>% filter(Cluster == cl)
  
  plot_ly(df, x = ~Month_date, y = ~Percent*100, color = ~Vaccine,
          type = 'scatter', mode = 'lines+markers',
          hoverinfo = 'text',
          text = ~paste("Vaccine:", Vaccine,
                        "<br>Coverage:", round(Percent*100,1), "%",
                        "<br>Month:", format(Month_date, "%b %Y"))) %>%
    layout(title = list(text = paste0("Cluster ", cl),
                        x = 0.5,       # Center the title
                        xanchor = "center"),
           xaxis = list(title = "Month",
                        tickvals = all_months,
                        ticktext = format(all_months, "%b %Y"),
                        tickangle = -45,
                        range = c(min(all_months), max(all_months))),
           yaxis = list(title = "Coverage (%)", range = c(0, 100)),
           showlegend = FALSE,
           margin = list(l = 60, r = 60, t = 80, b = 60))  # Add margins
})

p_final <- subplot(plots_list, nrows = length(plots_list),
                   shareX = TRUE, shareY = TRUE, titleY = TRUE) %>%
  layout(title = list(text = "Time Series of Vaccine Coverage by Cluster (BCG, OPV0, HepB0) Apr 2020 - Mar 2022",
                      x = 0.5, xanchor = "center"),  # Center main title
         showlegend = TRUE,
         legend = list(title = list(text = "Vaccine")),
         margin = list(l = 70, r = 70, t = 100, b = 70))  # Global margins

p_final


####################################################################################
# CLUSTER NAME LOOKUP (DERIVED FROM AVG COVERAGE & RISK)
####################################################################################

cluster_names <- c(
  "1" = "Lower Coverage (Medium Risk)",
  "2" = "Moderate Coverage (Medium Risk)",
  "3" = "Very High Coverage (Low Risk)",
  "4" = "High Coverage (Low Risk)"
)

####################################################################################
# PLOT INTERACTIVE COVERAGE PER CLUSTER (WITH CLUSTER NAMES)
####################################################################################
plots_list <- lapply(sort(unique(avg_ts_full_complete$Cluster)), function(cl) {
  
  df <- avg_ts_full_complete %>% filter(Cluster == cl)
  
  plot_ly(
    df,
    x = ~Month_date,
    y = ~Percent * 100,
    color = ~Vaccine,
    type = "scatter",
    mode = "lines+markers",
    hoverinfo = "text",
    text = ~paste(
      "Cluster:", cluster_names[as.character(cl)],
      "<br>Vaccine:", Vaccine,
      "<br>Coverage:", round(Percent * 100, 1), "%",
      "<br>Month:", format(Month_date, "%b %Y")
    )
  ) %>%
    layout(
      title = list(
        text = paste0("Cluster ", cl, ": ", cluster_names[as.character(cl)]),
        x = 0.5,
        xanchor = "center"
      ),
      xaxis = list(
        title = "Month",
        tickvals = all_months,
        ticktext = format(all_months, "%b %Y"),
        tickangle = -45,
        range = c(min(all_months), max(all_months))
      ),
      yaxis = list(
        title = "Coverage (%)",
        range = c(0, 100)
      ),
      showlegend = FALSE,
      margin = list(l = 60, r = 60, t = 80, b = 60)
    )
})

p_final <- subplot(
  plots_list,
  nrows = length(plots_list),
  shareX = TRUE,
  shareY = TRUE,
  titleY = TRUE
) %>%
  layout(
    title = list(
      text = "Time Series of Vaccine Coverage by Cluster (BCG, OPV0, HepB0), Apr 2020 – Mar 2022",
      x = 0.5,
      xanchor = "center"
    ),
    showlegend = TRUE,
    legend = list(title = list(text = "Vaccine")),
    margin = list(l = 70, r = 70, t = 100, b = 70)
  )

p_final


#library(htmlwidgets)

saveWidget(
  p_final,
  "DTW_Cluster_Vaccine_Coverage.html",
  selfcontained = TRUE
)

orca(
  p_final,
  "DTW_Cluster_Vaccine_Coverage.png",
  width = 1400,
  height = 2000
)


####################################################################################
# ASSIGN RISK LABELS TO CLUSTERS
####################################################################################


cluster_avg <- avg_ts %>%
  group_by(Cluster) %>%
  summarise(Avg_Coverage=mean(Percent, na.rm=TRUE)) %>%
  mutate(Risk=case_when(
    Avg_Coverage < 0.7 ~ "High Risk",
    Avg_Coverage >= 0.7 & Avg_Coverage < 0.9 ~ "Medium Risk",
    Avg_Coverage >= 0.9 ~ "Low Risk"
  ))

district_risk <- cluster_assignment %>%
  left_join(cluster_avg, by="Cluster") %>%
  arrange(Risk, Cluster)

print(district_risk)
write.csv(district_risk, "District_Risk.csv", row.names = FALSE)


###################################################################################

###################################################################################
# STATISTICAL TESTS & CLUSTER VALIDATION
###################################################################################

library(dtwclust)
library(pvclust)
library(cluster)
library(dplyr)

# ----------------------------------------
# Cophenetic Correlation Coefficient (CCC)
# ----------------------------------------

# Hierarchical clustering from DTW distance
hc <- hclust(dtw_dist, method = "average")

# Cophenetic distances
cophe_dist <- cophenetic(hc)

# Compute CCC
ccc <- cor(as.vector(dtw_dist), as.vector(cophe_dist))
cat("Cophenetic Correlation Coefficient (CCC):", round(ccc, 3), "\n\n")

# ----------------------------------------
# Cluster Stability via Bootstrap (pvclust)
# ----------------------------------------

# pvclust- rows = features (time points), columns = observations (districts)
ts_matrix_t <- t(ts_matrix)

set.seed(123)
pv <- pvclust(
  ts_matrix_t,
  method.hclust = "average", 
  method.dist = "euclidean",  
  nboot = 1000
)

# Plot dendrogram with AU values
plot(pv)
pvrect(pv, alpha = 0.95)

# Extract AU values and interpret stability
au_values <- data.frame(
  Cluster = rownames(pv$edges),
  AU_pval = pv$edges[, "au"]
) %>%
  mutate(
    Stability = case_when(
      AU_pval > 0.95 ~ "Very stable",
      AU_pval >= 0.8 & AU_pval <= 0.95 ~ "Moderately stable",
      AU_pval >= 0.6 & AU_pval < 0.8 ~ "Weakly stable",
      TRUE ~ "Unstable"
    )
  )

cat("\nCluster AU p-values and Stability:\n")
print(au_values)

####################################################################################
# Silhouette Analysis
####################################################################################

library(cluster)

# Convert DTW distance to matrix 
dtw_dist_mat <- as.matrix(dtw_dist)

dist_districts <- as.dist(
  dtw_dist_mat[rownames(ts_matrix), rownames(ts_matrix)]
)

# Cluster labels
cluster_labels <- ts_clusters@cluster

# Compute silhouette
sil <- silhouette(cluster_labels, dist_districts)

# Per-district silhouette scores
sil_df <- data.frame(
  District = rownames(ts_matrix),
  Cluster = cluster_labels,
  Silhouette_Width = sil[, "sil_width"]
)

# Cluster-wise summary
cluster_scores <- sil_df %>%
  group_by(Cluster) %>%
  summarise(
    Num_Districts  = n(),
    Avg_Silhouette = mean(Silhouette_Width),
    Min_Silhouette = min(Silhouette_Width),
    Max_Silhouette = max(Silhouette_Width)
  )

print(sil_df)
print(cluster_scores)

# Plot
plot(sil, main = "Silhouette Plot of DTW Clustering")

par(mar = c(5, 4, 4, 2) + 0.1) 

plot(
  sil,
  main = "",
  xlab = "Silhouette width s\u1D62",
  ylab = "Districts"
)

title(
  main = "Silhouette Plot of DTW Clustering",
  font.main = 2,
  cex.main = 1.3
)



cat("\nPer-district Silhouette Scores:\n")
print(sil_df)

cat("\nPer-cluster Silhouette Summary:\n")
print(cluster_scores)



####################################################################################################
#Other Clustering plots

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(dtwclust)
library(ggplot2)
library(plotly)
library(zoo)
library(Rtsne)
library(umap)
library(pheatmap)
library(cluster)
library(RColorBrewer)

####################################################################################
# t-SNE visualization
####################################################################################
library(Rtsne)
library(ggplot2)
library(ggrepel)

# Scale matrix
ts_matrix_scaled <- scale(ts_matrix)

set.seed(123)
tsne_res <- Rtsne(
  ts_matrix_scaled,
  dims = 2,
  perplexity = 5,
  pca = TRUE,
  verbose = TRUE
)

tsne_df <- data.frame(
  District = rownames(ts_matrix_scaled),
  Dim1 = tsne_res$Y[,1],
  Dim2 = tsne_res$Y[,2],
  Cluster = factor(ts_clusters@cluster)
)

tsne_df

p_tsne <- ggplot(tsne_df, aes(Dim1, Dim2, color = Cluster)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = District), size = 3) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "t-SNE of Districts by Birth-Dose Immunisation Trajectories",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

print(p_tsne)


####################################################################################
# UMAP visualization
####################################################################################

library(ggplot2)
library(ggrepel)

set.seed(123)
umap_res <- umap(ts_matrix_scaled)

umap_df <- data.frame(
  District = rownames(ts_matrix_scaled),
  Dim1 = umap_res$layout[,1],
  Dim2 = umap_res$layout[,2],
  Cluster = factor(ts_clusters@cluster)
)

umap_df

p_umap <- ggplot(umap_df, aes(Dim1, Dim2, color = Cluster)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = District), size = 3) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "UMAP of Districts by Birth-Dose Immunisation Trajectories",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

print(p_umap)

############################################################################################
