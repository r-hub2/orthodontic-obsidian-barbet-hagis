## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

required <- c("ape", "vegan", "dplyr")

if (!all(sapply(required, requireNamespace, quietly = TRUE)))
  knitr::opts_chunk$set(eval = FALSE)
data.table::setDTthreads(1L)

## ----libraries, message=FALSE, warning=FALSE----------------------------------
library("ape")
library("vegan")
library("dplyr")
library("hagis")
library("ggplot2")

## ----load-data----------------------------------------------------------------
head(P_sojae_survey) # survey sample data

head(sample_meta) # metatada about the sample collection locations

## ----clean-data---------------------------------------------------------------
P_sojae_survey$Isolate <-
  gsub(pattern = "MPS17_",
       replacement = "",
       x = P_sojae_survey$Isolate)
P_sojae_survey$Rps <-
  gsub(pattern = "Rps ",
       replacement = "",
       x = P_sojae_survey$Rps)

## ----hagis-arguments----------------------------------------------------------
hagis_args <- list(
  x = P_sojae_survey,
  cutoff = 60,
  control = "susceptible",
  sample = "Isolate",
  gene = "Rps",
  perc_susc = "perc.susc"
)

## ----create-matrix------------------------------------------------------------
P_sojae_survey.matrix <- do.call(create_binary_matrix, hagis_args)

P_sojae_survey.matrix

## ----distance-matrix----------------------------------------------------------
P_sojae_survey.matrix.jaccard <-
  vegdist(P_sojae_survey.matrix, "jaccard", na.rm = TRUE)

## ----barplot------------------------------------------------------------------
princoor.pathotype <- pcoa(P_sojae_survey.matrix.jaccard) 

barplot(princoor.pathotype$values$Relative_eig[1:10])

## ----axis-percent-------------------------------------------------------------
# Dimension (i.e., Axis 1 (PCOA1))
Axis1.percent <-
  princoor.pathotype$values$Relative_eig[[1]] * 100

# Dimension (i.e., Axis 2 (PCOA2))
Axis2.percent <-
  princoor.pathotype$values$Relative_eig[[2]] * 100

Axis1.percent

Axis2.percent

## ----create-df----------------------------------------------------------------
princoor.pathotype.data <-
  data.frame(
    Sample = as.integer(rownames(princoor.pathotype$vectors)),
    X = princoor.pathotype$vectors[, 1],
    Y = princoor.pathotype$vectors[, 2]
  )

## ----merge-metadata-----------------------------------------------------------
princoor.pathotype.data <-
  left_join(princoor.pathotype.data, sample_meta, by = "Sample")

princoor.pathotype.data

## ----ggplot-pca---------------------------------------------------------------
ggplot(data = princoor.pathotype.data, aes(x = X, y = Y)) +
  geom_point(aes(colour = Locale)) +
  xlab(paste("PCOA1 - ", round(Axis1.percent, 2), "%", sep = "")) +
  ylab(paste("PCOA2 - ", round(Axis2.percent, 2), "%", sep = "")) +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.key.size = unit(1, 'lines')
  ) +
  stat_ellipse(data = princoor.pathotype.data, aes(x = X, y = Y),
               level = 0.95) +
  ggtitle("Pathotype Jaccard Distances PCOA")

## ----create-group-lists-------------------------------------------------------
groups <- factor(c(rep("Michigan_1", 11), rep("Michigan_2", 10)))

# this number shows how many isolates are in all "groups" lists combined
length(groups)

# this shows the number of isolates within your data set, these numbers should
# match for downstream analyses to work!! 
length(unique(P_sojae_survey$Isolate))

## ----beta-dispersion----------------------------------------------------------
 # calculates the beta-dispersion for each group, when comparing 2 or more
pathotype.disp <-
  betadisper(P_sojae_survey.matrix.jaccard, groups)

# tests if centroid distances are significantly different from each other
pathotype.disp.anova <- anova(pathotype.disp) 
pathotype.disp.anova

# test significance between each group
pathotype.disp.TukeyHSD <- TukeyHSD(pathotype.disp)
pathotype.disp.TukeyHSD

# plot showing the dispersion for each group
plot(pathotype.disp, hull = FALSE, ellipse = TRUE)

## ----adonis-------------------------------------------------------------------
pathotype.adonis <- adonis2(P_sojae_survey.matrix.jaccard ~ groups)

pathotype.adonis

## ----anosim-------------------------------------------------------------------
pathotype.anosim <- anosim(P_sojae_survey.matrix.jaccard, groups)

pathotype.anosim

