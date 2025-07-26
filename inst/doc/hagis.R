## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)
data.table::setDTthreads(1L)

## ----load_data----------------------------------------------------------------
library("hagis")
head(P_sojae_survey)

## ----remove-gene--------------------------------------------------------------
P_sojae_survey$Rps <-
  gsub(
    pattern = "Rps ",
    replacement = "",
    x = P_sojae_survey$Rps
  )
head(P_sojae_survey)

## ----example-function, eval=FALSE---------------------------------------------
# Rps.summary <- summarize_gene(
#   x = P_sojae_survey,
#   cutoff = 60,
#   control = "susceptible",
#   sample = "Isolate",
#   gene = "Rps",
#   perc_susc = "perc.susc"
# )

## ----shared-args--------------------------------------------------------------
hagis_args <- list(
  x = P_sojae_survey,
  cutoff = 60,
  control = "susceptible",
  sample = "Isolate",
  gene = "Rps",
  perc_susc = "perc.susc"
)

## ----echo=TRUE----------------------------------------------------------------
Rps.summary <- do.call(summarize_gene, hagis_args)

Rps.summary

## ----pander-print-Rps, echo=TRUE----------------------------------------------
library(pander)

pander(Rps.summary)

## ----plot-summary, echo=TRUE--------------------------------------------------
autoplot(Rps.summary, type = "percentage")

autoplot(Rps.summary, type = "count")

## ----complexities, echo=TRUE, message=FALSE, warning=FALSE--------------------
complexities <- do.call(calculate_complexities, hagis_args)

complexities

## ----pander-print-complexities------------------------------------------------
pander(complexities$grouped_complexities)

pander(complexities$individual_complexities)

## ----summary-complexities-----------------------------------------------------
pander(summary(complexities))

## ----complexities-plot--------------------------------------------------------
autoplot(complexities, type = "percentage")

autoplot(complexities, type = "count")

## ----calculate-diversities, echo=TRUE-----------------------------------------
diversity <- do.call(calculate_diversities, hagis_args)
diversity

## ----diversity-pander---------------------------------------------------------
pander(diversity)

## ----diversities-table--------------------------------------------------------
diversities_table(diversity)

## ----individual-pathotypes----------------------------------------------------
individual_pathotypes(diversity)

## ----set-up-adv.plot----------------------------------------------------------
library(ggplot2)

Rps.plot <- autoplot(Rps.summary, type = "percentage")

Rps.plot

## ----change-plot-theme--------------------------------------------------------
Rps.plot <- Rps.plot +
  theme_minimal()

Rps.plot

## ----change-plot-font---------------------------------------------------------
Rps.plot <- Rps.plot +
  theme(text = element_text(face = "bold", family = "serif"))

Rps.plot

## ----horizontal-plot----------------------------------------------------------
Rps.plot <- Rps.plot +
  coord_flip()

Rps.plot

## ----use-Colors---------------------------------------------------------------
autoplot(Rps.summary, type = "percentage", color = "#18453b") +
  theme_bw() +
  theme(text = element_text(face = "bold", family = "serif"))

## ----sort-axis----------------------------------------------------------------
autoplot(
  Rps.summary,
  type = "percentage",
  color = "#18453b",
  order = "ascending"
) +
  theme_bw() +
  theme(text = element_text(face = "bold", family = "serif"))

