% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RQdeltaCT.R
\name{RQ_dCt}
\alias{RQ_dCt}
\title{RQ_dCt}
\usage{
RQ_dCt(
  data,
  group.study,
  group.ref,
  do.tests = TRUE,
  pairwise = FALSE,
  alternative = "two.sided",
  p.adjust.method = "BH",
  save.to.txt = FALSE,
  name.txt = "results_dCt"
)
}
\arguments{
\item{data}{Data object returned from delta_Ct() or filter_transformed_data() function.}

\item{group.study}{Character: name of study group (group of interest).}

\item{group.ref}{Character: name of reference group.}

\item{do.tests}{Logical: if TRUE, statistical significance of differences between compared groups will be calculated using Student's t test and Mann-Whitney U test. Default to TRUE.}

\item{pairwise}{Logical: if TRUE, a pairwise analysis will be performed (see details). Default to FALSE.}

\item{alternative}{Character: alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".}

\item{p.adjust.method}{Character: p value correction method for multiple testing, one of the "holm", "hochberg", "hommel",
"bonferroni", "BH" (default), "BY","fdr", or "none". See documentation for stats::p.adjust() function for details.}

\item{save.to.txt}{Logical: if TRUE, returned table with results will be saved to .txt file. Default to FALSE.}

\item{name.txt}{Character: name of saved .txt file, without ".txt" name of extension. Default to "results_dCt".}
}
\value{
Data frame with results (if pairwise = FALSE) or, if pairwise = TRUE, the list object
with two elements: a table with the results and the second table with fold change values calculated individually for each sample.
}
\description{
This function performs relative quantification of gene expression using 2^-dCt method.
}
\details{
This function calculates:
1. Means (returned in columns with the "_mean" pattern) and standard deviations (returned in columns with the "_sd" pattern)
of 2^-dCt transformed dCt values for each analyzed gene across compared groups.
2. P values of normality test (Shapiro_Wilk test) performed on 2^-dCt values across compared groups (returned in
columns with the "_norm_p" pattern).
3. Fold change values (returned in "FCh" column) calculated for each gene by dividing  mean of 2^-dCt values in study group
by mean of 2^-dCt values in reference group.
4. Statistics (returned in column with the "_test_stat" pattern) and p values (returned in column with "_test_p" pattern) of
differences in 2^-dCt values between study group and reference group using Student's t test and Mann-Whitney U test.
5. P values adjusted for multiple testing using a selected method.
}
\examples{
library(tidyverse)
library(coin)
data(data.Ct)
data.CtF <- filter_Ct(data.Ct,
                     remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
                     remove.Sample = c("Control08","Control16","Control22"))
data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8", transform = TRUE)
results.dCt <- RQ_dCt(data.dCt,
                           group.study = "Disease",
                           group.ref = "Control")
head(results.dCt)

}
