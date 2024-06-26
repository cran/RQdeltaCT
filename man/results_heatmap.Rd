% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RQdeltaCT.R
\name{results_heatmap}
\alias{results_heatmap}
\title{results_heatmap}
\usage{
results_heatmap(
  data,
  sel.Gene = "all",
  dist.row = "euclidean",
  dist.col = "euclidean",
  clust.method = "average",
  col.groups,
  colors = c("navy", "#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
    "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#C32B23", "#A50026", "#8B0000",
    "#7E0202", "#000000"),
  show.colnames = TRUE,
  show.rownames = TRUE,
  border.color = NA,
  fontsize = 10,
  fontsize.col = 10,
  fontsize.row = 10,
  angle.col = 0,
  cellwidth = NA,
  cellheight = NA,
  save.to.tiff = FALSE,
  dpi = 600,
  width = 15,
  height = 15,
  name.tiff = "heatmap_results"
)
}
\arguments{
\item{data}{Object returned from make_Ct_ready() or delta_Ct() functions.}

\item{sel.Gene}{Character vector with names of genes to include, or "all" (default) to use all genes.}

\item{dist.row, dist.col}{Character: name of method used for calculation of distances between rows or columns, derived from stats::dist() function,
must be one of "euclidean" (default) , "maximum", "manhattan", "canberra", "binary" or "minkowski".}

\item{clust.method}{Character: name of used method for agglomeration, derived from stats::hclust() function,
must be one of "ward.D", "ward.D2", "single", "complete", "average" (default), "mcquitty", "median" or "centroid".}

\item{col.groups}{A named list with colors for groups annotation (see example).}

\item{colors}{Vector with colors used to fill created heatmap.}

\item{show.colnames, show.rownames}{Logical: of TRUE, names of columns (sample names) and rows (gene names) will be shown. Both default to TRUE.}

\item{border.color}{Character: color of cell borders on heatmap. If set to NA (default) no border will be drawn.}

\item{fontsize}{Numeric: global fontsize of heatmap. Default to 10.}

\item{fontsize.col, fontsize.row}{Numeric: fontsize of colnames and rownames. Default to 10.}

\item{angle.col}{Integer: angle of the column labels, one of the 0, 45, 90, 270, and 315.}

\item{cellwidth, cellheight}{Numeric: width and height of individual cell. Both default to NA.
These parameters are useful in situations where margins are too small and the plot is cropped (column names and annotation legend are sometimes partially hidden).
Specification of this parameter allows to adjust size of the plot and solve this problem.}

\item{save.to.tiff}{Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.}

\item{dpi}{Integer: resolution of saved .tiff file. Default to 600.}

\item{width}{Numeric: width (in cm) of saved .tiff file. Default to 15.}

\item{height}{Numeric: height (in cm) of saved .tiff file. Default to 15.}

\item{name.tiff}{Character: name of saved .tiff file, without ".tiff" name of extension. Default to "heatmap_results".}
}
\value{
Heatmap with hierarchical clustering, displayed on the graphic device (if save.to.tiff = FALSE)
or saved to .tiff file (if save.to.tiff = TRUE).
}
\description{
This function creatse heatmap with hierarchical clustering.
}
\examples{
library(tidyverse)
library(pheatmap)
data(data.Ct)
data.CtF <- filter_Ct(data.Ct,
                      remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
                       remove.Sample = c("Control08","Control16","Control22"))
data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
# Remember to firstly create named list with colors for groups annotation:
colors.for.groups = list("Group" = c("Disease"="firebrick1","Control"="green3"))
# Vector of colors to fill the heatmap can be also specified to fit the user needings:
colors <- c("navy","navy","#313695","#4575B4","#74ADD1","#ABD9E9",
            "#E0F3F8","#FFFFBF","#FEE090","#FDAE61","#F46D43",
            "#D73027","#C32B23","#A50026","#8B0000",
            "#7E0202","#000000")
results_heatmap(data.dCt,
                sel.Gene = "all",
                col.groups = colors.for.groups,
                colors = colors,
                show.colnames = FALSE,
                show.rownames = TRUE,
                fontsize = 11,
                fontsize.row = 11)

}
