## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  cache = FALSE
)

## ----eval=FALSE---------------------------------------------------------------
#  remotes::install_github("Donadelnal/RQdeltaCT")
#  library(RQdeltaCT)

## ----echo=FALSE, out.width="650px", dpi = 600, fig.align="center", warning=FALSE, message=FALSE, cache=FALSE----
knitr::include_graphics("figure1ok.png")

## ----message=FALSE, cache=FALSE-----------------------------------------------
# Set path to file:
path <- system.file("extdata",
                    "data_Ct_long.txt",
                    package = "RQdeltaCT")

# Import file using path; remember to specify proper separator, decimal character, and numbers of necessary columns:
library(RQdeltaCT)
library(tidyverse)
data.Ct <- read_Ct_long(path = path,
                        sep = "\t",
                        dec = ".",
                        skip = 0,
                        add.column.Flag = TRUE,
                        column.Sample = 1,
                        column.Gene = 2,
                        column.Ct = 5,
                        column.Group = 9,
                        column.Flag = 4)

## ----cache=FALSE--------------------------------------------------------------
str(data.Ct)

## ----cache=FALSE--------------------------------------------------------------
library(tidyverse)
data.Ct <- mutate(data.Ct,
                  Flag = ifelse(Flag < 1, "Undetermined", "OK"))
str(data.Ct)

## ----cache=FALSE--------------------------------------------------------------
# Set paths to required files: 
path.Ct.file <- system.file("extdata",
                            "data_Ct_wide.txt",
                            package = "RQdeltaCT")
path.design.file <- system.file("extdata",
                                "data_design.txt",
                                package = "RQdeltaCT")

# Import files:
library(tidyverse)
data.Ct <- read_Ct_wide(path.Ct.file = path.Ct.file,
                        path.design.file = path.design.file,
                        sep ="\t",
                        dec = ".")

# Look at the structure:
str(data.Ct)

## ----cache=FALSE--------------------------------------------------------------
# Import file, be aware to specify parameters that fit to imported data:
data.Ct.wide <- read.csv(file = "data/data.Ct.wide.vign.txt",
                         header = TRUE,
                         sep = ",")
str(data.Ct.wide)

# The imported table is now transformed to a long-format structure. The "X" column is unnecessary and is removed. All variables also are converted to a character to unify the class of variables.
library(tidyverse)
data.Ct <- data.Ct.wide %>%
             select(-X) %>%
             mutate(across(everything(), as.character)) %>%
             pivot_longer(cols = -c(Group, Sample), names_to = "Gene", values_to = "Ct")
str(data.Ct)

## ----cache=FALSE--------------------------------------------------------------
data(data.Ct)
str(data.Ct)

data(data.Ct.pairwise)
str(data.Ct.pairwise)

## ----fig.dim=c(7.1,7), cache=FALSE--------------------------------------------
sample.Ct.control <- control_Ct_barplot_sample(data = data.Ct,
                                               flag.Ct = "Undetermined",
                                               maxCt = 35,
                                               flag = c("Undetermined"),
                                               axis.title.size = 9,
                                               axis.text.size = 7,
                                               plot.title.size = 9,
                                               legend.title.size = 9,
                                               legend.text.size = 9)


## ----fig.dim=c(7.1,5.5), cache=FALSE------------------------------------------
gene.Ct.control <- control_Ct_barplot_gene(data = data.Ct,
                                               flag.Ct = "Undetermined",
                                               maxCt = 35,
                                               flag = c("Undetermined"),
                                               axis.title.size = 9,
                                               axis.text.size = 9,
                                               plot.title.size = 9,
                                               legend.title.size = 9,
                                               legend.text.size = 9)

## ----cache=FALSE--------------------------------------------------------------
head(sample.Ct.control[[2]])

## ----cache=FALSE--------------------------------------------------------------
head(gene.Ct.control[[2]])

## ----cache=FALSE--------------------------------------------------------------
# Finding samples with more than half of the unreliable Ct values.
low.quality.samples <- filter(sample.Ct.control[[2]], Not.reliable.fraction > 0.5)$Sample
low.quality.samples <- as.vector(low.quality.samples)                        
low.quality.samples

## -----------------------------------------------------------------------------
# Finding genes with more than half of the unreliable Ct values in given group.
low.quality.genes <- filter(gene.Ct.control[[2]], Not.reliable.fraction > 0.5)$Gene
low.quality.genes <- unique(as.vector(low.quality.genes))                        
low.quality.genes

## ----cache=FALSE--------------------------------------------------------------
# Objects returned from the `low_quality_samples()` and `low_quality_genes()`functions can be used directly:
data.CtF <- filter_Ct(data = data.Ct,
                      flag.Ct = "Undetermined",
                      maxCt = 35,
                      flag = c("Undetermined"),
                      remove.Gene = low.quality.genes,
                      remove.Sample = low.quality.samples)

# Check dimensions of data before and after filtering:
dim(data.Ct)
dim(data.CtF)

## ----cache=FALSE--------------------------------------------------------------
# Without imputation:
data.CtF.ready <- make_Ct_ready(data = data.CtF,
                                imput.by.mean.within.groups = FALSE)
# A part of the data with missing values:
as.data.frame(data.CtF.ready)[19:25,]

# With imputation:
data.CtF.ready <- make_Ct_ready(data = data.CtF,
                                imput.by.mean.within.groups = TRUE)
# Missing values were imputed:
as.data.frame(data.CtF.ready)[19:25,]

## ----fig.dim=c(7.1,4), cache=FALSE--------------------------------------------
library(ctrlGene)
# Remember that the number of colors in col parameter should be equal to the number of tested genes:
ref <- find_ref_gene(data = data.CtF.ready,
                     groups = c("Disease","Control"),
                     candidates = c("Gene4", "Gene8","Gene10","Gene16","Gene17", "Gene18"),
                     col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "#1F77B4", "black"),
                     angle = 60,
                     axis.text.size = 7,
                     norm.finder.score = TRUE,
                     genorm.score = TRUE)
ref[[2]]

## ----cache=FALSE--------------------------------------------------------------
# For 2^-dCt^ method:
data.dCt.exp <- delta_Ct(data = data.CtF.ready,
                     normalise = TRUE,
                     ref = "Gene8",
                     transform = TRUE)

## ----cache=FALSE--------------------------------------------------------------
# For 2^-ddCt^ method:
data.dCt <- delta_Ct(data = data.CtF.ready,
                     normalise = TRUE,
                     ref = "Gene8",
                     transform = FALSE)

## ----fig.dim=c(7.1,6), cache=FALSE--------------------------------------------
control_boxplot_sample <- control_boxplot_sample(data = data.dCt,
                                                 axis.text.size = 7)

## ----fig.dim=c(7.1,4)---------------------------------------------------------
control_boxplot_gene <- control_boxplot_gene(data = data.dCt,
                                                 by.group = TRUE,
                                                 axis.text.size = 10)

## ----fig.dim=c(7.1,4), cache=FALSE--------------------------------------------
control_cluster_sample(data = data.dCt,
                       method.dist = "euclidean",
                       method.clust = "average",
                       label.size = 0.6)
control_cluster_gene(data = data.dCt,
                     method.dist = "euclidean",
                     method.clust = "average",
                     label.size = 0.8)

## ----fig.dim=c(5,5.5), fig.align='center', cache=FALSE------------------------
control.pca.sample <- control_pca_sample(data = data.dCt,
                                         point.size = 3,
                                         label.size = 2.5,
                                         legend.position = "top")

## ----fig.dim=c(5,5), fig.align='center', cache=FALSE--------------------------
control.pca.gene <- control_pca_gene(data = data.dCt)

## ----cache=FALSE--------------------------------------------------------------
data.dCtF <- filter_transformed_data(data = data.dCt,
                                     remove.Sample = c("Control11"))

## ----cache=FALSE--------------------------------------------------------------
data.dCt.exp <- delta_Ct(data = data.CtF.ready,
                         ref = "Gene8",
                         transform = TRUE)
library(coin)
results.dCt <- RQ_dCt(data = data.dCt.exp,
                            do.tests = TRUE,
                            group.study = "Disease",
                            group.ref = "Control")

# Obtained table can be sorted by, e.g. p values from the Mann-Whitney U test:
head(as.data.frame(arrange(results.dCt, MW_test_p)))

## ----cache=FALSE--------------------------------------------------------------
data.dCt <- delta_Ct(data = data.CtF.ready,
                     ref = "Gene8",
                     transform = FALSE)
library(coin)
results.ddCt <- RQ_ddCt(data = data.dCt,
                   group.study = "Disease",
                   group.ref = "Control",
                   do.tests = TRUE)

# Obtained table can be sorted by, e.g. p values from the Mann-Whitney U test:
head(as.data.frame(arrange(results.ddCt, MW_test_p)))

## ----fig.dim=c(7.1,4), cache=FALSE--------------------------------------------
# Variant with p values depending on the normality of the data:
library(ggsignif)
# Specifying vector with significance labels: 
signif.labels <- c("****","**","ns."," ns. ","  ns.  ","   ns.   ","    ns.    ","     ns.     ","      ns.      ","       ns.       ","        ns.        ","         ns.         ","          ns.          ","***")
# Genes with p < 0.05 and 2-fold changed expression between compared groups are considered significant: 
FCh.plot <- FCh_plot(data = results.ddCt,
                   use.p = TRUE,
                   mode = "depends",
                   p.threshold = 0.05,
                   use.FCh = TRUE,
                   FCh.threshold = 2,
                   signif.show = TRUE,
                   signif.labels = signif.labels,
                   angle = 20)
# Access the table with results:
head(as.data.frame(FCh.plot[[2]]))

## ----fig.dim=c(7.1,4), cache=FALSE--------------------------------------------
user <- data.dCt %>%
          pivot_longer(cols = -c(Group, Sample),
                       names_to = "Gene",
                       values_to = "dCt") %>%
          group_by(Gene) %>%
          summarise(MW_test_p = wilcox.test(dCt ~ Group)$p.value,
                    .groups = "keep")
# The stats::wilcox.test() functions is limited to cases without ties; therefore, a warning "cannot compute exact p-value with ties" will appear when ties occur.

FCh.plot <- FCh_plot(data = results.ddCt,
                   use.p = TRUE,
                   mode = "user",
                   p.threshold = 0.05,
                   use.FCh = TRUE,
                   FCh.threshold = 2,
                   signif.show = TRUE,
                   signif.labels = signif.labels,
                   angle = 30)
# Access the table with results:
head(as.data.frame(FCh.plot[[2]]))

## ----fig.dim=c(5,5.5), fig.align='center', cache=FALSE------------------------
# Genes with p < 0.05 and 2-fold changed expression between compared groups are considered significant: 
volcano <- results_volcano(data = results.ddCt,
                         mode = "depends",
                         p.threshold = 0.05,
                         FCh.threshold = 2)
# Access the table with results:
head(as.data.frame(volcano[[2]]))

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final_boxplot <- results_boxplot(data = data.dCtF,
                                 sel.Gene = c("Gene1","Gene12", "Gene16","Gene19"),
                                 by.group = TRUE,
                                 signif.show = TRUE,
                                 signif.labels = c("****","*","***"," * "),
                                 signif.dist = 1.05,
                                 faceting = TRUE,
                                 facet.row = 1,
                                 facet.col = 4,
                                 y.exp.up = 0.1,
                                 angle = 20,
                                 y.axis.title = "dCt")

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final_barplot <- results_barplot(data = data.dCtF,
                                 sel.Gene = c("Gene1","Gene12","Gene19","Gene20"),
                                 signif.show = TRUE,
                                 signif.labels = c("****","*","***"," * "),
                                 angle = 30,
                                 signif.dist = 1.05,
                                 faceting = TRUE,
                                 facet.row = 1,
                                 facet.col = 4,
                                 y.exp.up = 0.1,
                                 y.axis.title = "dCt")

## ----fig.dim=c(7.1,5), cache=FALSE--------------------------------------------
# Create named list with colors for groups annotation:
colors.for.groups = list("Group" = c("Disease"="firebrick1","Control"="green3"))
# Vector of colors for heatmap can be also specified to fit the user needings:
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
                fontsize.row = 11,
                cellwidth = 4)
# Cellwidth parameter was set to 4 to avoid cropping the image on the right side.

## ----fig.dim=c(5,5.5), fig.align='center', cache=FALSE------------------------
pca.kmeans <- pca_kmeans(data.dCt, 
                           sel.Gene = c("Gene1","Gene16","Gene19","Gene20"), 
                           legend.position = "top")
# Access to the confusion matrix:
pca.kmeans[[2]]

## ----fig.dim=c(5,6), fig.align='center', cache=FALSE--------------------------
pca.kmeans[[1]] + theme(legend.box = "vertical")

## ----fig.dim=c(6.5,6.5), fig.align='center', cache=FALSE----------------------
library(Hmisc)
library(corrplot)
# To make the plot more readable, only part of the data was used:
corr.samples <- corr_sample(data = data.dCt[15:30, ],
                            method = "pearson",
                            order = "hclust",
                            size = 0.7,
                            p.adjust.method = "BH")

## ----fig.dim=c(6.5,6.5), fig.align='center', cache=FALSE----------------------
library(Hmisc)
library(corrplot)
corr.genes <- corr_gene(data = data.dCt,
                            method = "spearman",
                            order = "FPC",
                            size = 0.7,
                            p.adjust.method = "BH")

## ----fig.dim=c(4.5,4.5), fig.align='center', cache=FALSE----------------------
library(ggpmisc)
Disease6_Control17 <- single_pair_sample(data = data.dCt,
                                         x = "Disease6",
                                         y = "Control17",
                                         point.size = 3,
                                         labels = TRUE,
                                         label = c("eq", "R2", "p"),
                                         label.position.x = 0.05)

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
library(ggpmisc)
Gene16_Gene17 <- single_pair_gene(data.dCt,
                                    x = "Gene16",
                                    y = "Gene17",
                                    by.group = TRUE,
                                    point.size = 3,
                                    labels = TRUE,
                                    label = c("eq", "R2", "p"),
                                    label.position.x = c(0.05),
                                    label.position.y = c(1,0.95))

## ----cache=FALSE--------------------------------------------------------------
library(pROC)
# Remember to specify the numbers of rows (panels.row parameter) and columns (panels.col parameter) to be sufficient to arrange panels: 
roc_parameters <- ROCh(data = data.dCt,
                       sel.Gene = c("Gene1","Gene12","Gene16","Gene19"),
                       groups = c("Disease","Control"),
                       panels.row = 2,
                       panels.col = 2)
roc_parameters

## ----echo=FALSE, out.width="500px", fig.align="center", warning=FALSE, message=FALSE, cache=FALSE----
knitr::include_graphics("ROC_plot_ind.png")

## ----fig.dim=c(5,4), fig.align='center', cache=FALSE--------------------------
library(oddsratio)

# Remember to set the increment parameter.
log.reg.results <- log_reg(data = data.dCt,
                           increment = 1,
                           sel.Gene = c("Gene1","Gene12","Gene16","Gene19"),
                           group.study = "Disease",
                           group.ref = "Control")
log.reg.results[[2]]

## ----cache=FALSE--------------------------------------------------------------
data(data.Ct.pairwise)
str(data.Ct.pairwise)

## ----fig.dim=c(7.1,5), cache=FALSE--------------------------------------------
library(tidyverse)
sample.Ct.control.pairwise <- control_Ct_barplot_sample(data = data.Ct.pairwise,
                                               flag.Ct = "Undetermined",
                                               maxCt = 35,
                                               flag = c("Undetermined"),
                                               axis.title.size = 9,
                                               axis.text.size = 9,
                                               plot.title.size = 9,
                                               legend.title.size = 9,
                                               legend.text.size = 9)


## ----fig.dim=c(7.1,5.5), cache=FALSE------------------------------------------
gene.Ct.control.pairwise <- control_Ct_barplot_gene(data = data.Ct.pairwise,
                                               flag.Ct = "Undetermined",
                                               maxCt = 35,
                                               flag = c("Undetermined"),
                                               axis.title.size = 9,
                                               axis.text.size = 9,
                                               plot.title.size = 9,
                                               legend.title.size = 9,
                                               legend.text.size = 9)

## ----cache=FALSE--------------------------------------------------------------
head(sample.Ct.control.pairwise[[2]])

## ----cache=FALSE--------------------------------------------------------------
head(gene.Ct.control.pairwise[[2]], 10)

## ----cache=FALSE--------------------------------------------------------------
# Finding samples with more than half of the unreliable Ct values.
low.quality.samples.pairwise <- filter(sample.Ct.control.pairwise[[2]],
                                       Not.reliable.fraction > 0.5)$Sample
low.quality.samples.pairwise <- as.vector(low.quality.samples.pairwise)                        
low.quality.samples.pairwise

## -----------------------------------------------------------------------------
# Finding genes with more than half of the unreliable Ct values in at least one group.
low.quality.genes.pairwise <- filter(gene.Ct.control.pairwise[[2]],
                                     Not.reliable.fraction > 0.5)$Gene
low.quality.genes.pairwise <- unique(as.vector(low.quality.genes.pairwise))                        
low.quality.genes.pairwise

## ----cache=FALSE--------------------------------------------------------------
# Objects returned from the `low_quality_samples()` and 
# `low_quality_genes()`functions can be used directly:
data.Ct.pairwiseF <- filter_Ct(data = data.Ct.pairwise,
                      flag.Ct = "Undetermined",
                      maxCt = 35,
                      flag = c("Undetermined"),
                      remove.Gene = low.quality.genes.pairwise,
                      remove.Sample = low.quality.samples.pairwise)

# Check dimensions of data before and after filtering:
dim(data.Ct.pairwise)
dim(data.Ct.pairwiseF)

## ----cache=FALSE--------------------------------------------------------------
# Without imputation:
data.Ct.pairwiseF.ready <- make_Ct_ready(data = data.Ct.pairwiseF,
                                         imput.by.mean.within.groups = FALSE)
# A part of the data with missing values:
as.data.frame(data.Ct.pairwiseF.ready)[9:19,10:15]

# With imputation:
data.Ct.pairwiseF.ready <- make_Ct_ready(data = data.Ct.pairwiseF,
                                         imput.by.mean.within.groups = TRUE)
# Missing values were imputed:
as.data.frame(data.Ct.pairwiseF.ready)[9:19,10:15]

## ----fig.dim=c(7.1,4), cache=FALSE--------------------------------------------
library(ctrlGene)
# Remember that the number of colors in col parameter should be equal to the number of tested genes:
ref.pairwise <- find_ref_gene(data = data.Ct.pairwiseF.ready,
                     groups = c("After","Before"),
                     candidates = c("Gene4","Gene13","Gene20"),
                     col = c("#66c2a5", "#fc8d62","#6A6599"),
                     angle = 90,
                     axis.text.size = 7,
                     norm.finder.score = TRUE,
                     genorm.score = TRUE)
ref.pairwise[[2]]

## ----cache=FALSE--------------------------------------------------------------
data.dCt.pairwise <- delta_Ct(data = data.Ct.pairwiseF.ready, ref = "Gene4", transform = FALSE)

## ----cache=FALSE--------------------------------------------------------------
data.dCt.exp.pairwise <- delta_Ct(data = data.Ct.pairwiseF.ready,
                     ref = "Gene4",
                     transform = TRUE)

## ----cache=FALSE--------------------------------------------------------------
data.dCt.pairwise <- delta_Ct(data = data.Ct.pairwiseF.ready, ref = "Gene4", transform = FALSE)
library(coin)
results.dCt.pairwise <- RQ_dCt(data = data.dCt.pairwise,
                               do.tests = TRUE,
                               pairwise = TRUE,
                               group.study = "After",
                               group.ref = "Before")

# Obtained table can be sorted by, e.g. p values from the Mann-Whitney U test:
results <- as.data.frame(arrange(results.dCt.pairwise[[1]], MW_test_p))
head(results)
# Access to the table with fold change values calculated individually for each pair of sampleS:
FCh <- results.dCt.pairwise[[2]]
head(FCh)

## ----cache=FALSE--------------------------------------------------------------
data.dCt.pairwise <- delta_Ct(data = data.Ct.pairwiseF.ready, ref = "Gene4", transform = FALSE)
library(coin)
# Remember to set pairwise = TRUE:
results.ddCt.pairwise <- RQ_ddCt(data = data.dCt.pairwise,
                   group.study = "After",
                   group.ref = "Before",
                   pairwise = TRUE,
                   do.tests = TRUE)

# Obtained table can be sorted by, e.g. p values from the Mann-Whitney U test:
results <- as.data.frame(arrange(results.ddCt.pairwise[[1]], MW_test_p))
head(results)
# Access to the table with fold change values calculated individually for each pair of samples:
FCh <- results.ddCt.pairwise[[2]]
head(FCh)

## ----fig.dim=c(7.1,6), cache=FALSE--------------------------------------------
control_boxplot_sample <- control_boxplot_sample(data = data.dCt.pairwise,
                                                 axis.text.size = 9,
                                                 y.axis.title = "dCt")

## ----fig.dim=c(7.1,4)---------------------------------------------------------
control_boxplot_gene <- control_boxplot_gene(data = data.dCt.pairwise,
                                             by.group = TRUE,
                                             axis.text.size = 10,
                                             y.axis.title = "dCt")

## ----fig.dim=c(7.1,6), cache=FALSE--------------------------------------------
# Remember to set pairwise.FCh to TRUE:
FCh <- results.dCt.pairwise[[2]]
control.boxplot.sample.pairwise <- control_boxplot_sample(data = FCh,
                                                          pairwise.FCh = TRUE,
                                                          axis.text.size = 9,
                                                          y.axis.title = "Fold change")
# There are some very high values, we can identify them using:
head(arrange(FCh, -FCh))

## ----fig.dim=c(7.1,4)---------------------------------------------------------
control.boxplot.gene.pairwise <- control_boxplot_gene(data = FCh,
                                                 by.group = FALSE,
                                                 pairwise.FCh = TRUE,
                                                 axis.text.size = 10,
                                                 y.axis.title = "Fold change")

## ----fig.dim=c(7.1,4), cache=FALSE--------------------------------------------
# Hierarchical clustering of samples:
control_cluster_sample(data = data.dCt.pairwise,
                       method.dist = "euclidean",
                       method.clust = "average",
                       label.size = 0.6)
# Hierarchical clustering of genes:
control_cluster_gene(data = data.dCt.pairwise,
                     method.dist = "euclidean",
                     method.clust = "average",
                     label.size = 0.8)

## ----fig.dim=c(7.1,4), cache=FALSE--------------------------------------------
# Remember to set pairwise.FCh = TRUE:
control_cluster_sample(data = FCh,
                       pairwise.FCh = TRUE,
                       method.dist = "euclidean",
                       method.clust = "average",
                       label.size = 0.7)
control_cluster_gene(data = FCh,
                     pairwise.FCh = TRUE,
                     method.dist = "euclidean",
                     method.clust = "average",
                     label.size = 0.8)

## ----fig.dim=c(5,5.5), fig.align='center', cache=FALSE------------------------
control.pca.sample.pairwise <- control_pca_sample(data = data.dCt.pairwise,
                                         point.size = 3,
                                         label.size = 2.5,
                                         hjust = 0.5,
                                         legend.position = "top")

## ----fig.dim=c(5,5), fig.align='center', cache=FALSE--------------------------
control.pca.gene.pairwise <- control_pca_gene(data = data.dCt.pairwise,
                                              hjust = 0.5)

## ----fig.dim=c(5,5.5), fig.align='center', cache=FALSE------------------------
control.pca.sample.pairwise <- control_pca_sample(data = FCh,
                                         pairwise.FCh = TRUE,
                                         colors = "black",
                                         point.size = 3,
                                         label.size = 2.5,
                                         hjust = 0.5)

## ----fig.dim=c(5,5), fig.align='center', cache=FALSE--------------------------
control.pca.gene.pairwise <- control_pca_gene(data = FCh,
                                     pairwise.FCh = TRUE,
                                     color = "black",
                                     hjust = 0.5)

## ----cache=FALSE--------------------------------------------------------------
data.dCt.pairwise.F <- filter_transformed_data(data = data.dCt.pairwise,
                                     remove.Sample = c("Sample22", "Sample23", "Sample15","Sample03"))
dim(data.dCt.pairwise)
dim(data.dCt.pairwise.F)

## -----------------------------------------------------------------------------
# Remember to set pairwise = TRUE:
results.ddCt.pairwise <- RQ_ddCt(data = data.dCt.pairwise.F,
                   group.study = "After",
                   group.ref = "Before",
                   pairwise = TRUE,
                   do.tests = TRUE)

# Obtained table can be sorted by, e.g. p values from the Mann-Whitney U test:
results <- as.data.frame(arrange(results.ddCt.pairwise[[1]], MW_test_p))
head(results)

## ----fig.dim=c(7.1,4), cache=FALSE--------------------------------------------
library(ggsignif)
#  Remember to use the first element of list object returned by `RQ_dCt()` or RQ_ddCt() function:
FCh.plot.pairwise <- FCh_plot(data = results.ddCt.pairwise[[1]],
                   use.p = TRUE,
                   mode = "depends.adj",
                   p.threshold = 0.05,
                   use.FCh = TRUE,
                   FCh.threshold = 1.5,
                   angle = 20)
# Access the table with results:
head(as.data.frame(FCh.plot.pairwise[[2]]))

## ----fig.dim=c(7.1,4), cache=FALSE--------------------------------------------
# Firstly prepare a vector with significance labels specified according to the user needings:

signif.labels <- c("ns.",
                   " ns. ",
                   "  ns.  ",
                   "   ns.   ",
                   "    ns.    ",
                   "     ns.     ",
                   "      ns.      ",
                   "       ns.       ",
                   "        ns.        ",
                   "         ns.         ",
                   "**",
                   "***")
# Remember to set signif.show = TRUE:
FCh.plot.pairwise <- FCh_plot(data = results.ddCt.pairwise[[1]],
                   use.p = TRUE,
                   mode = "depends.adj",
                   p.threshold = 0.05,
                   use.FCh = TRUE,
                   FCh.threshold = 1.5,
                   use.sd = TRUE,
                   signif.show = TRUE,
                   signif.labels = signif.labels,
                   signif.dist = 0.2,
                   angle = 20)

## ----fig.dim=c(7.1,4), cache=FALSE--------------------------------------------
# Variant with user p values - the used p values are calculated using the stats::wilcox.test() function:
user <- data.dCt.pairwise %>%
          pivot_longer(cols = -c(Group, Sample),
                       names_to = "Gene",
                       values_to = "dCt") %>%
          group_by(Gene) %>%
          summarise(MW_test_p = wilcox.test(dCt ~ Group)$p.value,
                    .groups = "keep")
# The stats::wilcox.test() functions is limited to cases without ties; therefore, 
# a warning "cannot compute exact p-value with ties" will appear when ties occur.
signif.labels <- c("ns.",
                   " ns. ",
                   "  ns.  ",
                   "   ns.   ",
                   "    ns.    ",
                   "     ns.     ",
                   "      ns.      ",
                   "       ns.       ",
                   "        ns.        ",
                   "         ns.         ",
                   " * ",
                   "***")
FCh.plot <- FCh_plot(data = results.ddCt.pairwise[[1]],
                   use.p = TRUE,
                   mode = "user",
                   p.threshold = 0.05,
                   use.FCh = TRUE,
                   FCh.threshold = 1.5,
                   use.sd = TRUE,
                   signif.show = TRUE,
                   signif.labels = signif.labels,
                   signif.dist = 0.4,
                   angle = 30)

## ----fig.dim=c(6,4.5), fig.align='center', cache=FALSE------------------------
# Genes with p < 0.05 and 2-fold changed expression between compared groups 
# are considered significant. Remember to use the first element of list object 
# returned by RQ_ddCt() function: 
RQ.volcano.pairwise <- results_volcano(data = results.ddCt.pairwise[[1]],
                         mode = "depends.adj",
                         p.threshold = 0.05,
                         FCh.threshold = 1.5)
# Access the table with results:
head(as.data.frame(RQ.volcano.pairwise[[2]]))

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final.boxplot.pairwise <- results_boxplot(data = data.dCt.pairwise.F,
                                 sel.Gene = c("Gene8","Gene19"),
                                 by.group = TRUE,
                                 signif.show = TRUE,
                                 signif.labels = c("***","*"),
                                 signif.dist = 1.2,
                                 faceting = TRUE,
                                 facet.row = 1,
                                 facet.col = 2,
                                 y.exp.up = 0.1,
                                 y.axis.title = "dCt")

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
data.dCt.pairwise.F$Group <- factor(data.dCt.pairwise.F$Group, levels = c("Before", "After"))

final.boxplot.pairwise <- results_boxplot(data = data.dCt.pairwise.F,
                                 sel.Gene = c("Gene8","Gene19"),
                                 by.group = TRUE,
                                 signif.show = TRUE,
                                 signif.labels = c("***","*"),
                                 signif.dist = 1.2,
                                 faceting = TRUE,
                                 facet.row = 1,
                                 facet.col = 2,
                                 y.exp.up = 0.1,
                                 y.axis.title = "dCt")

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final.barplot.pairwise <- results_barplot(data = data.dCt.pairwise.F,
                                 sel.Gene = c("Gene8","Gene19"),
                                 signif.show = TRUE,
                                 signif.labels = c("***","*"),
                                 angle = 30,
                                 signif.dist = 1.05,
                                 faceting = TRUE,
                                 facet.row = 1,
                                 facet.col = 4,
                                 y.exp.up = 0.1,
                                 y.axis.title = "dCt")

## ----fig.dim=c(7.1,5), cache=FALSE--------------------------------------------
# Create named list with colors for groups annotation:
colors.for.groups = list("Group" = c("After"="firebrick1", "Before"="green3"))
# Vector of colors for heatmap can be also specified to fit the user needings:
colors <- c("navy","navy","#313695","#313695","#4575B4","#4575B4","#74ADD1","#74ADD1",
            "#ABD9E9","#E0F3F8","#FFFFBF","#FEE090","#FDAE61","#F46D43",
            "#D73027","#C32B23","#A50026","#8B0000", 
            "#7E0202","#000000")
library(pheatmap)
results_heatmap(data.dCt.pairwise.F,
                sel.Gene = "all",
                col.groups = colors.for.groups,
                colors = colors,
                show.colnames = TRUE,
                show.rownames = TRUE,
                fontsize = 11,
                fontsize.row = 10,
                cellwidth = 8,
                angle.col = 90)
# Cellwidth parameter was set to 4 to avoid cropping the image on the right side.

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
parallel.plot <- parallel_plot(data = data.dCt.pairwise.F,
                                sel.Gene = c("Gene19","Gene8"),
                               order = c(4,3))

## ----fig.dim=c(5,5.5), fig.align='center', cache=FALSE------------------------
pca.kmeans <- pca_kmeans(data.dCt.pairwise.F, 
                           sel.Gene = c("Gene8","Gene19"), 
                           legend.position = "top")

## ----fig.dim=c(5,6), fig.align='center', cache=FALSE--------------------------
pca.kmeans[[1]] + theme(legend.box = "vertical")

## -----------------------------------------------------------------------------
pca.kmeans[[2]]

## ----fig.dim=c(6.5,6.5), fig.align='center', cache=FALSE----------------------
library(Hmisc)
library(corrplot)
# To make the plot more readable, only part of the data was used:
corr.samples <- corr_sample(data = data.dCt.pairwise.F[1:10, ],
                            method = "pearson",
                            order = "hclust",
                            size = 0.7,
                            p.adjust.method = "BH",
                            add.coef = "white")

## ----fig.dim=c(6.5,6.5), fig.align='center', cache=FALSE----------------------
library(Hmisc)
library(corrplot)
corr.genes <- corr_gene(data = data.dCt.pairwise.F,
                            method = "pearson",
                            order = "FPC",
                            size = 0.7,
                            p.adjust.method = "BH")

## ----fig.dim=c(7.1,5), fig.align='center', cache=FALSE------------------------
library(ggpmisc)
Sample08_Sample26 <- single_pair_sample(data = data.dCt.pairwise,
                                        pairwise.data = TRUE,
                                        by.group = TRUE,
                                         x = "Sample08",
                                         y = "Sample26",
                                         point.size = 3,
                                         labels = TRUE,
                                         label = c("eq", "R2", "p"),
                                         label.position.x = 0.05,
                                        label.position.y = c(1, 0.95))

## ----fig.dim=c(7.1,5), fig.align='center', cache=FALSE------------------------
library(ggpmisc)
Gene16_Gene17 <- single_pair_gene(data.dCt.pairwise.F,
                                    x = "Gene16",
                                    y = "Gene17",
                                    by.group = TRUE,
                                    point.size = 3,
                                    labels = TRUE,
                                    label = c("eq", "R2", "p"),
                                    label.position.x = c(0.05),
                                    label.position.y = c(1,0.95))

## ----cache=FALSE--------------------------------------------------------------
library(pROC)
# Remember to specify the numbers of rows (panels.row parameter) and columns (panels.col parameter) to provide sufficient place to arrange panels: 
roc_parameters <- ROCh(data = data.dCt.pairwise.F,
                       sel.Gene = c("Gene8","Gene19"),
                       groups = c("After","Before"),
                       panels.row = 1,
                       panels.col = 2)
# Access to calculated parameters:
roc_parameters

## ----echo=FALSE, out.width="500px", fig.align="center", warning=FALSE, message=FALSE, cache=FALSE----
knitr::include_graphics("ROC_plot.png")

## -----------------------------------------------------------------------------
# Filter data:
data <- data.dCt.pairwise.F[, colnames(data.dCt.pairwise.F) %in% c("Group", "Sample", "Gene19")]
# Perform analysis:
data_roc <- roc(response = data$Group,
            predictor = as.data.frame(data)$Gene19,
            levels = c("Before","After"),
            smooth = FALSE,
            auc = TRUE,
            plot = FALSE,
            ci = TRUE,
            of = "auc",
            quiet = TRUE)
# Gain parameters:
parameters <- coords(data_roc,
                    "best",
                    ret = c("threshold",
                            "specificity",
                            "sensitivity",
                            "accuracy",
                            "ppv",
                            "npv",
                            "youden"))
parameters
# Gain AUC
data_roc$auc

## ----fig.dim=c(5,4), fig.align='center', cache=FALSE--------------------------
library(oddsratio)

# Remember to set the increment parameter.
log.reg.results <- log_reg(data = data.dCt.pairwise.F,
                           increment = 1,
                           sel.Gene = c("Gene8","Gene19"),
                           group.study = "After",
                           group.ref = "Before",
                           log.axis = TRUE)
log.reg.results[[2]]

## -----------------------------------------------------------------------------
sessionInfo()

