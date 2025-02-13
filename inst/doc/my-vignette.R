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

# Import file using path; remember to specify proper separator, decimal character, and number of necessary columns:
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
# Import file, be aware of specifying parameters that fit the imported data:
data.Ct.wide <- read.csv(file = "data/data.Ct.wide.vign.txt",
                         header = TRUE,
                         sep = ",")
str(data.Ct.wide)

# The imported table is now transformed into a long-format structure.  
library(tidyverse)
data.Ct <- data.Ct.wide %>%
             select(-X) %>% # The "X" column is unnecessary and is removed.
             mutate(across(everything(), as.character)) %>% # All variables also are converted to a character to unify the class of variables.
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

## ----fig.dim=c(7.1,8), cache=FALSE--------------------------------------------
library(tidyverse)
library(pheatmap)
data(data.Ct)

# Vector of colors to fill the heatmap can be specified to fit the user's needs:
colors <- c("#4575B4","#FFFFBF","#C32B23")
control_heatmap(data.Ct,
                sel.Gene = "all",
                colors = colors,
                show.colnames = TRUE,
                show.rownames = TRUE, 
                fontsize = 9,
                fontsize.row = 9,
                angle.col = 45)

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
                     groups = c("AAA","Control"),
                     candidates = c("CCL5", "IL1B","GAPDH","TGFB","TNF", "VEGFA"),
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
                         ref = "GAPDH",
                         transform = TRUE)

## ----cache=FALSE--------------------------------------------------------------
# For 2^-ddCt^ method:
data.dCt <- delta_Ct(data = data.CtF.ready,
                     normalise = TRUE,
                     ref = "GAPDH",
                     transform = FALSE)

## ----fig.dim=c(7.1,6), cache=FALSE--------------------------------------------
control_boxplot_sample <- control_boxplot_sample(data = data.dCt,
                                                 y.axis.title = "dCt",
                                                 axis.text.size = 7)

## ----fig.dim=c(7.1,4)---------------------------------------------------------
control_boxplot_gene <- control_boxplot_gene(data = data.dCt,
                                             by.group = TRUE,
                                             y.axis.title = "dCt",
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

## ----fig.dim=c(4,4.5), fig.align='center', cache=FALSE------------------------
control.pca.sample <- control_pca_sample(data = data.dCt,
                                         point.size = 3,
                                         label.size = 2.5,
                                         legend.position = "top")

## ----fig.dim=c(4,4), fig.align='center', cache=FALSE--------------------------
control.pca.gene <- control_pca_gene(data = data.dCt)

## ----cache=FALSE--------------------------------------------------------------
data.dCtF <- filter_transformed_data(data = data.dCt,
                                     remove.Sample = c("Control11"))

## ----cache=FALSE--------------------------------------------------------------
data.dCt.exp <- delta_Ct(data = data.CtF.ready,
                         ref = "GAPDH",
                         transform = TRUE)
library(coin)
results.dCt <- RQ_dCt(data = data.dCt.exp,
                      do.tests = TRUE,
                      group.study = "AAA",
                      group.ref = "Control")

# Obtained table can be sorted by, e.g. p values from the Mann-Whitney U test:
head(as.data.frame(arrange(results.dCt, MW_test_p)))

## ----cache=FALSE--------------------------------------------------------------
data.dCt <- delta_Ct(data = data.CtF.ready,
                     ref = "GAPDH",
                     transform = FALSE)
library(coin)
results.ddCt <- RQ_ddCt(data = data.dCt,
                        group.study = "AAA",
                        group.ref = "Control",
                        do.tests = TRUE)

# Obtained table can be sorted by, e.g. p values from the Mann-Whitney U test:
head(as.data.frame(arrange(results.ddCt, MW_test_p)))

## ----fig.dim=c(7.1,4), cache=FALSE--------------------------------------------
# Variant with p values depending on the normality of the data:
library(ggsignif)
# Specifying vector with significance labels: 
signif.labels <- c("****",
                   "**",
                   "ns.",
                   " ns. ",
                   "  ns.  ",
                   "   ns.   ",
                   "    ns.    ",
                   "     ns.     ",
                   "      ns.      ",
                   "       ns.       ",
                   "        ns.        ",
                   "         ns.         ",
                   "          ns.          ",
                   "***")
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
# The stats::wilcox.test() functions is limited to cases without ties; 
# therefore, a warning "cannot compute exact p-value with ties" will appear when ties occur.

FCh.plot <- FCh_plot(data = results.ddCt,
                   use.p = TRUE,
                   mode = "user",
                   p.threshold = 0.05,
                   use.FCh = TRUE,
                   FCh.threshold = 2,
                   signif.show = TRUE,
                   signif.labels = signif.labels,
                   angle = 30)
# Access the table with results (p.used column was changed):
head(as.data.frame(FCh.plot[[2]]))

## ----fig.dim=c(4,4.5), fig.align='center', cache=FALSE------------------------
# Genes with p < 0.05 and 2-fold changed expression between compared groups are considered significant: 
volcano <- results_volcano(data = results.ddCt,
                           mode = "depends",
                           p.threshold = 0.05,
                           FCh.threshold = 2)
# Access the table with results:
head(as.data.frame(volcano[[2]]))

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final_boxplot <- results_boxplot(data = data.dCtF,
                                 sel.Gene = c("ANGPT1","IL8", "VEGFB"),
                                 by.group = TRUE,
                                 signif.show = TRUE,
                                 signif.labels = c("****","**","***"),
                                 signif.dist = 1.05,
                                 faceting = TRUE,
                                 facet.row = 1,
                                 facet.col = 4,
                                 y.exp.up = 0.1,
                                 angle = 20,
                                 y.axis.title = "dCt")

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final_boxplot_no_fac <- results_boxplot(data = data.dCtF,
                                        sel.Gene = c("ANGPT1","IL8", "VEGFB"),
                                        by.group = TRUE,
                                        signif.show = FALSE,  # Disable significance labels
                                        faceting = FALSE, # Disable faceting
                                        y.axis.title = "dCt") +
                                        theme(axis.text.x = element_text(size = 5, colour = "black"))

# Add x axis annotations and ticks:
final_boxplot_no_fac <- final_boxplot_no_fac + 
                         theme(axis.text.x = element_text(size = 11, colour = "black", face="italic"), # Use italic font for human gene symbols
                               axis.ticks.x = element_line(colour = "black"))
final_boxplot_no_fac

## ----cache=FALSE--------------------------------------------------------------
data.label <- data.frame(matrix(nrow = 3, ncol = 4)) # Number of rows should be equal to number of genes 
rownames(data.label) <- c("ANGPT1","IL8","VEGFB")
colnames(data.label) <- c("x", "xend", "y", "annotation")
data.label$Gene <- rownames(data.label)

data.label$y <- 1 + c(max(data.dCtF$ANGPT1), max(data.dCtF$IL8), max(data.dCtF$VEGFB))
data.label$x <- c(0.81,1.81,2.81)
data.label$xend <- c(1.19,2.19,3.19)
data.label$annotation <- c("****","**","***")

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final_boxplot_no_fac_ok <- final_boxplot_no_fac +
                            geom_signif(annotation = data.label$annotation, 
                                        y_position = data.label$y, 
                                        xmin = data.label$x, 
                                        xmax = data.label$xend,
                                        tip_length = 0.01,
                                        textsize = 5)
final_boxplot_no_fac_ok

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final_boxplot_no_fac_ok <- final_boxplot_no_fac_ok +
                            scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) # Make room for the first label
                            
final_boxplot_no_fac_ok

## ----fig.dim=c(5,3.5), fig.align='center', cache=FALSE------------------------
library(tidyverse)
colnames(data.dCtF)
data.dCtF.slim <- pivot_longer(data.dCtF, cols = ANGPT1:VEGFC, names_to = "gene", values_to = "exp")

# Select genes
data.dCtF.slim_sel <- data.dCtF.slim[data.dCtF.slim$gene %in% c("ANGPT1","IL8","VEGFB"), ]

# Change order of groups if needed
data.dCtF.slim_sel$Group <- factor(data.dCtF.slim_sel$Group, levels = c("Control","AAA"))

# Create plot
final_boxplot_no_colors <- ggplot(data.dCtF.slim_sel, aes(x = Group, y = exp)) + 
                            geom_boxplot(outlier.shape = NA, coef = 2) +
                            theme_bw() +
                            ylab("dCt") +
                            xlab("") +
                            theme(axis.text = element_text(size = 10, color = "black")) + 
                            theme(axis.title = element_text(size = 10, color = "black")) +
                            theme(panel.grid.major.x = element_blank()) +
                            facet_wrap(vars(gene), nrow = 1, dir = "h", scales = "free")

final_boxplot_no_colors

## ----fig.dim=c(5,3.5), fig.align='center', cache=FALSE------------------------
data.label <- data.frame(matrix(nrow = 3, ncol = 4)) # Number of rows is equal to number of genes 
rownames(data.label) <- c("ANGPT1","IL8","VEGFB")
colnames(data.label) <- c("x", "xend", "y", "annotation")
data.label$gene <- rownames(data.label) # Name of column with gene symbols in this table must be 
# the same as name of the column with gene symbols in data used for create the plot.

data.label$y <- 0.5 + c(max(data.dCtF$ANGPT1), max(data.dCtF$IL8), max(data.dCtF$VEGFB))
data.label$x <- c(1,1,1)
data.label$xend <- c(1.98,1.98,1.98)
data.label$annotation <- c("****","**","***")

final_boxplot_no_colors_labels <- final_boxplot_no_colors +
 geom_signif(
    stat = "identity",
    data = data.label,
    aes(x = x,
        xend = xend,
        y = y,
        yend = y,
        annotation = annotation),
    color = "black",
    manual = TRUE) +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

final_boxplot_no_colors_labels

## ----fig.dim=c(5,3.5), fig.align='center', cache=FALSE------------------------
final_boxplot_no_colors_labels_points <- final_boxplot_no_colors_labels +
                                          geom_point(position=position_jitter(w=0.1,h=0), 
                                                     alpha = 0.7, 
                                                     size = 1.5)
                            
final_boxplot_no_colors_labels_points

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final_barplot <- results_barplot(data = data.dCtF,
                                 sel.Gene = c("ANGPT1","IL8", "VEGFB"),
                                 signif.show = TRUE,
                                 signif.labels = c("****","**","***"),
                                 angle = 30,
                                 signif.dist = 1.05,
                                 faceting = TRUE,
                                 facet.row = 1,
                                 facet.col = 4,
                                 y.exp.up = 0.1,
                                 y.axis.title = "dCt")

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final_barplot <- results_barplot(data = data.dCtF,
                                 sel.Gene = c("ANGPT1","IL8", "VEGFB"),
                                 signif.show = TRUE,
                                 signif.labels = c("****","**","***"),
                                 angle = 0,
                                 signif.dist = 1.05,
                                 faceting = FALSE,
                                 y.exp.up = 0.1,
                                 y.axis.title = "dCt")

# Add italic font to the x axis:
final_barplot <- final_barplot + 
                  theme(axis.text.x = element_text(face="italic")) 

final_barplot

## ----fig.dim=c(5,3.5), fig.align='center', cache=FALSE------------------------
library(tidyverse)
colnames(data.dCtF)
data.dCtF.slim <- pivot_longer(data.dCtF, cols = ANGPT1:VEGFC, names_to = "gene", values_to = "exp")

# Select genes
data.dCtF.slim_sel <- data.dCtF.slim[data.dCtF.slim$gene %in% c("ANGPT1","IL8","VEGFB"), ]

# Change order of groups if needed
data.dCtF.slim_sel$Group <- factor(data.dCtF.slim_sel$Group, levels = c("Control","AAA"))

data.mean <- data.dCtF.slim_sel %>%
              group_by(Group, gene) %>%
              summarise(mean = mean(exp, na.rm = TRUE), .groups = "keep")

data.sd <- data.dCtF.slim_sel %>%
            group_by(Group, gene) %>%
            summarise(sd = sd(exp, na.rm = TRUE), .groups = "keep")

data.mean$sd <- data.sd$sd

final_barplot_no_colors <- ggplot(data.mean, aes(x = Group, y = mean)) +
                              geom_errorbar(aes(group = Group,
                                                y = mean,
                                                ymin = ifelse(mean < 0, mean - abs(sd), mean),
                                                ymax = ifelse(mean > 0, mean + abs(sd), mean)),
                                            width = .2,
                                            position = position_dodge(0.9)) +
                
                              geom_col(aes(group = Group),
                                       position = position_dodge(0.9),
                                       width = 0.7,
                                       color = "black") +
                  xlab("") +
                  ylab("dCt") +
                  theme_bw() +
                  #theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
                  #theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
                  #theme(legend.text = element_text(size = legend.text.size, colour ="black")) +
                  #theme(legend.title = element_text(size = legend.title.size, colour =  "black")) +
                  theme(panel.grid.major.x = element_blank()) +
  facet_wrap(vars(gene), scales = "free", nrow = 1)

final_barplot_no_colors

## ----fig.dim=c(5,3.5), fig.align='center', cache=FALSE------------------------
data.label <- data.frame(matrix(nrow = 3, ncol = 4)) # Number of rows is equal to number of genes 
rownames(data.label) <- c("ANGPT1","IL8","VEGFB")
colnames(data.label) <- c("x", "xend", "y", "annotation")
data.label$gene <- rownames(data.label) # Name of column with gene symbols in this table 
# must be the same as name of the column with gene symbols in data used for create the plot.

data.mean <- data.mean %>%
               mutate(max = mean + sd) %>%
               group_by(gene) %>%
               summarise(height = max(max, na.rm = TRUE), .groups = "keep")
data.label$y <- 0.5 + data.mean$height
data.label$x <- c(1,1,1)
data.label$xend <- c(1.98,1.98,1.98)
data.label$annotation <- c("****","**","***")

final_barplot_no_colors_labels <- final_barplot_no_colors +
 geom_signif(
    stat = "identity",
    data = data.label,
    aes(x = x,
        xend = xend,
        y = y,
        yend = y,
        annotation = annotation,
        textsize = 5),
    color = "black",
    manual = TRUE) +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

final_barplot_no_colors_labels

## ----fig.dim=c(7.1,5), cache=FALSE--------------------------------------------
# Create named list with colors for groups annotation:
colors.for.groups = list("Group" = c("AAA"="firebrick1","Control"="green3"))
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
                cellwidth = 4) # It avoids cropping the image on the right side.

## ----fig.dim=c(5,5.5), fig.align='center', cache=FALSE------------------------
pca.kmeans <- pca_kmeans(data.dCt, 
                           sel.Gene = c("ANGPT1","IL8", "VEGFB"), 
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
                            p.adjust.method = "BH",
                            add.coef = "white")

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
AAA6_Control17 <- single_pair_sample(data = data.dCt,
                                     x = "AAA6",
                                     y = "Control17",
                                     point.size = 3,
                                     labels = TRUE,
                                     label = c("eq", "R2", "p"),
                                     label.position.x = 0.05)

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
library(ggpmisc)
PDGFB_TGFB <- single_pair_gene(data.dCt,
                               x = "PDGFB",
                               y = "TGFB",
                               by.group = TRUE,
                               point.size = 3,
                               labels = TRUE,
                               label = c("eq", "R2", "p"),
                               label.position.x = c(0.05),
                               label.position.y = c(1,0.95))

## ----cache=FALSE--------------------------------------------------------------
library(pROC)
# Remember to specify the numbers of rows (panels.row parameter) and columns (panels.col parameter) 
# to be sufficient to arrange panels: 
roc_parameters <- ROCh(data = data.dCt,
                       sel.Gene = c("ANGPT1","IL8", "VEGFB"),
                       groups = c("AAA","Control"),
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
                           sel.Gene = c("ANGPT1","IL8", "VEGFB"),
                           group.study = "AAA",
                           group.ref = "Control")
log.reg.results[[2]]

## ----fig.dim=c(5,4), fig.align='center', cache=FALSE--------------------------
log.reg.results.sorted <- log.reg.results[[1]] +
                           scale_y_discrete(limits = rev(sort(log.reg.results[[2]]$Gene)))
log.reg.results.sorted

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
data.dCt.pairwise <- delta_Ct(data = data.Ct.pairwiseF.ready, 
                              ref = "Gene4", 
                              transform = FALSE)

## ----cache=FALSE--------------------------------------------------------------
data.dCt.exp.pairwise <- delta_Ct(data = data.Ct.pairwiseF.ready,
                                  ref = "Gene4",
                                  transform = TRUE)

## ----cache=FALSE--------------------------------------------------------------
data.dCt.pairwise <- delta_Ct(data = data.Ct.pairwiseF.ready, 
                              ref = "Gene4", 
                              transform = FALSE)
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
data.dCt.pairwise <- delta_Ct(data = data.Ct.pairwiseF.ready, 
                              ref = "Gene4", 
                              transform = FALSE)
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

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
control.pca.gene.pairwise <- control_pca_gene(data = data.dCt.pairwise,
                                              hjust = 0.5)

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
control.pca.sample.pairwise <- control_pca_sample(data = FCh,
                                                  pairwise.FCh = TRUE,
                                                  colors = "black",
                                                  point.size = 3,
                                                  label.size = 2.5,
                                                  hjust = 0.5)

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
control.pca.gene.pairwise <- control_pca_gene(data = FCh,
                                              pairwise.FCh = TRUE,
                                              color = "black",
                                              hjust = 0.5)

## ----cache=FALSE--------------------------------------------------------------
data.dCt.pairwise.F <- filter_transformed_data(data = data.dCt.pairwise,
                                               remove.Sample = c("Sample22", 
                                                                 "Sample23", 
                                                                 "Sample15",
                                                                 "Sample03"))
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
final.boxplot.pairwise <- results_boxplot(data = data.dCt.pairwise.F,
                                          sel.Gene = c("Gene8","Gene19"),
                                          by.group = TRUE,
                                          signif.show = TRUE,
                                          signif.labels = c("***","*"),
                                          signif.dist = 1.2,
                                          faceting = FALSE,
                                          y.exp.up = 0.1,
                                          y.axis.title = "dCt")
# Add x axis annotations and ticks:
final.boxplot.pairwise <- final.boxplot.pairwise + 
                            theme(axis.text.x = element_text(size = 11, colour = "black", face="italic"),
                                  axis.ticks.x = element_line(colour = "black"))

final.boxplot.pairwise

## ----fig.dim=c(4,3.5), fig.align='center', cache=FALSE------------------------
library(tidyverse)
colnames(data.dCt.pairwise.F)
data.dCt.pairwise.F.slim <- pivot_longer(data.dCt.pairwise.F, 
                                         cols = Gene10:Gene8, 
                                         names_to = "gene", 
                                         values_to = "exp")

# Select genes
data.dCt.pairwise.F.slim.sel <- data.dCt.pairwise.F.slim[data.dCt.pairwise.F.slim$gene %in% c("Gene19","Gene8"), ]

# Change order of groups if needed
data.dCt.pairwise.F.slim.sel$Group <- factor(data.dCt.pairwise.F.slim.sel$Group,
                                             levels = c("Before","After"))

# Create plot
final_boxplot_no_colors <- ggplot(data.dCt.pairwise.F.slim.sel, aes(x = Group, y = exp)) + 
                            geom_boxplot(outlier.shape = NA, coef = 2) +
                            theme_bw() +
                            ylab("dCt") +
                            xlab("") +
                            theme(axis.text = element_text(size = 8, color = "black")) + 
                            theme(axis.title = element_text(size = 10, color = "black")) +
                            theme(panel.grid.major.x = element_blank()) +
                            facet_wrap(vars(gene), nrow = 1, dir = "h", scales = "free")

final_boxplot_no_colors

## ----fig.dim=c(5,3.5), fig.align='center', cache=FALSE------------------------
data.label <- data.frame(matrix(nrow = 2, ncol = 4)) # Number of rows is equal to number of genes 
rownames(data.label) <- c("Gene19","Gene8")
colnames(data.label) <- c("x", "xend", "y", "annotation")
data.label$gene <- rownames(data.label) # Name of column with gene symbols in this table 
# must be the same as name of the column with gene symbols in data used for create the plot.

data.label$y <- 0.5 + c(max(data.dCt.pairwise.F$Gene19), max(data.dCt.pairwise.F$Gene8))
data.label$x <- c(1,1)
data.label$xend <- c(1.98,1.98)
data.label$annotation <- c("***","*")

final_boxplot_no_colors_labels <- final_boxplot_no_colors +
 geom_signif(
    stat = "identity",
    data = data.label,
    aes(x = x,
        xend = xend,
        y = y,
        yend = y,
        annotation = annotation),
    color = "black",
    manual = TRUE) +
   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

final_boxplot_no_colors_labels

## ----fig.dim=c(5,3.5), fig.align='center', cache=FALSE------------------------
final_boxplot_no_colors_labels_points <- final_boxplot_no_colors_labels +
                                          geom_point(position=position_jitter(w=0.1,h=0), alpha = 0.7, size = 1.5)
                            
final_boxplot_no_colors_labels_points

## ----fig.dim=c(4,4.5), fig.align='center', cache=FALSE------------------------
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

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final.barplot.pairwise <- results_barplot(data = data.dCt.pairwise.F,
                                          sel.Gene = c("Gene8","Gene19"),
                                          signif.show = TRUE,
                                          signif.labels = c("***","*"),
                                          angle = 0,
                                          signif.dist = 1.05,
                                          faceting = FALSE,
                                          y.exp.up = 0.1,
                                          y.axis.title = "dCt")

# Add italic font to the x axis:
final.barplot.pairwise <- final.barplot.pairwise + 
                            theme(axis.text.x = element_text(face="italic")) 

final.barplot.pairwise

## ----fig.dim=c(5,3.5), fig.align='center', cache=FALSE------------------------
library(tidyverse)
colnames(data.dCt.pairwise.F)
data.dCt.pairwise.F.slim <- pivot_longer(data.dCt.pairwise.F, 
                                         cols = Gene10:Gene8, 
                                         names_to = "gene", 
                                         values_to = "exp")

# Select genes
data.dCt.pairwise.F.slim.sel <- data.dCt.pairwise.F.slim[data.dCt.pairwise.F.slim$gene %in% c("Gene19","Gene8"), ]

# Change order of groups if needed
data.dCt.pairwise.F.slim.sel$Group <- factor(data.dCt.pairwise.F.slim.sel$Group, 
                                             levels = c("Before","After"))

data.mean <- data.dCt.pairwise.F.slim.sel %>%
              group_by(Group, gene) %>%
              summarise(mean = mean(exp, na.rm = TRUE), .groups = "keep")

data.sd <- data.dCt.pairwise.F.slim.sel %>%
            group_by(Group, gene) %>%
            summarise(sd = sd(exp, na.rm = TRUE), .groups = "keep")

data.mean$sd <- data.sd$sd

final_barplot_no_colors <- ggplot(data.mean, aes(x = Group, y = mean)) +
                              geom_errorbar(aes(group = Group,
                                                y = mean,
                                                ymin = ifelse(mean < 0, mean - abs(sd), mean),
                                                ymax = ifelse(mean > 0, mean + abs(sd), mean)),
                                            width = .2,
                                            position = position_dodge(0.9)) +
                
                              geom_col(aes(group = Group),
                                       position = position_dodge(0.9),
                                       width = 0.7,
                                       color = "black") +
                              xlab("") +
                              ylab("dCt") +
                              theme_bw() +
                              theme(panel.grid.major.x = element_blank()) +
                              facet_wrap(vars(gene), scales = "free", nrow = 1)

final_barplot_no_colors

## ----fig.dim=c(5,3.5), fig.align='center', cache=FALSE------------------------
data.label <- data.frame(matrix(nrow = 2, ncol = 4)) # Number of rows is equal to number of genes 
rownames(data.label) <- c("Gene19","Gene8")
colnames(data.label) <- c("x", "xend", "y", "annotation")
data.label$gene <- rownames(data.label) # Name of column with gene symbols 
# in this table must be the same as name of the column with gene symbols 
# in data used for create the plot.

data.mean <- data.mean %>%
               mutate(max = mean + sd) %>%
               group_by(gene) %>%
               summarise(height = max(max, na.rm = TRUE), .groups = "keep")
data.label$y <- 0.5 + data.mean$height
data.label$x <- c(1,1)
data.label$xend <- c(1.98,1.98)
data.label$annotation <- c("***","*")

final_barplot_no_colors_labels <- final_barplot_no_colors +
                                   geom_signif(
                                   stat = "identity",
                                   data = data.label,
                                   aes(x = x,
                                      xend = xend,
                                      y = y,
                                      yend = y,
                                      annotation = annotation,
                                      textsize = 5),
                                   color = "black",
                                   manual = TRUE) +
                                   scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

final_barplot_no_colors_labels

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
# Cellwidth parameter was set to 8 to avoid cropping the image on the right side.

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
# Remember to specify the numbers of rows (panels.row parameter) and columns (panels.col parameter) 
# to provide sufficient place to arrange panels: 
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
                           log.axis = TRUE,
                           p.adjust = FALSE)
log.reg.results[[2]]

## ----fig.dim=c(5,4), fig.align='center', cache=FALSE--------------------------
log.reg.results.sorted <- log.reg.results[[1]] +
                           scale_y_discrete(limits = rev(sort(log.reg.results[[2]]$Gene)))
log.reg.results.sorted

## ----cache=FALSE--------------------------------------------------------------
data("data.Ct.3groups")
str(data.Ct.3groups)
table(data.Ct.3groups$Group)

## ----fig.dim=c(7.1,9), cache=FALSE--------------------------------------------
sample.Ct.control.3groups <- control_Ct_barplot_sample(data = data.Ct.3groups,
                                                       flag.Ct = "Undetermined",
                                                       maxCt = 35,
                                                       flag = c("Undetermined"),
                                                       axis.title.size = 9,
                                                       axis.text.size = 7,
                                                       plot.title.size = 9,
                                                       legend.title.size = 9,
                                                       legend.text.size = 9)


## ----fig.dim=c(7.1,5.5), cache=FALSE------------------------------------------
gene.Ct.control.3groups <- control_Ct_barplot_gene(data = data.Ct.3groups,
                                                   flag.Ct = "Undetermined",
                                                   maxCt = 35,
                                                   flag = c("Undetermined"),
                                                   axis.title.size = 9,
                                                   axis.text.size = 9,
                                                   plot.title.size = 9,
                                                   legend.title.size = 9,
                                                   legend.text.size = 9)

## ----cache=FALSE--------------------------------------------------------------
head(sample.Ct.control.3groups[[2]])

## ----cache=FALSE--------------------------------------------------------------
head(gene.Ct.control.3groups[[2]])

## ----fig.dim=c(7.1,9), cache=FALSE--------------------------------------------
library(tidyverse)
library(pheatmap)
data("data.Ct.3groups")

# Vector of colors to fill the heatmap can be specified to fit the user needs:
colors <- c("#4575B4","#FFFFBF","#C32B23")
control_heatmap(data.Ct.3groups,
                sel.Gene = "all",
                colors = colors,
                show.colnames = TRUE,
                show.rownames = TRUE, 
                fontsize = 9,
                fontsize.row = 9,
                angle.col = 45)

## ----cache=FALSE--------------------------------------------------------------
# Finding samples with more than half of the unreliable Ct values.
low.quality.samples.3groups <- filter(sample.Ct.control.3groups[[2]], Not.reliable.fraction > 0.5)$Sample
low.quality.samples.3groups <- as.vector(low.quality.samples.3groups)                        
low.quality.samples.3groups

## -----------------------------------------------------------------------------
# Finding genes with more than half of the unreliable Ct values in given group.
low.quality.genes.3groups <- filter(gene.Ct.control.3groups[[2]], Not.reliable.fraction > 0.5)$Gene
low.quality.genes.3groups <- unique(as.vector(low.quality.genes.3groups))                        
low.quality.genes.3groups

## ----cache=FALSE--------------------------------------------------------------
# Data filtering
data.CtF.3groups <- filter_Ct(data = data.Ct.3groups,
                              flag.Ct = "Undetermined",
                              maxCt = 35,
                              flag = c("Undetermined"),
                              remove.Gene = low.quality.genes.3groups,
                              remove.Sample = low.quality.samples.3groups)

# Collapsing technical replicates without imputation:
data.CtF.ready.3groups <- make_Ct_ready(data = data.CtF.3groups,
                                        imput.by.mean.within.groups = FALSE)
# A part of the data with missing values:
as.data.frame(data.CtF.ready.3groups)[25:30,]

# Collapsing technical replicates with imputation:
data.CtF.ready.3groups <- make_Ct_ready(data = data.CtF.3groups,
                                        imput.by.mean.within.groups = TRUE)
# Missing values were imputed:
as.data.frame(data.CtF.ready.3groups)[25:30,]



## ----fig.dim=c(7.1,4), cache=FALSE--------------------------------------------
library(ctrlGene)
# Remember that the number of colors in col parameter should be equal to the number of tested genes:
ref.3groups <- find_ref_gene(data = data.CtF.ready.3groups,
                             groups = c("AAA","Control","VV"),
                             candidates = c("CCL5", "GAPDH","IL1B","TGFB", "VEGFA"),
                             col = c("#66c2a5", "#fc8d62","#6A6599", "#1F77B4", "black"),
                             angle = 60,
                             axis.text.size = 7,
                             norm.finder.score = TRUE,
                             genorm.score = TRUE)
ref.3groups[[2]]

## ----cache=FALSE--------------------------------------------------------------
# For 2-dCt method:
data.dCt.exp.3groups <- delta_Ct(data = data.CtF.ready.3groups,
                                 normalise = TRUE,
                                 ref = "VEGFA",
                                 transform = TRUE)

# For 2-ddCt method:
data.dCt.3groups <- delta_Ct(data = data.CtF.ready.3groups,
                             normalise = TRUE,
                             ref = "VEGFA",
                             transform = FALSE)

## ----fig.dim=c(7.1,8), cache=FALSE--------------------------------------------
control_boxplot_sample_3groups <- control_boxplot_sample(data = data.dCt.3groups,
                                                         colors = c("#66c2a5", "#fc8d62", "#8DA0CB"),
                                                         axis.text.size = 7)

## ----fig.dim=c(7.1,4)---------------------------------------------------------
control_boxplot_gene_3groups <- control_boxplot_gene(data = data.dCt.3groups,
                                                     by.group = TRUE,
                                                     colors = c("#66c2a5", "#fc8d62", "#8DA0CB"),
                                                     axis.text.size = 10)

## ----fig.dim=c(5,5.5), fig.align='center', cache=FALSE------------------------
control.pca.sample.3groups <- control_pca_sample(data = data.dCt.3groups,
                                                 colors = c("#66c2a5", "#fc8d62", "#8DA0CB"),
                                                 point.size = 3,
                                                 label.size = 2.5,
                                                 legend.position = "top")

## ----fig.dim=c(5,5), fig.align='center', cache=FALSE--------------------------
control.pca.gene.3groups <- control_pca_gene(data = data.dCt.3groups)

## ----fig.dim=c(7.1,4), cache=FALSE--------------------------------------------
control_cluster_sample(data = data.dCt.3groups,
                       method.dist = "euclidean",
                       method.clust = "average",
                       label.size = 0.5)
control_cluster_gene(data = data.dCt.3groups,
                     method.dist = "euclidean",
                     method.clust = "average",
                     label.size = 0.5)

## ----cache=FALSE--------------------------------------------------------------
data.dCtF.3groups <- filter_transformed_data(data = data.dCt.3groups,
                                             remove.Sample = c("AAA14"))

## ----cache=FALSE--------------------------------------------------------------
data.dCt.exp.3groups <- delta_Ct(data = data.CtF.ready.3groups,
                                 ref = "VEGFA",
                                 transform = TRUE)
library(coin)
results.dCt.3groups <- RQ_dCt(data = data.dCt.exp.3groups,
                              do.tests = TRUE,
                              group.study = "VV",
                              group.ref = "Control")

# Obtained table can be sorted by, e.g. p values from the Mann-Whitney U test:
head(as.data.frame(arrange(results.dCt.3groups, MW_test_p)))

## ----cache=FALSE--------------------------------------------------------------
data.dCt.3groups <- delta_Ct(data = data.CtF.ready.3groups,
                             ref = "VEGFA",
                             transform = FALSE)
library(coin)
results.ddCt.3groups <- RQ_ddCt(data = data.dCt.3groups,
                                group.study = "VV",
                                group.ref = "Control",
                                do.tests = TRUE)

# Obtained table can be sorted by, e.g. p values from the Mann-Whitney U test:
head(as.data.frame(arrange(results.ddCt.3groups, MW_test_p)))

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final_boxplot_3groups <- results_boxplot(data = data.dCtF.3groups,
                                 sel.Gene = c("ANGPT1","VEGFB", "VEGFC"),
                                 colors = c("#66c2a5", "#fc8d62", "#8DA0CB"),
                                 by.group = TRUE,
                                 signif.show = TRUE,
                                 signif.labels = c("****","***","*"),
                                 signif.dist = 1.05,
                                 faceting = TRUE,
                                 facet.row = 1,
                                 facet.col = 4,
                                 y.exp.up = 0.1,
                                 angle = 20,
                                 y.axis.title = "dCt")

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final_barplot_3groups <- results_barplot(data = data.dCtF.3groups,
                                 sel.Gene = c("ANGPT1","VEGFB","VEGFC"),
                                 colors = c("#66c2a5", "#fc8d62", "#8DA0CB"),
                                 signif.show = TRUE,
                                 signif.labels = c("****","***","*"),
                                 angle = 30,
                                 signif.dist = 1.05,
                                 faceting = TRUE,
                                 facet.row = 1,
                                 facet.col = 4,
                                 y.exp.up = 0.1,
                                 y.axis.title = "dCt")

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
# Draw plot without statistical significance labels:
final_boxplot_3groups <- results_boxplot(data = data.dCtF.3groups,
                                         sel.Gene = c("ANGPT1","VEGFB", "VEGFC"),
                                         colors = c("#66c2a5", "#fc8d62", "#8DA0CB"),
                                         by.group = TRUE,
                                         signif.show = FALSE, # It avoids drawing labels
                                         faceting = TRUE,
                                         facet.row = 1,
                                         facet.col = 4,
                                         y.exp.up = 0.2,
                                         angle = 20,
                                         y.axis.title = "dCt")

## ----cache=FALSE--------------------------------------------------------------
# Prepare expression data:
data <- pivot_longer(data.dCtF.3groups,
                     !c(Sample, Group),
                     names_to = "Gene" ,
                     values_to = "value")

# filter for genes:
data <- filter(data, Gene %in% c("ANGPT1","VEGFB", "VEGFC"))

# Find maximum value in each group:
label.height <- data %>%
                 group_by(Gene) %>%
                 summarise(height = max(value), .groups = "keep")

# Prepare empty data frame:
data.label.empty <- data.frame(matrix(nrow = length(unique(label.height$Gene)), ncol = 4))
rownames(data.label.empty) <- label.height$Gene
colnames(data.label.empty) <- c("x", "xend", "y", "annotation")
data.label.empty$Gene <- rownames(data.label.empty)

# Fill a data frame with coordinates for right pair:
data.label.right <- data.label.empty
data.label.right$x <- rep(1.01, nrow(data.label.right))
data.label.right$xend <- rep(1.25, nrow(data.label.right))
data.label.right$y <- label.height$height + 0.5
data.label.right$annotation <- c("right1","right2","right3")

# Fill a data frame with coordinates for left pair:
data.label.left <- data.label.empty
data.label.left$x <- rep(0.98, nrow(data.label.left))
data.label.left$xend <- rep(0.75, nrow(data.label.left))
data.label.left$y <- label.height$height + 0.5
data.label.left$annotation <- c("left1","left2","left3")

# Fill a data frame with coordinates for edge pair:
data.label.edge <- data.label.empty
data.label.edge$x <- rep(0.75, nrow(data.label.edge))
data.label.edge$xend <- rep(1.25, nrow(data.label.edge))
data.label.edge$y <- label.height$height + 1.2
data.label.edge$annotation <- c("edge1","edge2","edge3")

## ----fig.dim=c(6.5,4.5), fig.align='center', cache=FALSE----------------------
final_boxplot_3groups +
  geom_signif(
          stat = "identity",
          data = data.label.right,
          aes(x = x,
              xend = xend,
              y = y,
              yend = y,
              annotation = annotation),
          color = "black",
          manual = TRUE) +
  
  geom_signif(
          stat = "identity",
          data = data.label.left,
          aes(x = x,
              xend = xend,
              y = y,
              yend = y,
              annotation = annotation),
          color = "black",
          manual = TRUE) +
  
   geom_signif(
          stat = "identity",
          data = data.label.edge,
          aes(x = x,
              xend = xend,
              y = y,
              yend = y,
              annotation = annotation),
          color = "black",
          manual = TRUE) +

  scale_y_continuous(expand = expansion(mult = c(0.1, 0.12))) # it makes space for labels

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
final_barplot_3groups <- results_barplot(data = data.dCtF.3groups,
                                         sel.Gene = c("ANGPT1","VEGFB","VEGFC"),
                                         colors = c("#66c2a5", "#fc8d62", "#8DA0CB"),
                                         signif.show = FALSE,
                                         angle = 30,
                                         faceting = TRUE,
                                         facet.row = 1,
                                         facet.col = 4,
                                         y.exp.up = 0.1,
                                         y.axis.title = "dCt")

## ----cache=FALSE--------------------------------------------------------------
# Prepare expression data:
data <- pivot_longer(data.dCtF.3groups,
                     !c(Sample, Group),
                     names_to = "Gene" ,
                     values_to = "value")

# filter for genes:
data <- filter(data, Gene %in% c("ANGPT1","VEGFB", "VEGFC"))

# Calculate mean and standard deviation for each group:
data.mean <- data %>%
             group_by(Group, Gene) %>%
             summarise(mean = mean(value, na.rm = TRUE), .groups = "keep")

data.sd <- data %>%
           group_by(Group, Gene) %>%
           summarise(sd = sd(value, na.rm = TRUE), .groups = "keep")

data.mean$sd <- data.sd$sd

#Find the highest values:
label.height <- data.mean %>%
                mutate(max = mean + sd) %>%
                group_by(Gene) %>%
                summarise(height = max(max, na.rm = TRUE), .groups = "keep")

# Prepare empty data frame:
data.label.empty <- data.frame(matrix(nrow = length(unique(data.mean$Gene)), ncol = 4))
rownames(data.label.empty) <- unique(data.mean$Gene)
colnames(data.label.empty) <- c("x", "xend", "y", "annotation")
data.label.empty$Gene <- rownames(data.label.empty)

# Fill a data frame with coordinates for left pair:
data.label.left <- data.label.empty
data.label.left$x <- rep(0.97, nrow(data.label.left))
data.label.left$xend <- rep(0.7, nrow(data.label.left))
data.label.left$y <- label.height$height + 0.3
data.label.left$annotation <- c("left1","left2","left3")

# Fill a data frame with coordinates for left pair:
data.label.right <- data.label.empty
data.label.right$x <- rep(1.01, nrow(data.label.right))
data.label.right$xend <- rep(1.28, nrow(data.label.right))
data.label.right$y <- label.height$height + 0.3
data.label.right$annotation <- c("right1","right2","right3")

# Fill a data frame with coordinates for edge pair:
data.label.edge <- data.label.empty
data.label.edge$x <- rep(0.7, nrow(data.label.edge))
data.label.edge$xend <- rep(1.28, nrow(data.label.edge))
data.label.edge$y <- label.height$height + 1
data.label.edge$annotation <- c("edge1","edge2","edge3")

## ----fig.dim=c(6,4.5), fig.align='center', cache=FALSE------------------------
final_barplot_3groups +
  geom_signif(
          stat = "identity",
          data = data.label.left,
          aes(x = x,
              xend = xend,
              y = y,
              yend = y,
              annotation = annotation),
          color = "black",
          manual = TRUE) +
  
   geom_signif(
          stat = "identity",
          data = data.label.right,
          aes(x = x,
              xend = xend,
              y = y,
              yend = y,
              annotation = annotation),
          color = "black",
          manual = TRUE) +
  
   geom_signif(
          stat = "identity",
          data = data.label.edge,
          aes(x = x,
              xend = xend,
              y = y,
              yend = y,
              annotation = annotation),
          color = "black",
          manual = TRUE) +
  
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.12)))

## ----fig.dim=c(9,5), cache=FALSE----------------------------------------------
# Create named list with colors for groups annotation:
colors.for.groups = list("Group" = c("AAA"="#f98517","Control"="#33b983", "VV"="#bf8cfc"))
# Vector of colors for heatmap:
colors <- c("navy","navy","#313695","#4575B4","#74ADD1","#ABD9E9",
            "#E0F3F8","#FFFFBF","#FEE090","#FDAE61","#F46D43",
            "#D73027","#C32B23","#A50026","#8B0000", 
            "#7E0202","#000000")
results_heatmap(data.dCtF.3groups,
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
pca.kmeans <- pca_kmeans(data.dCtF.3groups, 
                         sel.Gene = c("ANGPT1","VEGFB", "VEGFC"),
                         k.clust = 3,
                         clust.names = c("Cluster1", "Cluster2", "Cluster3"),
                         point.shape = c(19, 17, 18),
                         point.color = c("#66c2a5", "#fc8d62", "#8DA0CB"),
                         legend.position = "top")
# Access to the confusion matrix:
pca.kmeans[[2]]

## ----fig.dim=c(5,6), fig.align='center', cache=FALSE--------------------------
pca.kmeans[[1]] + theme(legend.box = "vertical")

## ----fig.dim=c(6.5,6.5), fig.align='center', cache=FALSE----------------------
library(Hmisc)
library(corrplot)
# To make the plot more readable, only part of the data was used:
corr.samples <- corr_sample(data = data.dCtF.3groups[15:30, ],
                            method = "pearson",
                            order = "hclust",
                            size = 0.7,
                            p.adjust.method = "BH",
                            add.coef = "white")

## ----fig.dim=c(6.5,6.5), fig.align='center', cache=FALSE----------------------
library(Hmisc)
library(corrplot)
corr.genes <- corr_gene(data = data.dCtF.3groups,
                        method = "spearman",
                        order = "FPC",
                        size = 0.7,
                        p.adjust.method = "BH")

## ----fig.dim=c(4.5,4.5), fig.align='center', cache=FALSE----------------------
library(ggpmisc)
AAA6_AAA43 <- single_pair_sample(data = data.dCtF.3groups,
                                         x = "AAA6",
                                         y = "AAA43",
                                         point.size = 3,
                                         labels = TRUE,
                                         label = c("eq", "R2", "p"),
                                         label.position.x = 0.05)

## ----fig.dim=c(5,4.5), fig.align='center', cache=FALSE------------------------
library(ggpmisc)
# Without stratification by groups:
PDGFB_TGFB <- single_pair_gene(data.dCtF.3groups,
                               x = "PDGFB",
                               y = "TGFB",
                               by.group = FALSE,
                               point.size = 3,
                               labels = TRUE,
                               label = c("eq", "R2", "p"),
                               label.position.x = c(0.05),
                               label.position.y = c(1,0.95))

## ----fig.dim=c(6,5.5), fig.align='center', cache=FALSE------------------------
library(ggpmisc)
# With stratification by groups:
PDGFB_TGFB <- single_pair_gene(data.dCtF.3groups,
                               x = "PDGFB",
                               y = "TGFB",
                               by.group = TRUE,
                               colors = c("#66c2a5", "#fc8d62", "#8DA0CB"), # Vector of colors
                               point.size = 3,
                               labels = TRUE,
                               label = c("eq", "R2", "p"),
                               label.position.x = c(0.05),
                               label.position.y = c(1,0.95,0.9)) # Labels position

## ----cache=FALSE--------------------------------------------------------------
library(pROC)
# Remember to specify the numbers of rows (panels.row parameter) and columns (panels.col parameter) to be sufficient to arrange panels: 
roc_parameters <- ROCh(data = data.dCtF.3groups,
                       sel.Gene = c("ANGPT1","VEGFB", "VEGFC"),
                       groups = c("Control","VV"),
                       panels.row = 2,
                       panels.col = 2)
roc_parameters

## ----echo=FALSE, out.width="500px", fig.align="center", warning=FALSE, message=FALSE, cache=FALSE----
knitr::include_graphics("ROC_plot_3groups.png")

## ----fig.dim=c(5,4), fig.align='center', cache=FALSE--------------------------
library(oddsratio)
# Remember to set the increment parameter.
log.reg.results <- log_reg(data = data.dCtF.3groups,
                           increment = 1,
                           sel.Gene = c("ANGPT1","VEGFB", "VEGFC"),
                           group.study = "VV",
                           group.ref = "Control")
log.reg.results[[2]]

## ----fig.dim=c(5,4), fig.align='center', cache=FALSE--------------------------
log.reg.results.sorted <- log.reg.results[[1]] +
                           scale_y_discrete(limits = rev(sort(log.reg.results[[2]]$Gene)))
log.reg.results.sorted

## -----------------------------------------------------------------------------
sessionInfo()

