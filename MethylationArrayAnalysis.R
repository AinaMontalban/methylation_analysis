###################################################
#           Methylation array analysis            #
###################################################
#
# 1. Data loading
#     1.1 Targets
#     1.2 IDAT files
# 2. Quality control
#     2.1 Heatmap SNP
#     2.2 Failed samples with det-P -> remove them
#     2.3 Proportion of failed Probes 
#     2.4 Density plots
#     2.5 BeanPlots
#     2.6 Control Type
#     2.7 Sex Prediction
#     2.8 QC report
# 3. Normalization
#     3.1 Normalization method: ssNoob 
# 4. Filtering
#     4.1 Remove failed probes (CpG; or mark as NA)
#     4.2 Remove probes on the sex chromosomes
#     4.3 Remove probes with SNPs at CpG site and at a single nucleotide extension
#     4.4 Exclude cross reactive probes
#     4.5 Optional (?) Filtering by variance
# 5. Data Exploration
#     5.1 MDS
#     5.2 t-SNE
#     5.3 UMAP
#     5.4 Examine higher dimensions
# 6. Calculte beta-values and m-values filtered
# 7. Probe-wise differential methylation analysis
#     7.1 Paired-test check
#     7.2 Table with the differential methylated CpGs (DMPs)
#     7.3 Regions differentially methylated (DMRs)
#     7.4 Heatmap with the significant analysis
# 8. Gene Ontology Enrichment (cluster profiler)
#
###################################################

###################################################
## Packages needed
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(Rtsne)
library(missMethyl)
library(matrixStats)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(gplots)
library(ggplot2)
library(svglite)

#=====================
#### Master module

# Full path to the directory
#path <- "/home/aina/Internship/methylation-minfi"
path <- "/mnt/ElRaid/amontalban/PROJECTS/methylation"

# Path to the samplesheet
#samplesheet_path <- file.path(path, "targets_subset")
samplesheet_path <- file.path(path, "data")

# Path to the IDAT files
#IDAT_path <- file.path(path, "data_subset")
IDAT_path <- file.path(path, "tmp")

# Path to the folder to save the results
#results_path <- file.path(path, "results")
results_path <- "/home/amontalban/Documents/methylation-minfi/results"
setwd(results_path)
#=====================
# 1. Data loading
## 1.1 Targets
targets <- read.metharray.sheet(samplesheet_path, pattern="sample_sheet_PDX.csv")
head(targets)

## Notice that the Basename column is not correct. It musts specify the location
## of each individual IDAT file in the experiment. (w/ the path)
targets$Basename <- file.path(IDAT_path, targets$Sample_Plate, targets$Slide, 
                              paste(targets$Slide, targets$Array, sep = "_"))
# Example: "/home/aina/Internship/methylation-minfi/data_subset/SMET0224/201496850086/201496850086_R02C01"
# In my case, we need: "/home/aina/Internship/methylation-minfi/data_subset/201496850086_R02C01"

targets$Basename <- file.path(IDAT_path, 
                              paste(targets$Slide, targets$Array, sep = "_"))

# Enforce Sample_Name as character.
targets$Sample_Name <- as.character(targets$Sample_Name)
print(paste("Full sample sheet rows (samples):", nrow(targets)))

### Let's study the sample sheet column by column 

colnames(targets)

## The data come from PDXs:
# Tumor tissue that has been taken from a patient and implanted into mice for research purposes. 
# Cancer drugs and other types of treatment may be tested on xenografts to see how well 
#they work before they are given to the patient. 
# Patient-derived xenografts may be used to help plan treatment and learn what the best treatment may be for a patient. 
# They are also being used in the development of new cancer drugs. Also called PDX. 
# NCI Dictionary of Cancer Terms. https://www.cancer.gov/publications/dictionaries/cancer-terms/def/patient-derived-xenograft 
####========================================================
## Sample_Name: it contains the name of the sample. (i.e  B01_Plate3)
## Sample Well: the well containing the specific sample in the plate (i.e B01)
## Sample_Plate: name for the plate containing bisulfite-converted DNA samples (i.e SMET0224)
## Sample_Group: name of the sample group (i.e NonResponder)
## Pool_ID: (i.e EPIC)
## Investigator: (i.e Livio)
## Project: in our data we have two projects (i.e [1] "TS01030062" "TS01030102") unique(targets$Project)
## Site.Primary: ???? (ie not specified)
## Site.Subtype: ???? (ie not specified)
## Type: ???? (ie not specified)
## Cell_Line: No (paired-test?) no
## Wildtype: no
## Normal: no
## Scan_Date: 04-Nov-2017
## Sample_Type: ????????
## Barcode
## Array: position on the sample on the Beadchip
## Slide
## Basename
## Sentrix_ID: BeadChip ID
####========================================================

### cgenalogy: CRC0291LMX0A02204
## The first 7 characters indicate the material obtained from a surgery. CRC0291
## If the patient has submitted to different surgeries the number will be different.

## LMX indicate liver metastasis

## --H: human, --X: xenograft

# Remove samples without group information
targets <- targets[targets$Sample_Group != "", ]
print(paste("Sample sheet rows (samples) after removing samples w/o info:", 
            nrow(targets)))

## 1.2 IDAT files
# Read in the raw data from the IDAT files.
rgSet <- read.metharray.exp(targets=targets, force = TRUE)
#annot <- paste(rgSet@annotation[[1]], rgSet@annotation[[2]], sep=".")
annot <- getAnnotation(rgSet)

#=====================
# 2. Quality control
# 2.1 Heatmap SNP

## In order to check if we have samples from the same individual, we can
## study the beta-values of the SNPs. Those samples that are equal, in other words, 
## they are from the same patient, will have the same SNPs.
png(file = file.path(results_path, "SNP_heatmap.png"),
    width = 4000,
    height = 4000,
    res = 250)
snps <- getSnpBeta(rgSet)
colnames(snps) <- targets$Sample_Name
my_palette <- colorRampPalette(c("yellow", "orange", "darkred"))(n = 299)
heatmap.2(snps, col = my_palette, trace = "none")
dev.off()

## d3heatmap
require(d3heatmap)
d3heatmap(snps, col = my_palette, trace = "none")
# There are same individuals


corrmatrix <- cor(snps)
require(reshape2)
require(plotly)
long <- melt(corrmatrix)

pg <- ggplot(long, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())

p<-ggplot(as.data.frame(long), aes(x=Var1, y=Var2, fill=value)) + geom_tile() + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
  scale_fill_gradient2(low="white", high ="darkred", mid ="lightgrey" , midpoint = 0.5)
ggplotly(p) 

corrmatrix[upper.tri(corrmatrix)] <- 0
subset_high_correlation <- melt(corrmatrix)
paired <- subset(subset_high_correlation, value>0.90 & value != 1.0)
nrow(paired)
paired1 <- subset(long, value>0.90 & value != 1.0)

require(corrplot)
?cor
corrplot(corrmatrix,  tl.pos='n', insig = "label_sig", pch.col = "white")
res1 <- cor.mtest(snps, conf.level = .95)
corrplot(corrmatrix,  diag = FALSE, order = "FPC",
         tl.pos = "n", tl.cex = 0.5, method = "color", type = "upper")

write.csv(corrmatrix, "correlation_matrix_SNPs.csv")
predictedSex <- getSex(GRSetSq, cutoff = -2)$predictedSex
names(colData(GRSetSq))
addSex(GRSetSq)
plotSex(addSex(GRSetSq))
same_patients_idx <- unique(c((paired$Var1), (paired$Var2)))
patGRSet <- GRSetSq[,same_patients_idx]
same_patients <- (c(as.character(paired$Var1), as.character(paired$Var2)))

same_patients_idx_nu <- (c((paired$Var1), (paired$Var2)))


create_idx_pairs(paired)







# 2.2 Failed samples with det-P -> remove them
detP <- detectionP(rgSet)
head(detP)
pval_means <- colMeans(detP)
names(pval_means) <- targets$Sample_Name
df_pval_means <- as.data.frame(pval_means)
colnames(df_pval_means) <- "pvals"
p <- ggplot(df_pval_means, aes(x=rownames(df_pval_means), y=pvals)) +
  geom_col(show.legend = FALSE, color="darkgrey") +
  theme_classic() + 
  scale_y_continuous(limits=c(0,0.08)) + 
  geom_hline(yintercept = 0.05, color="red") + 
  theme(axis.text.x = element_text(hjust = 1, angle=90, size = 3)) +
  labs(y="P-values", x="")+
  geom_hline(yintercept = 0.01, color="green")
#ggsave(filename = "detection_pvals.svg", plot = p, width = 6, height = 4)
ggsave(filename = "detection_pvals.png", plot = p, width = 6, height = 4)

## In the plot, we can observe any failed samples


# 2.3 Proportion of failed Probes 
## David Piñeyro code
probe_detP_max = 0.01
# Also examine the proportion of failed probes (detP >= probe_detP_max).
failed_proportion <- apply(detP, 2, function(x) {
  # All probes detP == 0 means totally failed array.
  if (sum(x) == 0) {
    return(NA)
  } else {
    failed_probes <- sum(x >= probe_detP_max)
    return(failed_probes / length(x))
  }
})
names(failed_proportion) <- targets$Sample_Name
png(filename = "failedPropProbes.png",  height = 1000,
    width = 2000,
    res = 200)
barplot(failed_proportion,  
        ylim = c(0, 0.12),
        las=2, 
        cex.names=0.8, ylab="Proportion of failed probes")
abline(h=0.1,col="red")
abline(h=0.05,col="blue")
abline(h= 0.01, col = "green")
legend("topleft", legend=levels(factor(targets$Sample_Group)), 
       bg="white", cex = 0.4)
dev.off()
# Output the sample sheet, but with an added column with the failed cpg
# proportion.
targets$Prop_failed_probes <- failed_proportion
write.csv(targets, file = file.path(results_path, 
                                    "Sample_sheet_with_prop_failed_cpgs.csv"))


#  2.4 Density plots
# B-values
bVals_raw <- getBeta(rgSet)
png(filename = "raw_bvals_densityPlot.png")
densityPlot(bVals_raw, sampGroups = targets$Sample_Group, legend=FALSE, main="Raw")
legend("top", legend = levels(factor(targets$Sample_Group)), text.col = brewer.pal(8,"Dark2"))
dev.off()

mset <- preprocessRaw(rgSet)
mVals_raw <- getM(mset)
png(filename = "raw_mvals_densityPlot.png")
densityPlot(bVals_raw, sampGroups = targets$Sample_Group, legend=FALSE, xlab = "M Values", main = "Raw")
legend("top", legend = levels(factor(targets$Sample_Group)), text.col = brewer.pal(8,"Dark2"))
dev.off()


#  2.5 BeanPlots

numPosition = 10000
idx <- sample(nrow(bVals_raw), numPosition)
b_subset <- as.matrix(bVals_raw[idx, ])
x <- melt(b_subset, varnames = c("cpg", "sample"))
o <- order(colnames(b))
sb <- x[1:30, ]
ggplot(data=sb, aes(x=sample, y=value)) + geom_violin() + coord_flip()

png(filename = "bvals_densityBeanPlot.png",   height = 2000,
                                                width = 2000,
                                                res = 200)
densityBeanPlot((bVals_raw), sampGroups = targets$Sample_Group, sampNames = targets$Sample_Name)
axis(1,cex.axis=1)
dev.off()


head(bVals_raw)

#     2.6 Control Type
controlType <- c("BISULFITE CONVERSION I",
                 "BISULFITE CONVERSION II",
                 "EXTENSION",
                 "HYBRIDIZATION",
                 "NEGATIVE",
                 "NON-POLYMORPHIC",
                 "SPECIFICITY I",
                 "SPECIFICITY II",
                 "TARGET REMOVAL",
                 "STAINING")

redControls <- getRed(rgSet)
greenControls <- getGreen(rgSet)

log2_subset_GC <- log2(subset)
df_subset_GC <- melt(log2_subset_GC)
ggplot(data=as.data.frame(df_subset_GC), aes(x=Var2, y=value)) + 
  geom_point(color="red", size=1.5) + scale_y_continuous(limits = c(-1, 20)) + 
  theme(axis.text.x = element_text(hjust = 1, angle=45)) +
  geom_hline(yintercept =threshold, linetype="dashed") + ylab("Log2 Intensity") + 
  scale_x_discrete(labels=groupNames) + xlab("Samples") + ggtitle(paste("Red Channel", title))

#     2.7 Sex Prediction

predictedSex <- getSex(GRSetSq)$predictedSex
targets$predSex <- predictedSex
nrow(paired)
write.csv(paired, "same_individuals.csv")
for (row in paired){
  print(row)
}

f <- function(x, output){
  var1 <- x[1]
  idx1 <- which(targets[,] == var1)
  sex1 <- targets[idx1,]$predSex
  gen1 <- targets[idx1,]$genealogy
  cgen1 <- targets[idx1,]$cgenalogy
  methid1 <- targets[idx1,]$meth_id
  bl <- c()
  var2 <- x[2]
  idx2 <- which(targets[,] == var2)
  sex2 <- targets[idx2,]$predSex
  gen2 <- targets[idx2,]$genealogy
  cgen2 <- targets[idx2,]$cgenalogy
  methid2 <- targets[idx1,]$meth_id
  
  cat(paste(var1, sex1,gen1, cgen1, methid1, bl, var2, sex2, gen2, cgen2, methid2, sep = ","), file=output, append = T, fill = T)
}

apply(paired, 1, f, output = "paired.csv")



# f <- function(x, output){
#   for (i in 1:2){
#     var1 <- x[i]
#     idx1 <- which(targets[,] == var1)
#     sex1 <- targets[idx1,]$predSex
#     gen1 <- targets[idx1,]$genealogy
#     cgen1 <- targets[idx1,]$cgenalogy
#     methid1 <- targets[idx1,]$meth_id
#     cat(paste(var1, sex1,gen1, cgen1, methid1, sep = ","), file=output, append = T, fill = T)
#   }
#   cat(paste(c(), sep = ","), file=output, append = T, fill = T)
# }
# 
# apply(paired, 1, f, output = "paired_for.csv")


#     2.8 Report All QC pdf
qcReport(rgSet,  sampNames = targets$Sample_Name, 
         sampGroups = factor(targets$Sample_Group), 
         pdf = file.path(results_path, "qcReport.pdf"))

# 3. Normalization
#     3.1 Normalization method: ssNoob 

norm_method <- "ssNoob"
if (norm_method == "Quantile"){
  GRSetSq <- preprocessQuantile(rgSet)
  
} else if (norm_method == "Funnorm") {
  GRSetSq <- preprocessFunnorm(rgSet)
  
} else if (norm_method == "SWAN") {
  mSetSq <- preprocessSWAN(rgSet)
  # Convert to RatioSet object.
  RatioSetSq <- ratioConvert(mSetSq)
  # Convert to GenomicRatioSet object.
  GRSetSq <- mapToGenome(RatioSetSq)
} else if (norm_method == "ssNoob") {
  mSetSq <- preprocessNoob(rgSet)
  # Convert to RatioSet object.
  RatioSetSq <- ratioConvert(mSetSq)
  # Convert to GenomicRatioSet object.
  GRSetSq <- mapToGenome(RatioSetSq)
} else if (norm_method == "Illumina") {
  mSetSq <- preprocessIllumina(rgSet, bg.correct = TRUE, normalize = "controls",
                               reference = 1)
  # Convert to RatioSet object.
  RatioSetSq <- ratioConvert(mSetSq)
  # Convert to GenomicRatioSet object.
  GRSetSq <- mapToGenome(RatioSetSq)
} else {
  stop("[ERROR] normalization method not correctly specified.")
}

GRSetSq

bVals_normalized <- getBeta(GRSetSq)
png(filename = "normalized_bvals_densityPlot.png")
densityPlot(bVals_normalized, sampGroups = targets$Sample_Group, legend=FALSE, main="Normalized")
legend("top", legend = levels(factor(targets$Sample_Group)), text.col = brewer.pal(8,"Dark2"))
dev.off()

mVals_normalized <- getM(GRSetSq)
png(filename = "normalized_mvals_densityPlot.png")
densityPlot(mVals_normalized, sampGroups = targets$Sample_Group, legend=FALSE, main="Normalized", xlab = "M values")
legend("top", legend = levels(factor(targets$Sample_Group)), text.col = brewer.pal(8,"Dark2"))
dev.off()


write.csv(bVals_normalized, file = "normalized_bVals.csv")
write.csv(mVals_normalized, file = "normalized_mVals.csv")

############################## Predicted Sex

predictedSex <- getSex(GRSetSq)
same_patients_idx
plotSex(addSex(GRSetSq))
rownames(predictedSex) <- targets$Sample_Name
df_eq_patients = predictedSex[same_patients_idx, ]

plot(yMed ~ xMed, 
     data = df_eq_patients,
     pch="",
     ylab = "Y chr; median total intensity (log2)",
     xlab = "X chr; median total intensity (log2)", 
     col = as.numeric(factor(df_eq_patients$predictedSex)))
text(yMed ~ xMed, 
     data = df_eq_patients, 
     labels = rownames(df_eq_patients),
     col = as.numeric(factor(df_eq_patients$predictedSex)), cex=0.5)
legend("topleft", 
       legend = levels(factor(df_eq_patients$predictedSex)),
       text.col = labels(factor(df_eq_patients$predictedSex)))



# 4. Filtering
#     4.1 Remove failed probes (CpG; or mark as NA)
detP <- detectionP(rgSet)
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(GRSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)
# keep
# FALSE   TRUE 
# 44026 821833

GRSetSqFlt <- GRSetSq[keep,]
bVals <- bVals[keep, ]
mVals <- mVals[keep, ]

#     4.2 Remove probes on the sex chromosomes

# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(GRSetSqFlt) %in% annot$Name[annot$chr %in% 
                                                     c("chrX","chrY")])
table(keep)
GRSetSqFlt <- GRSetSqFlt[keep,]

#     4.3 Remove probes with SNPs at CpG site and at a single nucleotide extension

# remove probes with SNPs at CpG site
GRSetSqFlt <- dropLociWithSnps(GRSetSqFlt)
GRSetSqFlt

#     4.4 Exclude cross reactive probes
"/mnt/ElRaid/amontalban/PROJECTS/methylation/illumina450k_filtering-master/"
# exclude cross reactive probes 
xReactiveProbes <- read.csv(file=paste("/mnt/ElRaid/amontalban/PROJECTS/methylation/illumina450k_filtering-master",
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)
xReactiveProbes <- read.csv(file=paste("/mnt/ElRaid/amontalban/PROJECTS/methylation/illumina450k_filtering-master",
                                       "EPIC/13059_2016_1066_MOESM1_ESM.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(GRSetSqFlt) %in% xReactiveProbes$TargetID)

table(keep)
GRSetSqFlt <- GRSetSqFlt[keep,]

#  4.5 Optional (?) Filtering by variance

# 5. Data Exploration
#     5.1 MDS

# MDS - raw data
# Beta-values
# M-values
png(filename = "MDSplot_raw.png", height = 1000,
    width = 2000,
    res = 200)
par(mfrow=c(1,2))
pal <- brewer.pal(8,"Dark2")
plotMDS(getM(mset), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], cex=0.8, pch=19, main="M values")
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.65, bg="white")
plotMDS(getBeta(mset), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], cex=0.8, pch=19)
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.65, bg="white")
dev.off()


# MDS - normalized
png(filename = "MDS_plot_normalized.png", height = 1000,
width = 2000,
res = 200)
par(mfrow=c(1,2))
pal <- brewer.pal(8,"Dark2")
plotMDS(getM(GRSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], cex=0.8, pch=19)
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.65, bg="white")

plotMDS(getBeta(GRSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], cex=0.8, pch=19)
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")
dev.off()

# MDS - filtering
png(filename = "MDSplot_filtered.png", height = 1000,
    width = 2000,
    res = 200)
par(mfrow=c(1,2))
pal <- brewer.pal(8,"Dark2")
plotMDS(getM(GRSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], cex=0.8, pch=19)
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.65, bg="white")

plotMDS(getBeta(GRSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], cex=0.8, pch=19)
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")
dev.off()

var1_idx <- c(paired$Var1)
var2_idx <- c(paired$Var2)

# Create replcates
for (i in 1:nrow(paired)){
  
}

rep1 <- GRSetSq[,c(9,98)]
pca_rep1 <- prcomp(getM(rep1))
pca_rep1$sdev


#     5.2 t-SNE

# t-sne (of 50000 random CpGs)
set.seed(1234)
bVals_tsne_subset <- sample(which(complete.cases(bVals_normalized)), 50000)
tsne <- Rtsne(t(bVals_normalized[bVals_tsne_subset, ]), 
              dims = 2, 
              perplexity=6, 
              verbose=TRUE, 
              max_iter = 5000)
png(filename = file.path(results_path, "tsne_by_sample_group.png"),
    res = 300,
    height = 2700,
    width = 2700)
plot(tsne$Y, t='p', pch = 19, main="tsne", col=pal[factor(targets$Sample_Group)],
     xlab = "First dimension", ylab = "Second dimension")
legend("bottomright", legend=levels(factor(targets$Sample_Group)), 
       text.col=pal[factor(targets$Sample_Group)], cex=0.7, bg="white")
dev.off()

#     5.3 UMAP
require(umap)
?umap
bVals_umap_subset <- sample(which(complete.cases(bVals_normalized)), 50000)
x <- t(bVals_normalized[bVals_umap_subset, ])
class(x)
u = umap(x)
plot(u$layout)
head(u$layout)


#     5.4 Examine higher dimensions
png(filename = "MDS_higherDim.png", height = 1000,
    width = 2000,
    res = 200)
par(mfrow=c(1,3))
# Examine higher dimensions to look at other sources of variation
plotMDS(getM(GRSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(1,3), pch=19)
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(GRSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(2,3), pch=19)
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(GRSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(3,4), pch=19)
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")
dev.off()

# Calculte beta-values and m-values filtered
# calculate M-values for statistical analysis
mVals <- getM(GRSetSqFlt)
head(mVals[,1:5])

bVals <- getBeta(GRSetSqFlt)
head(bVals[,1:5])

png(filename = "m_beta_values_densityPlot_filtered.png", height = 1000,
    width = 2000,
    res = 200)
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()

# 7. Probe-wise differential methylation analysis
#     7.1 Paired-test check

contr <- "Responder-NonResponder"
# this is the factor of interest
factorOfInterest <- factor(targets$Sample_Group)
# this is the individual effect that we need to account for
#individual <- factor(targets$Sample_Source) 

# use the above to create a design matrix
design <- model.matrix(~0+factorOfInterest, data=targets)
colnames(design) <- levels(factorOfInterest)
# fit the linear model 
fit <- lmFit(mVals, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(contr,
                            levels=design)
contMatrix
factorOfInterest <- factor(targets$Sample_Group)
design <- model.matrix(~0+factorOfInterest, data = targets)
colnames(design) <- levels(factorOfInterest)
fit <- lmFit(mVals, design)
print(unique(factorOfInterest)[[1]])
contr <- paste(as.character(unique(factorOfInterest)[2]), as.character(unique(factorOfInterest)[1]), sep = "-")
print(contr)
contMatrix <- makeContrasts(contrasts = contr, levels = design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
annSub <- annot[match(rownames(mVals), annot$Name), c(1:4, 22:24)]
DMPs <- topTable(fit2, num=Inf, coef = 1, genelist = annSub)
summary(decideTests(fit2))
# > summary(decideTests(fit2))
# Responder-NonResponder
# Down                    53749
# NotSig                 703881
# Up                      21027

#     7.2 Table with the differential methylated CpGs (DMPs)
DMPs <- topTable(fit2, num=Inf, coef = 1, genelist = annSub)
head(DMPs)
dim(DMPs)
#[1] 778657     13
write.table(DMPs, 
            file = file.path(
              results_path, 
              paste0("DMPs_", contr[1], ".csv")), 
            sep=",", 
            row.names=FALSE)

# plot the top 4 most significantly differentially methylated CpGs (Optional).
pdf(file = file.path(results_path, "Top_4_cpgs.pdf"))
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bVals, 
          cpg = cpg, 
          pheno = factor(targets$Sample_Group), 
          ylab = "Beta values")
})
dev.off()
#     7.3 Regions differentially methyl ated (DMRs)

myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "Responder-NonResponder", arraytype = "EPIC")
#Your contrast returned 74776 individually significant probes. We recommend the default setting of pcutoff in dmrcate().
str(myAnnotation)
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges
# GRanges object with 11693 ranges and 8 metadata columns:
#   seqnames              ranges strand |   no.cpgs     min_smoothed_fdr             Stouffer                HMFDR               Fisher
# <Rle>           <IRanges>  <Rle> | <integer>            <numeric>            <numeric>            <numeric>            <numeric>
#   [1]    chr10 119300783-119306659      * |        42 1.50082861557861e-42 1.72214622000015e-46  0.00143592797806496 1.13026080662784e-40
# [2]     chr6   29585579-29594481      * |        39 5.68435466944129e-33 5.16749569562118e-41 0.000455942927211095 7.42874135638148e-38
# [3]     chr6   29570008-29575145      * |        28 4.00012768795888e-34  4.6884935401062e-41 0.000543452632307404 1.33740939380158e-37
# [4]     chr6   33130696-33149047      * |       122 1.26703981822703e-23 2.47510216338683e-38   0.0150505499767834 2.02912983810536e-32
# [5]    chr10     8092775-8097689      * |        57 9.00528748329111e-62 1.79549865165171e-40   0.0177435358020792 6.51578273423315e-32
# ...      ...                 ...    ... .       ...                  ...                  ...                  ...                  ...
# [11689]     chr6   31370532-31371614      * |        30 2.08419011800205e-05    0.306678233336363    0.195250680115013    0.403825601555806
# [11690]    chr18   11908320-11908331      * |         4  0.00037317411112869    0.277602387708461    0.367348423576868    0.453513434860284
# [11691]    chr13   20875915-20876028      * |         2 0.000317042155334718    0.784423309260349     0.31320737710662    0.492559220960415
# [11692]     chr8 126441997-126442060      * |         2 0.000279111550869505    0.801201605160385    0.647051114114112    0.809184196357049
# [11693]    chr19     2783719-2784040      * |         2 0.000343583621180086    0.977032404511928     0.90987560251073    0.984616667113857
# maxdiff              meandiff                              overlapping.genes
# <numeric>             <numeric>                                    <character>
#   [1]    -0.177704110179173    -0.112423686724337                                   EMX2, EMX2OS
# [2]    -0.219420464546638    -0.109988431154654                                GABBR1, SNORA20
# [3]    -0.210726487048355     -0.12629397810295                                GABBR1, SNORA20
# [4]     -0.18698814602606   -0.0607566872588069                      SNORA38, COL11A2, SNORA20
# [5]    -0.177509459612241    -0.095727898283236 RP11-379F12.4, GATA3, GATA3-AS1, RP11-379F12.3
# ...                   ...                   ...                                            ...
# [11689]    -0.033628510027584  -0.00565432061880964                            HCP5, MICA, SNORA20
# [11690]    0.0119258390790291   0.00701114990732444                                          MPPE1
# [11691]   -0.0675533721747508    -0.033425532862982                               SNORA16, SNORD37
# [11692] -0.000406109611949663 -0.000190348769967191                                           <NA>
#   [11693]   -0.0105135956484278  -0.00433851870450386                                           <NA>
#   -------
#   seqinfo: 22 sequences from an unspecified genome; no seqlengths
# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]

# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 2, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "EPIC", genome = "hg19")

myAnnotation_1 <- DMRcate::cpg.annotate(object = mVals, datatype = "array", what = "M", 
                                        analysis.type = "differential", design = design, 
                                        contrasts = TRUE, cont.matrix = contMatrix, 
                                        coef = contr[1], arraytype = "EPIC", fdr = 0.05)

str(myAnnotation_1)


DMRs_1 <- DMRcate::dmrcate(myAnnotation_1, lambda=1000, C=2)

print(DMRs_1)

# convert the regions to annotated genomic ranges
#data(dmrcatedata)  # This was used in old versions of DMRcate
results.ranges.1 <- DMRcate::extractRanges(DMRs_1, genome = "hg19")

# Save results tables
write.csv(results.ranges.1, file = file.path(
  results_path,
  paste0("DMRs_results_", contr[1],".csv")
)
)

# PLOT
# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]
samps <- 1:nrow(targets)

# draw the plot for the top DMR

png(file = file.path(results_path, paste0("Top_DMR_", contr[1], ".png")),
    res = 300,
    height = 6000,
    width = 4000)
DMRcate::DMR.plot(ranges=results.ranges.1, dmr=1, CpGs=bVals, phen.col=cols, 
                  what = "Beta", arraytype = "EPIC", pch=16, toscale=TRUE, 
                  plotmedians=TRUE, genome="hg19", samps=samps)
dev.off()


# 7.4 Heatmap with the significant analysis
sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
print(sigCpGs[1:10])
print(length(sigCpGs))

png(file = file.path(
  results_path, 
  paste0("Heatmap_sigCpGs", contr[1], ".png")),
  width = 2000,
  height = 2000,
  res = 260)
heatmap.2(bVals[sigCpGs[1:1000], ], col = greenred(75), trace = "none", 
          labRow = FALSE, ColSideColors = pal[factor(targets$Sample_Group)])
dev.off()

d3heatmap(bVals[sigCpGs[1:1000], ], col = greenred(75), trace = "none")



# 8. Gene Ontology Enrichment (cluster profiler)

#### prova

patients <- targets$cgenalogy
patients <- substr(patients, 1, 7)
length(patients)
length(unique(patients))
