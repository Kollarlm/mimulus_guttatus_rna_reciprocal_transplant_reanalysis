# Set working directory and load necessary packages
setwd("/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis")

## Modified from Eleanore Ritters

library("DESeq2")
library( "gplots" )
library( "RColorBrewer" )
library("genefilter")

##################################### DIAGNOSTIC PLOTS FOR COL#####################################
##Preparing Col files for DESeq 
ff <- list.files( path = "./data/S1", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
counts.files <- lapply( ff, read.table, skip = 4 )
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 4 ] ) )
ff <- gsub( "ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "[.]/star/", "", ff )
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1
colData <- read.csv(file="matrix.csv")
colData$genotype <- as.factor(colData$genotype)

##Making Col DESeq Data Set
dds = DESeqDataSetFromMatrix(countData = counts, colData = colData,  design = ~genotype + treatment + genotype:treatment)
#Make sure Control is the first level of treatment
dds$treatment <- relevel(dds$treatment, ref = "Control")
#Make sure that dds is properly formatted
as.data.frame(colData(dds))

##Run DESeq
dds <- DESeq(dds)
res <- results(dds)

##Basic diagnostic plots
plotMA( res, ylim = c(-1, 1) )
plotDispEsts( dds, ylim = c(1e-6, 1e1) )
hist( res$pvalue, breaks=20, col="grey" )

##Filtering
# create bins using the quantile function
qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
# "cut" the genes into the bins
bins <- cut( res$baseMean, qs )
# rename the levels of the bins using the middle point
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
# calculate the ratio of ?p? values less than .01 for each bin
ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
# plot these ratios
barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")

##rlog transformation and more diagnostic plots
rld <- rlog( dds )
head( assay(rld) )
par( mfrow = c( 1, 2 ) )
plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )
sampleDists <- dist( t( assay(rld) ) )
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$treatment,
                                     rld$patient, sep="-" )
colnames(sampleDistMatrix) <- NULL
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

plotPCA(rld, intgroup = c("treatment"))
plotPCA(rld, intgroup = c("genotype"))

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
