help("DGEList")
d <- DGEList(counts = mobData,group = factor(mobDataGroups))
d
dim(d)
d.full <- d
d$counts
head(cpm(d))
help("cpm")
apply(d$counts, 2, sum)
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)
head(d)
d$samples$lib.size <- colSums(d$counts)
d$samples
d <- calcNormFactors(d)
d
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)
d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)
# GLM extimates of dispersion
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)
#Comparing Deseq and EdgeR
require(DESeq2)
cds <- newCountDataSet( data.frame(d$counts), d$samples$group )
cds <- DESeqDataSetFromMatrix(
  countData = d$counts,
  colData = data.frame(group = d$samples$group),
  design = ~ group)
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
cds <- estimateDispersions(cds)
plotDispEsts(cds)
#Differential Expression
et12 <- exactTest(d1, pair=c(1,2)) 
et13 <- exactTest(d1, pair=c(1,3))
et23 <- exactTest(d1, pair=c(2,3))
topTags(et12, n=10)
topTags(et13, n=5)
topTags(et23, n=5)
cds <- DESeq(cds)
de1 <- results(cds, contrast = c("group", "MM", "WM"))
summary(de1)
de1 <- rownames(de1) [which (de1$padj < 0.05)]
plotMA(de1)
abline(h = c(-2, 2), col = "blue")
#GLM testing for differential expression
design.mat
fit <- glmFit(d2, design.mat)
lrt12 <- glmLRT(fit, contrast=c(1,-1,0))
lrt13 <- glmLRT(fit, contrast=c(1,0,-1))
lrt23 <- glmLRT(fit, contrast=c(0,1,-1))
topTags(lrt12, n=10)
# RNA sequencing for Cancer research- R section
metadatah<- read.csv(file= "metadata.csv", header=TRUE, sep = ",")
counts <- read.csv("hisat2.csv", header=TRUE, sep = ",")
head(metadatah)
head(counts)
countdata = counts[,c(7:32)]
countdata <- countdata[,-c(2)]
countdata <- countdata[,-c(1)]
rownames(countdata) <- counts[,1]
table(colnames(countdata)==metadatah$SampleID)
dds <- DESeqDataSetFromMatrix(countData = round(countdata), colData = metadatah, design = ~Condition)
nrow(dds)
head(dds)
dds1 <- dds[ rowSums(counts(dds)) >= 10, ]
nrow(dds1)
PCAdata <- prcomp(t(assay(dds1)))
autoplot(PCAdata, data = metadatah,colour = "Condition", label = FALSE)
plotDensity(assay(dds1), col=1:24,lwd=2,lty=1,xlab("Density"),ylab("Counts"))
vst =vst(dds1, blind=FALSE)
PCAdata <- prcomp(t(assay(vst)))
autoplot(PCAdata, data = metadatah,colour = "Condition",label = FALSE, main="PCA")
clusters2 <- hclust(dist( t( assay(vst) ) ),method ="ward.D")
plot(clusters2, main = "Dendrogram", label = metadatah$Condition)
plotDensity(assay(vst), lwd=2,lty=1,xlab("Density"),ylab("Counts"), main = "Density plot")
sampleDists <- dist(t(assay(vst)))
sampleDistMatrix <- as.matrix( sampleDists)
rownames(sampleDistMatrix) <- paste( vst$SampleID,vst$Condition, sep="-" )
colnames(sampleDistMatrix) <- paste( dds1$SampleID,dds1$Condition, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, main = "Sample Distance Matrix ")
# Running DeSeq
dds1 <- DESeq(dds1)
head(assay(dds1))
res_table <- results(dds1)
summary(res_table)
res2 <- results(dds1, alpha=0.05)
sum(res2$padj <0.05, na.rm=TRUE)
res_sig <- subset(res2, padj <0.05)
dim(res_sig)
res_small_p <- res_sig[order(res_sig$pvalue),]
dim(res_small_p)
write.csv(as.data.frame(res_small_p), "DESeq2_hisat.csv")
