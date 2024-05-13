library(affy) ### Load your library every time
library(GEOquery)
library(limma)
library(WGCNA)
library(sva)
library(biomaRt)
library(Biobase)

setwd("C:/Users/hanzo/OneDrive/Desktop/Microarray Course/Microarray course")

###### Steps in Analysis
### (1) Get Data
### (2) Quality Control (QC) on Raw Data 
### (3) Normalization
### (4) Batch Correction
### (5) Outlier Removal
### (6) QC on Normalized Data
### (7) Covariate Analysis
### (8) Annotate Probes
### (9) Collapse Rows
### (10) Differential Expression Analysis

###### (1) Get Data

## Two Choices For Getting Data:
## (a) Load your data directly from a source into your R session (using GEOquery R library)
## (b) Download your data to your local computer or server before loading it into your R session

### (a) Load your data directly from a source into your R session (using GEOquery R library)

## GOAL: to get datMeta and datExpr for prefrontal cortex samples only 

# getGEO() to get phenotype data for entire experiment - datMeta
gse <- getGEO("GSE20295", GSEMatrix =TRUE,getGPL=FALSE)

datMeta <- pData(gse[[1]])
rownames(datMeta) <- datMeta[,2]
datMeta$title <- gsub(" ","_",datMeta$title)
idx <- which(datMeta$source_name_ch1 == "Postmortem brain prefrontal cortex")

# getGEOSuppFiles() to get .CEL files to create affy object with expression data for entire experiment - datExpr
getGEOSuppFiles("GSE20295")
# Use software such as 7-zip to extract into "GSE20295_RAW"
# manually change GSM506039_1_1364_BA9_Pm.CEL.gz to GSM506039_1364_BA9_Pm.CEL.gz first - this was an upload error
filesPFC <- paste(datMeta$geo_accession,"_",datMeta$title,".CEL.gz",sep="")[idx]
data.affy <- ReadAffy(celfile.path = "./GSE20295/GSE20295_RAW", filenames = filesPFC)
datExpr <- exprs(data.affy)
datMeta <- datMeta[idx,]

# check ordering
GSM <- rownames(pData(data.affy))
GSM <- substr(GSM,1,9)
idx <- match(GSM, datMeta$geo_accession)
datMeta <- datMeta[idx,] 

# reformat datMeta and datExpr
datMeta <- datMeta[,-c(3:7,14:36)]
colnames(datMeta)[5:8] <- c("Dx","Sex","Age","Region")  
datMeta$Dx <- gsub("disease state: control","CTL",datMeta$Dx)
datMeta$Dx <- gsub("disease state: Parkinson's disease","PKD",datMeta$Dx)
datMeta$Dx <- as.factor(datMeta$Dx)
datMeta$Sex <- gsub("gender: male","M",datMeta$Sex)
datMeta$Sex <- gsub("gender: female","F",datMeta$Sex)
datMeta$Age <- gsub("age: ","",datMeta$Age)
datMeta$Region <- gsub("brain region: ","",datMeta$Region)

### (b) Download your data to your local computer or server before loading it into your R session

## Although option (a) is more straightforward, sometimes with large datasets (> 500 Mb) it is better to 
## 'pre-filter', or only extract/download data that you need - in this experiment, we only want prefrontal cortex.
## Additionally, if you have your own experiment, you can simply use ReadAffy() to load your own expression data.

## You can download the prefrontal cortex expression data directly from the GEO page, 
## (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20295). For GSE20295_RAW.tar choose 'custom' download, and only
## select files with 'BA9' and '.CEL.gz' in the name. Download these files into a new directory 'GSE20295_PFC' 
## and load it into R from there (remember to extract the raw data first with 7-zip).

## GOAL: to get datMeta and datExpr for prefrontal cortex samples only 

# getGEO() to get phenotype data for entire experiment - datMeta
gse <- getGEO("GSE20295",GSEMatrix =TRUE,getGPL=FALSE)

datMeta <- pData(gse[[1]])
rownames(datMeta) <- datMeta[,2]
datMeta$title <- gsub(" ","_",datMeta$title)
idx <- which(datMeta$source_name_ch1 == "Postmortem brain prefrontal cortex")

# manually change GSM506039_1_1364_BA9_Pm.CEL.gz to GSM506039_1364_BA9_Pm.CEL.gz first - this was an upload error
data.affy <- ReadAffy(celfile.path = "./GSE20295_PFC/GSE20295_RAW")
datExpr <- exprs(data.affy)
datMeta <- datMeta[idx,]

# check ordering
GSM <- rownames(pData(data.affy))
GSM <- substr(GSM,1,9)
idx <- match(GSM, datMeta$geo_accession)
datMeta <- datMeta[idx,] 
colnames(datExpr)=rownames(datMeta)

# reformat datMeta and datExpr
datMeta <- datMeta[,-c(3:7,14:36)]
colnames(datMeta)[5:8] <- c("Dx","Sex","Age","Region")  
datMeta$Dx <- gsub("disease state: control","CTL",datMeta$Dx)
datMeta$Dx <- gsub("disease state: Parkinson's disease","PKD",datMeta$Dx)
datMeta$Dx <- as.factor(datMeta$Dx)
datMeta$Sex <- gsub("gender: male","M",datMeta$Sex)
datMeta$Sex <- gsub("gender: female","F",datMeta$Sex)
datMeta$Age <- gsub("age: ","",datMeta$Age)
datMeta$Region <- gsub("brain region: ","",datMeta$Region)

###### (2) Quality Control (QC) on Raw Data 

##primary analysis - examine distribution of expression data across samples

datExpr <- log2(datExpr)
dim(datExpr)

pdf(file="PreQC_Plots.pdf")

#boxplot

boxplot(datExpr,range=0, col = as.numeric(datMeta$Dx), xaxt='n', xlab = "Array", main = "Boxplot", ylab = "Intensity")
legend("topright",legend = levels(datMeta$Dx),fill = as.numeric(as.factor(levels(datMeta$Dx))))

#histogram

i=1; plot(density((datExpr[,i]),na.rm=T),col = as.numeric(datMeta$Dx)[i],
     main = "Histogram", xlab="log2 exp",xlim=c(4,16),ylim=c(0,0.5))
for(i in 2:dim(datExpr)[2]){
  lines(density((datExpr[,i]),na.rm=T), col = as.numeric(datMeta$Dx)[i],)
}
legend("topright",legend = levels(datMeta$Dx),fill = as.numeric(as.factor(levels(datMeta$Dx))))

#mdsplot

mds = cmdscale(dist(t(datExpr)),eig=TRUE)
plot(mds$points,col=as.numeric(datMeta$Dx),pch=19,main="MDS")
legend("topright",legend = levels(datMeta$Dx),fill = as.numeric(as.factor(levels(datMeta$Dx))))

dev.off()

###### (3) Normalization

datExpr <- rma(data.affy, background=T, normalize=T, verbose=T)
datExpr <- exprs(datExpr)

###### (4) Batch Correction

pdf("BatchPlots.pdf")

batch <- protocolData(data.affy)$ScanDate
batch <- substr(batch,1,8)
batch <- as.factor(batch)
table(batch)
datMeta$Batch <- batch

plot(mds$points,col = as.numeric(datMeta$Batch),pch=19,main="MDS Batch")
legend("topright",legend = levels(datMeta$Batch),fill = as.numeric(as.factor(levels(datMeta$Batch))))

# create a datAll expression set object which stores both datMeta and datExpr
# this object will be useful going forward when we have to remove samples from our analysis

datMeta_proc <- new("AnnotatedDataFrame",data=datMeta)
colnames(datExpr) <- rownames(datMeta)

datAll <- new("ExpressionSet",exprs=datExpr,phenoData=datMeta_proc)
rm(datExpr,datMeta,datMeta_proc)  ## we don't need them anymore!

# remove singular batches
to_remove <- (pData(datAll)$Batch == "09/04/03")
datAll <- datAll[,!to_remove]
pData(datAll)$Batch <- droplevels(pData(datAll)$Batch)

# remove batch contributions from expression
mod <- model.matrix(~pData(datAll)$Dx)   
batch <- as.factor(pData(datAll)$Batch)
datExpr.combat <- ComBat(dat = exprs(datAll),batch = batch,mod = mod)

exprs(datAll) <- datExpr.combat

mds = cmdscale(dist(t(exprs(datAll))),eig=TRUE)
plot(mds$points,col=as.numeric(pData(datAll)$Dx),pch=19,main="MDS Diagnosis Post-ComBat")
legend("topright",legend = levels(pData(datAll)$Dx),fill = as.numeric(as.factor(levels(pData(datAll)$Dx))))

plot(mds$points,col=as.numeric(pData(datAll)$Batch),pch=19,main="MDS Batch Post-ComBat")
legend("topright",legend = levels(pData(datAll)$Batch),fill = as.numeric(as.factor(levels(pData(datAll)$Batch))))

dev.off()

###### (5) Outlier Removal

pdf(file="Outliers.pdf")

tree <- hclust(dist(t(exprs(datAll))), method="average")
plot(tree)

normadj <- (0.5 + 0.5*bicor(exprs(datAll)))^2
netsummary <- fundamentalNetworkConcepts(normadj)
C <- netsummary$Connectivity
Z.C <- (C-mean(C))/sqrt(var(C))

datLabel <- pData(datAll)$Dx
plot(1:length(Z.C),Z.C,main="Outlier Plot",xlab = "Samples",ylab="Connectivity Z Score")
text(1:length(Z.C),Z.C,label=datLabel,pos=3,cex=0.6)
abline(h= -2, col="red")

to_keep <- abs(Z.C) < 2
table(to_keep)
colnames(exprs(datAll))[!to_keep]

datAll <- datAll[,to_keep]

dev.off()

###### (6) QC on Normalized Data

dim(datAll)

pdf(file="PostQC_Plots.pdf")

#boxplot

boxplot(exprs(datAll),range=0, col = as.numeric(pData(datAll)$Dx), xaxt='n', xlab = "Array", main = "Boxplot", ylab = "Intensity")
legend("topright",legend = levels(pData(datAll)$Dx),fill = as.numeric(as.factor(levels(pData(datAll)$Dx))))

#histogram

i=1; plot(density((exprs(datAll)[,i]),na.rm=T),col = as.numeric(pData(datAll)$Dx[i]),
          main = "Histogram", xlab="log2 exp")
for(i in 2:dim(exprs(datAll))[2]){
  lines(density((exprs(datAll)[,i]),na.rm=T), col = as.numeric(pData(datAll)$Dx)[i],)
}
legend("topright",legend = levels(pData(datAll)$Dx),fill = as.numeric(as.factor(levels(pData(datAll)$Dx))))

#mdsplot

mds = cmdscale(dist(t(exprs(datAll))),eig=TRUE)
plot(mds$points,col=as.numeric(pData(datAll)$Dx),pch=19,main="MDS")
legend("topright",legend = levels(pData(datAll)$Dx),fill = as.numeric(as.factor(levels(pData(datAll)$Dx))))

dev.off()

###### (7) Covariate Analysis

pdf(file="CovariatesANOVA.pdf")

plot(pData(datAll)$Dx, ylab="Number", main="Subjects")

for(i in c("Sex","Age","Batch")){
  
  if( i == "Sex" || i == "Batch" ){
    print(paste(i,"Character Graph",sep=" "))
    A = anova(lm(as.numeric(as.factor(pData(datAll)[,i])) ~ pData(datAll)$Dx)); p = A$"Pr(>F)"[1]   
    plot(as.factor(pData(datAll)[,i]) ~ pData(datAll)$Dx, main=paste(i," p=", signif(p,2)), ylab="", xlab="")
  }
  else{
    print(paste(i,"Number Graph",sep=" "))
    A = anova(lm(as.numeric(pData(datAll)[,i]) ~ pData(datAll)$Dx)); p = A$"Pr(>F)"[1]   
    plot(as.numeric(as.character(pData(datAll)[,i])) ~ pData(datAll)$Dx, main=paste(i," p=", signif(p,2)), ylab="", xlab="")
  }
}

dev.off()

###### (8) Annotate Probes

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")

f <- listFilters(ensembl)
a <- listAttributes(ensembl)

identifier <- "affy_hg_u133a"
getinfo <- c("affy_hg_u133a", "ensembl_gene_id", "entrezgene", "external_gene_name")
geneDat <- getBM(attributes = getinfo, filters=identifier, values = rownames(exprs(datAll)),mart=ensembl)

idx <- match(rownames(exprs(datAll)),geneDat$affy_hg_u133a)
geneDat <- geneDat[idx,]

table(is.na(geneDat$ensembl_gene_id))

to_keep <- (is.na(geneDat$ensembl_gene_id) == FALSE)
geneDat <- geneDat[to_keep,]
datAll <- datAll[to_keep,]

dim(datAll)
dim(geneDat)

###### (9) Collapse Rows

table(duplicated(geneDat$affy_hg_u133a))
table(duplicated(geneDat$ensembl_gene_id))

CR <- collapseRows(exprs(datAll), rowGroup = geneDat$ensembl_gene_id, rowID = geneDat$affy_hg_u133a)
exprs(datAll) <- CR$datETcollapsed
idx <- match(CR$group2row[,"selectedRowID"], geneDat$affy_hg_u133a)
geneDat <- geneDat[idx,]
rownames(geneDat) <- geneDat$ensembl_gene_id

dim(datAll)
dim(geneDat)

write.csv(exprs(datAll), file = "datExpr.csv")
write.csv(pData(datAll), file = "datMeta.csv")
write.csv(geneDat, file = "geneDat.csv")
save(datAll,file="datAll.RData")

datExpr<-read.csv("datExpr.csv",row.names=1)
datMeta<-read.csv("datMeta.csv",row.names=1)
geneDat<-read.csv("geneDat.csv")
load(file="datAll.RData")

###### (10) Differential Expression Analysis

pData(datAll)$Sex<-factor(pData(datAll)$Sex,levels=c("M","F"))
mod <- model.matrix(~pData(datAll)$Dx+pData(datAll)$Sex+as.numeric(pData(datAll)$Age))
fit <- lmFit(exprs(datAll),mod)
fit <- eBayes(fit)
tt <- topTable(fit,coef = 2,n = Inf,genelist = geneDat)

pdf(file="PvalueDistribution.pdf")
hist(tt$P.Value)
dev.off()


