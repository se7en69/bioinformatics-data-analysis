# Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35570

library(gplots)
library(multtest)
library(MASS)

#read in data annotation files
#setwd('path/to/files') #add path to directory with files and uncomment
setwd('C:/Users/hanzo/OneDrive/Desktop/Microarray Course')
dat <- read.delim('GSE35570_series_matrix.txt',header=T,row.names=1,skip=62)
dat <- head(dat,-1)
ann <- read.table('GSE35570_ann.txt',header=F,row.names=1)

#calculate pearson correlation values
pearson.cor <- cor(dat,method='pearson')

cols <- rev(colorpanel(25,'orange','black','blue')) #set colors to be used in matrix
#heatmap of correlation matrix
layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1),14,1,byrow=TRUE)) #format matrix image
par(oma=c(5,7,1,1),mar=c(8,8,4,1)) #set margin parameters
image(pearson.cor,main='Correlation Matrix, Papillary Thyroid Cancer Samples',
      axes=F,col=cols) #create image of colored plot
axis(1,at=seq(0,1,length=ncol(pearson.cor)),label=dimnames(pearson.cor)[[2]],
     cex.axis=0.9,las=2) #create x axis
axis(2,at=seq(0,1,length=ncol(pearson.cor)),label=dimnames(pearson.cor)[[2]],
     cex.axis=0.9,las=2) #create y axis
title(xlab='Samples',line=6)#adjust axis titles to not overlap sample labels
title(ylab='Samples',line=6)
###Outliers: GSM870796,GSM870800,GSM870804,GSM870852,GSM870877

#CV vs mean plot to identify any additional potential outliers
dat.mean <- apply(log2(dat),2,mean)
dat.sd <- sqrt(apply(log2(dat),2,var))
dat.cv <- dat.sd/dat.mean

plot(dat.mean,dat.cv,main='Papillary Thyroid Cancer Recurrence Dataset\nSample CV vs Mean',
     xlab='Mean',ylab='CV',type='n')
points(dat.mean,dat.cv,col='blue',pch=16)
text(dat.mean,dat.cv,label=dimnames(dat)[[2]],pos=3,cex=0.8)
###Outliers: GSM870784

#remove outlier samples
dat <- subset(dat, select=-c(GSM870784,GSM870796,GSM870800,GSM870804,GSM870852,
                             GSM870877))
#remove outliers from annotation
ann <- ann[!row.names(ann) %in% c('GSM870784','GSM870796','GSM870800','GSM870804',
                                  'GSM870852','GSM870877'),,drop=F]
colnames(ann)[1] <- 'Class'

#CV vs Mean plot of genes pre-filter
gene.mean <- apply(log2(dat),1,mean)
gene.sd <- sqrt(apply(log2(dat),1,var))
gene.cv <- gene.sd/gene.mean

par(oma=c(1,1,1,1))
plot(gene.mean,gene.cv,main='Papillary Thyroid Cancer Dataset\nGene CV vs Mean',
     xlab='Mean',ylab='CV',type='n')
points(gene.mean,gene.cv,col='blue',pch=16)

#filter out genes in bottom quantile of means
q <- quantile(gene.mean)
genes.top <- gene.mean[gene.mean >= q['25%']]
dat.filter <- dat[names(genes.top),]

#3-level ANOVA to return p-values
aov.all.genes <- function(x,s1,s2,s3) {
  x1 <- as.numeric(x[s1])
  x2 <- as.numeric(x[s2])
  x3 <- as.numeric(x[s3])
  fac <- c(rep('A',length(x1)), rep('B',length(x2)), rep('C',length(x3)))
  a.dat <- data.frame(as.factor(fac),c(x1,x2,x3))
  names(a.dat) <- c('factor','express')
  p.out <- summary(aov(express~factor, a.dat))[[1]][1,5]
  #p.out <- summary(aov(express~factor, a.dat))[[1]][1,4]	# use to get F-statistic
  return(p.out)
}

dat.p <- apply(dat.filter,1, aov.all.genes,s1=ann$Class=='NRE',
                 s2=ann$Class=='RE',s3=ann$Class=='N')
adj.p <- mt.rawp2adjp(dat.p,proc='Holm') #holms p-value adjustment
row.names(adj.p$adjp) <- names(dat.p[adj.p$index])

#extract significant genes using p < 0.001 threshold
sig.genes <- adj.p$adjp[adj.p$adjp[,2] < 0.001,2]
#subset data matrix by significant genes only
dat <- dat.filter[names(sig.genes),]

#Principal component analysis, plot first two eigenvectors
dat.pca <- prcomp(t(dat))
plot(range(dat.pca$x[,1]),range(dat.pca$x[,2]),main='PCA Plot of Papillary Thyroid Cancer
     PC1 vs PC2',type='n',xlab='pc1',ylab='pc2')
points(dat.pca$x[,1][ann$Class=='NRE'],dat.pca$x[,2][ann$Class=='NRE'],col='red',
       pch=16)
points(dat.pca$x[,1][ann$Class=='RE'],dat.pca$x[,2][ann$Class=='RE'],col='blue',
       pch=16)
points(dat.pca$x[,1][ann$Class=='N'],dat.pca$x[,2][ann$Class=='N'],col='green',
       pch=16)
legend('topright',legend=c('Sporadic Cancer','Radiation Exposure',
                           'Normal'), pch=16,col=c('red','blue','green'))

#visually shows the expression differences of significantly differentially
#expressed genes between normal vs tumor samples
clrs <- c("#FF0000","#CC0000","#990000","#660000","#330000","#000000","#000000",
    "#0A3300","#146600","#1F9900","#29CC00","#33FF00")
heatmap(as.matrix(dat),col=clrs,main='Heatmap of Significant Genes vs Samples',
             xlab='Sample Classes',ylab='Genes',labCol=ann$Class)
legend('topleft',c('Upregulated','Downregulated'),fill=c('#33FF00','#FF0000'))

#now I want to find the most up/down-regulated genes between classes
#calculate log2 fold change of each gene between sample classes
dat.log <- log2(dat)
norm <- row.names(ann)[which(ann$Class == 'N')]
no.rad <- row.names(ann)[which(ann$Class == 'NRE')]
rad.exp <- row.names(ann)[which(ann$Class == 'RE')]
cancer <- row.names(ann)[which(ann$Class %in% c('RE','NRE'))]

norm.mean <- apply(dat.log[,norm],1,mean,na.rm=T)
no.rad.mean <- apply(dat.log[,no.rad],1,mean,na.rm=T)
rad.exp.mean <- apply(dat.log[,rad.exp],1,mean,na.rm=T)
cancer.mean <- apply(dat.log[,cancer],1,mean,na.rm=T)

radexp.norad.fold <- rad.exp.mean - no.rad.mean #rad exp vs no rad exp
norm.cancer.fold <- norm.mean - cancer.mean #normal vs all cancer samples

#cancer + radiation exposure vs cancer with no radiation exposure
radexp.norad.sorted <- sort(radexp.norad.fold,decreasing=T)
#print 5 most up- and down-regulated genes b/w rad. exp. & non-rad. exp.
print(head(radexp.norad.sorted))
print(tail(radexp.norad.sorted))

#volcano plot, expression fold changes b/w NRE vs RE tumor samples
p.trans <- -1 * log10(dat.p)
par(oma=c(3,1,3,1),mfrow=c(1,1))
plot(range(p.trans),range(radexp.norad.fold),type='n',xlab='-1log10(p-values)',
     ylab='Fold Change, log2 Scale',
     main='Volcano Plot, Rad Exp vs No Rad Samples Fold Change')
points(p.trans[names(radexp.norad.fold)],radexp.norad.fold,pch=16)

#volcano plot to find largest fold changes between normal & cancerous cells
plot(range(p.trans),range(norm.cancer.fold),type='n',xlab='-1log10(p-values)',
     ylab='Fold Change, log2 Scale',
     main='Volcano Plot, Normal vs Tumor Samples Fold Change')
points(p.trans[names(norm.cancer.fold)],norm.cancer.fold,pch=16)
abline(h=log2(2)) #lines at +/- 2-fold changes
abline(h=-log2(2))

#sort to find min and max fold changes
#normal vs cancer cells
norm.cancer.sorted <- sort(norm.cancer.fold,decreasing=T)
#print 5 most up- and down- regulated genes b/w cancer cells and normal cells
cat('Top downregulated genes in cancer cells:\n',head(names(norm.cancer.sorted),n=5))
cat('Top upregulated genes in cancer cells:\n',tail(names(norm.cancer.sorted),n=5))



##additional dimensionality reduction if necessary
#classification using linear discriminant analysis
#train set: first 10 normal, first 8 NRE, first 8 RE
dat.trans <- t(dat)
dat.trans <- data.frame(ann,dat.trans)
dat.train <- rbind(dat.trans[which(dat.trans$Class=='N')[1:10],],
                   dat.trans[which(dat.trans$Class=='RE')[1:8],],
                   dat.trans[which(dat.trans$Class=='NRE')[1:8],])
dat.train <- as.data.frame(dat.train)
#test set: all samples, remove class labels (already stored in ann)
dat.trans <- dat.trans[,-1]

#run LDA using all genes
train.lda <- lda(dat.train$Class~.,dat.train)
dat.pred <- predict(train.lda,dat.trans)
table(dat.pred$class,ann$Class)

#return point color based on cluster assignment
class.colors.lda <- function(x) {
  if(dat.pred$class[x] == 'N'){
    return('red')
  }
  else if(dat.pred$class[x] == 'NRE') {
    return('blue')
  }
  else if(dat.pred$class[x] =='RE') {
    return('green')
  }
}
#plot first 2 discriminant functions against each other
par(mfrow=c(1,1))
plot(dat.pred$x,pch=16,xlab='Discriminant Function 1',type='n',
     ylab='Discriminant Function 2',main="Discriminant Functions for PTC dataset")
for (i in 1:nrow(dat.pred$x)){
  points(dat.pred$x[i,1],dat.pred$x[i,2],col=class.colors.lda(i),pch=class.point(i))
}
legend('top',c('Group 1','Group 2','Group 3','Sporadic Cancer',
                    'Radiation Exposure','Normal Tissue'),col=cols,pch=symbols)

