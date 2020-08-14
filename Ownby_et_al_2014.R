### Data sets
### minl <- two-way table containing mineralogical data
### chem <- two-way table containing log10 transformed chemical ppm data

### initialize required packages
library(cluster)
library(vegan)

### create matrix of Gower dissimilarities among samples for mineral data
min.gow <- as.matrix(daisy(minl,metric='gower'))
### standardize by dividing matrix by maximum observed value
min.gow <- min.gow/max(min.gow)

### create matrix of Euclidean distances among samples for chemical data
chem.ec <- as.matrix(dist(chem,method='euclidean'))
### standardize by dividing matrix by maximum observed value
chem.ec <- chem.ec/max(chem.ec)

### set value for Mu
mu <- 0.5

### create combined dissimilarity/distance matrix based on value of Mu
com.diss <- (chem.ec*mu)+(min.gow*(1-mu))

### subject combined matrix to mullti-dimensional scaling (MDS)
mds.results <- cmdscale(com.diss,k=2) ### k <- number of dimensions returned

### plot MDS results and label points
plot(mds.results,xlab='Dim 1',ylab='Dim 2',pch=16)
text(mds.results,labels=row.names(com.diss),pos=1,cex=0.5)

### Calculating Sibson's coefficient for values of Mu and plotting results
### create sequence of Mu values from 0 to 1 by 0.05
mu.list <- seq(from=0,to=1,by=0.05)

### Calculate combined dissimilarity/distance matrices for all values of Mu
com.list <- list()
for (i in 1:length(mu.list)) {
  com.diss <- (chem.ec*mu.list[i])+(min.gow*(1-mu.list[i]))
  com.list [[i]] <- cmdscale(com.diss)}

### Calculate Sibson's coefficient comparing combined dissimilarities to    ### mineral and chemical data for all values of Mu
sib.min <- list()
for(i in 1:length(mu.list)) {
  sib.min[[i]] <- procrustes(cmdscale(min.gow),com.list[[i]],symmetric=T)$ss}

sib.chem <- list()
for(i in 1:length(mu.list)) {
  sib.chem[[i]] <- procrustes(cmdscale(chem.ec),com.list[[i]],symmetric=T)$ss}

### Plot Sibson's coefficients for mineral and chemical data by values of Mu
plot(mu.list,as.vector(sib.chem),type='l',col='red',ylim=c(0,1))
points(mu.list,as.vector(sib.min),type='l',lty=2,col='blue')
