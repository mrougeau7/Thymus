setwd("~/GitHub/Thymus/Datasets/Data with 11 cell types/Mutant2")
setwd("~/Thymus/Datasets/Data with 11 cell types/Mutant2")  
load("Mutant2_Rarefaction_150.rdata")

library(gplots)

msz=150  #matrix size in hist2d
cluster=4  #number of clusters produced for K means
txtfiles_Mutant2_R=list.files(pattern="*.txt")  #loads in order all files within folder

mat_Mutant2_R<-matrix(data=NA,nrow=22500,ncol=33) #this will change with size of matrix and # of sections/cell types

toFill_Mutant2_R<-as.data.frame(mat_Mutant2_R)  #changes to data frame

toUse0_Mutant2_R<-seq(1,33,1) #uses all 33 columns
toUse1_Mutant2_R<-seq(1,33,3) #uses columns that correspond to Section1
toUse2_Mutant2_R<-seq(2,33,3) #Section2
toUse3_Mutant2_R<-seq(3,33,3) #Section3

for(i in toUse0_Mutant2_R){
  tmp_Mutant2_R = read.table(txtfiles_Mutant2_R[i], sep="\t", head=T)
  X_Mutant2_R=tmp_Mutant2_R[[1]]  #[,1]
  Y_Mutant2_R=tmp_Mutant2_R[[2]]
  my.xy_Mutant2_R<-hist2d(X_Mutant2_R,Y_Mutant2_R,nbins=c(msz,msz))
  obj_Mutant2_R <- my.xy_Mutant2_R$counts
  
  obj_Mutant2_R[1,1] = obj_Mutant2_R[1,1]-1
  obj_Mutant2_R[150,1] = obj_Mutant2_R[150,1]-1
  obj_Mutant2_R[1,150] = obj_Mutant2_R[1,150]-1
  obj_Mutant2_R[150,150]=obj_Mutant2_R[150,150]-1
  
  count_Mutant2_R<- as.vector(obj_Mutant2_R)
  toFill_Mutant2_R[,i]<-count_Mutant2_R
  
  result_Mutant2_R<-toFill_Mutant2_R[,toUse0_Mutant2_R]
  result_Mutant2_R->result0_Mutant2_R
}

for(i in toUse1_Mutant2_R){
  tmp_Mutant2_R = read.table(txtfiles_Mutant2_R[i], sep="\t", head=T)
  X_Mutant2_R=tmp_Mutant2_R[[1]]  #[,1]
  Y_Mutant2_R=tmp_Mutant2_R[[2]]
  my.xy_Mutant2_R<-hist2d(X_Mutant2_R,Y_Mutant2_R,nbins=c(msz,msz))
  obj_Mutant2_R <- my.xy_Mutant2_R$counts
  
  obj_Mutant2_R[1,1] = obj_Mutant2_R[1,1]-1
  obj_Mutant2_R[150,1] = obj_Mutant2_R[150,1]-1
  obj_Mutant2_R[1,150] = obj_Mutant2_R[1,150]-1
  obj_Mutant2_R[150,150]=obj_Mutant2_R[150,150]-1
  
  count_Mutant2_R<- as.vector(obj_Mutant2_R)
  toFill_Mutant2_R[,i]<-count_Mutant2_R
  
  result_Mutant2_R<-toFill_Mutant2_R[,toUse1_Mutant2_R]
  result_Mutant2_R->result1_Mutant2_R
}

for(i in toUse2_Mutant2_R){
  tmp_Mutant2_R = read.table(txtfiles_Mutant2_R[i], sep="\t", head=T)
  X_Mutant2_R=tmp_Mutant2_R[[1]]  #[,1]
  Y_Mutant2_R=tmp_Mutant2_R[[2]]
  my.xy_Mutant2_R<-hist2d(X_Mutant2_R,Y_Mutant2_R,nbins=c(msz,msz))
  obj_Mutant2_R <- my.xy_Mutant2_R$counts
  
  obj_Mutant2_R[1,1] = obj_Mutant2_R[1,1]-1
  obj_Mutant2_R[150,1] = obj_Mutant2_R[150,1]-1
  obj_Mutant2_R[1,150] = obj_Mutant2_R[1,150]-1
  obj_Mutant2_R[150,150]=obj_Mutant2_R[150,150]-1
  
  count_Mutant2_R<- as.vector(obj_Mutant2_R)
  toFill_Mutant2_R[,i]<-count_Mutant2_R
  
  result_Mutant2_R<-toFill_Mutant2_R[,toUse2_Mutant2_R]
  result_Mutant2_R->result2_Mutant2_R
}

for(i in toUse3_Mutant2_R){
  tmp_Mutant2_R = read.table(txtfiles_Mutant2_R[i], sep="\t", head=T)
  X_Mutant2_R=tmp_Mutant2_R[[1]]  #[,1]
  Y_Mutant2_R=tmp_Mutant2_R[[2]]
  my.xy_Mutant2_R<-hist2d(X_Mutant2_R,Y_Mutant2_R,nbins=c(msz,msz))
  obj_Mutant2_R <- my.xy_Mutant2_R$counts
  
  obj_Mutant2_R[1,1] = obj_Mutant2_R[1,1]-1
  obj_Mutant2_R[150,1] = obj_Mutant2_R[150,1]-1
  obj_Mutant2_R[1,150] = obj_Mutant2_R[1,150]-1
  obj_Mutant2_R[150,150]=obj_Mutant2_R[150,150]-1
  
  count_Mutant2_R<- as.vector(obj_Mutant2_R)
  toFill_Mutant2_R[,i]<-count_Mutant2_R
  
  result_Mutant2_R<-toFill_Mutant2_R[,toUse3_Mutant2_R]
  result_Mutant2_R->result3_Mutant2_R   ###Section3
}


library (vegan)
library(labdsv)

# This step normalizes your data and is optional.
#spe.std <- decostand(spe, "normalize")  #You can also use "standardize". See 'help' for details.

spe.std0 <- decostand(result0_Mutant2_R, "normalize")  #You can also use "standardize". See 'help' for details.
spe.std1 <- decostand(result1_Mutant2_R, "normalize") 
spe.std2 <- decostand(result2_Mutant2_R, "normalize") 
spe.std3 <- decostand(result3_Mutant2_R, "normalize") 

spe.kmeans_All_Mutant2_R <- kmeans(spe.std0, centers=cluster, nstart=1000)
spe.kmeans_S1_Mutant2_R <- kmeans(spe.std1, centers=cluster, nstart=1000)
spe.kmeans_S2_Mutant2_R <- kmeans(spe.std2, centers=cluster, nstart=1000)
spe.kmeans_S3_Mutant2_R <- kmeans(spe.std3, centers=cluster, nstart=1000)

# Do the k-means clustering [read 'help' for the 'kmeans' function to see what the arguments "centers" (clusters or k) and "nstart" (randomizations) mean].

v0_Mutant2_R=spe.kmeans_All_Mutant2_R$cluster
v1_Mutant2_R=spe.kmeans_S1_Mutant2_R$cluster  
v2_Mutant2_R=spe.kmeans_S2_Mutant2_R$cluster
v3_Mutant2_R=spe.kmeans_S3_Mutant2_R$cluster

v0.1_Mutant2_R=matrix(v0_Mutant2_R,nrow=150,ncol=150)
v1.1_Mutant2_R=matrix(v1_Mutant2_R,nrow=150,ncol=150)
v2.1_Mutant2_R=matrix(v2_Mutant2_R,nrow=150,ncol=150)
v3.1_Mutant2_R=matrix(v3_Mutant2_R,nrow=150,ncol=150)

heatmap.2( v0.1_Mutant2_R, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v0.1_Mutant2_R,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("red", "yellow", "blue", "green")) 

heatmap.2( v1.1_Mutant2_R, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v1.1_Mutant2_R,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("red", "yellow", "blue", "green")) 

heatmap.2( v2.1_Mutant2_R, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v2.1_Mutant2_R,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("red", "yellow", "blue", "green")) 

heatmap.2( v3.1_Mutant2_R, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v3.1_Mutant2_R,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("red", "yellow", "blue", "green")) 

my.order<-c(1,2,3,4) # define the order we want to plot panels
par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
for (i in 1:max(spe.kmeans_All_Mutant2$cluster)){ 
  barplot(colSums(result0_Mutant2[which(spe.kmeans_All_Mutant2$cluster==i),]),main=i,ylim=c(0,11500)) 
  # we pick out desired cluster and plot
}

my.order<-c(1,2,3,4) # define the order we want to plot panels
par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
for (i in 1:max(spe.kmeans_S1_Mutant2$cluster)){ 
  barplot(colSums(result1_Mutant2[which(spe.kmeans_S1_Mutant2$cluster==i),]),main=i,ylim=c(0,11500)) 
  #we pick out desired cluster and plot
}

my.order<-c(1,2,3,4) # define the order we want to plot panels
par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
for (i in 1:max(spe.kmeans_S2_Mutant2$cluster)){ 
  barplot(colSums(result2_Mutant2[which(spe.kmeans_S2_Mutant2$cluster==i),]),main=i,ylim=c(0,11500)) 
  # we pick out desired cluster and plot
}

my.order<-c(1,2,3,4) # define the order we want to plot panels
par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
for (i in 1:max(spe.kmeans_S3_Mutant2$cluster)){ 
  barplot(colSums(result3_Mutant2[which(spe.kmeans_S3_Mutant2$cluster==i),]),main=i,ylim=c(0,11500)) 
  # we pick out desired cluster and plot
}

b1_Mutant2<-colSums(result0_Mutant2[which(spe.kmeans_All_Mutant2$cluster==1),])
b2_Mutant2<-colSums(result0_Mutant2[which(spe.kmeans_All_Mutant2$cluster==2),])
b3_Mutant2<-colSums(result0_Mutant2[which(spe.kmeans_All_Mutant2$cluster==3),])
b4_Mutant2<-colSums(result0_Mutant2[which(spe.kmeans_All_Mutant2$cluster==4),])

b1_S1_Mutant2<-colSums(result1_Mutant2[which(spe.kmeans_S1_Mutant2$cluster==1),])
b2_S1_Mutant2<-colSums(result1_Mutant2[which(spe.kmeans_S1_Mutant2$cluster==2),])
b3_S1_Mutant2<-colSums(result1_Mutant2[which(spe.kmeans_S1_Mutant2$cluster==3),])
b4_S1_Mutant2<-colSums(result1_Mutant2[which(spe.kmeans_S1_Mutant2$cluster==4),])

b1_S2_Mutant2<-colSums(result2_Mutant2[which(spe.kmeans_S2_Mutant2$cluster==1),])
b2_S2_Mutant2<-colSums(result2_Mutant2[which(spe.kmeans_S2_Mutant2$cluster==2),])
b3_S2_Mutant2<-colSums(result2_Mutant2[which(spe.kmeans_S2_Mutant2$cluster==3),])
b4_S2_Mutant2<-colSums(result2_Mutant2[which(spe.kmeans_S2_Mutant2$cluster==4),])

b1_S3_Mutant2<-colSums(result3_Mutant2[which(spe.kmeans_S3_Mutant2$cluster==1),])
b2_S3_Mutant2<-colSums(result3_Mutant2[which(spe.kmeans_S3_Mutant2$cluster==2),])
b3_S3_Mutant2<-colSums(result3_Mutant2[which(spe.kmeans_S3_Mutant2$cluster==3),])
b4_S3_Mutant2<-colSums(result3_Mutant2[which(spe.kmeans_S3_Mutant2$cluster==4),])

###and normalize these so that each cluster type has the same total amount of cells (in other words, we're ###getting the proportion of cell types in each cluster)

b1_Mutant2<-b1_Mutant2/sum(b1_Mutant2)
b2_Mutant2<-b2_Mutant2/sum(b2_Mutant2)
b3_Mutant2<-b3_Mutant2/sum(b3_Mutant2)
b4_Mutant2<-b4_Mutant2/sum(b4_Mutant2)

b1_S1_Mutant2<-b1_S1_Mutant2/sum(b1_S1_Mutant2)
b2_S1_Mutant2<-b2_S1_Mutant2/sum(b2_S1_Mutant2)
b3_S1_Mutant2<-b3_S1_Mutant2/sum(b3_S1_Mutant2)
b4_S1_Mutant2<-b4_S1_Mutant2/sum(b4_S1_Mutant2)

b1_S2_Mutant2<-b1_S2_Mutant2/sum(b1_S2_Mutant2)
b2_S2_Mutant2<-b2_S2_Mutant2/sum(b2_S2_Mutant2)
b3_S2_Mutant2<-b3_S2_Mutant2/sum(b3_S2_Mutant2)
b4_S2_Mutant2<-b4_S2_Mutant2/sum(b4_S2_Mutant2)

b1_S3_Mutant2<-b1_S3_Mutant2/sum(b1_S3_Mutant2)
b2_S3_Mutant2<-b2_S3_Mutant2/sum(b2_S3_Mutant2)
b3_S3_Mutant2<-b3_S3_Mutant2/sum(b3_S3_Mutant2)
b4_S3_Mutant2<-b4_S3_Mutant2/sum(b4_S3_Mutant2)

#We can look at the difference between pairs of clusters to see which cell types are most different between #clusters

par(mfrow=c(3,2))  ##All
barplot(b1_Mutant2-b2_Mutant2,main="1-2")
barplot(b1_Mutant2-b3_Mutant2,main="1-3")
barplot(b1_Mutant2-b4_Mutant2,main="1-4")
barplot(b2_Mutant2-b3_Mutant2,main="2-3")
barplot(b2_Mutant2-b4_Mutant2,main="2-4")
barplot(b3_Mutant2-b4_Mutant2,main="3-4")

par(mfrow=c(3,2))
barplot(b1_S1_Mutant2-b2_S1_Mutant2,main="1-2")
barplot(b1_S1_Mutant2-b3_S1_Mutant2,main="1-3")
barplot(b1_S1_Mutant2-b4_S1_Mutant2,main="1-4")
barplot(b2_S1_Mutant2-b3_S1_Mutant2,main="2-3")
barplot(b2_S1_Mutant2-b4_S1_Mutant2,main="2-4")
barplot(b3_S1_Mutant2-b4_S1_Mutant2,main="3-4")

par(mfrow=c(3,2))
barplot(b1_S2_Mutant2-b2_S2_Mutant2,main="1-2")
barplot(b1_S2_Mutant2-b3_S2_Mutant2,main="1-3")
barplot(b1_S2_Mutant2-b4_S2_Mutant2,main="1-4")
barplot(b2_S2_Mutant2-b3_S2_Mutant2,main="2-3")
barplot(b2_S2_Mutant2-b4_S2_Mutant2,main="2-4")
barplot(b3_S2_Mutant2-b4_S2_Mutant2,main="3-4")

par(mfrow=c(3,2))
barplot(b1_S3_Mutant2-b2_S3_Mutant2,main="1-2")
barplot(b1_S3_Mutant2-b3_S3_Mutant2,main="1-3")
barplot(b1_S3_Mutant2-b4_S3_Mutant2,main="1-4")
barplot(b2_S3_Mutant2-b3_S3_Mutant2,main="2-3")
barplot(b2_S3_Mutant2-b4_S3_Mutant2,main="2-4")
barplot(b3_S3_Mutant2-b4_S3_Mutant2,main="3-4")

###Bray Curtis

library(vegan) 
library(dplyr) # used for chaining and manipulation  
# http://blog.rstudio.org/2014/01/17/introducing-dplyr/
#reads as: take only cluster 1 data from 'spe', then take column sums (i.e. total cell type count) 
#then transpose (swap rows/cols needed for bray curtis calc)

df_Mutant2<-result0_Mutant2[which(spe.kmeans_All_Mutant2$cluster==1),] %>% colSums() %>% t()   # example of chaining. 
for (j in 2:4){  
  tmp_Mutant2<-result0_Mutant2[which(spe.kmeans_All_Mutant2$cluster==j),] %>% colSums() %>% t() 
  df_Mutant2<-rbind(df_Mutant2,tmp_Mutant2)
}

df_S1_Mutant2<-result1_Mutant2[which(spe.kmeans_S1_Mutant2$cluster==1),] %>% colSums() %>% t()
for (j in 2:4){  
  tmp_S1_Mutant2<-result1_Mutant2[which(spe.kmeans_S1_Mutant2$cluster==j),] %>% colSums() %>% t() 
  df_S1_Mutant2<-rbind(df_S1_Mutant2,tmp_S1_Mutant2)
}

df_S2_Mutant2<-result2_Mutant2[which(spe.kmeans_S1_Mutant2$cluster==1),] %>% colSums() %>% t()
for (j in 2:4){  
  tmp_S2_Mutant2<-result2_Mutant2[which(spe.kmeans_S1_Mutant2$cluster==j),] %>% colSums() %>% t() 
  df_S1_Mutant2<-rbind(df_S1_Mutant2,tmp_S1_Mutant2)
}

df_S3_Mutant2<-result3_Mutant2[which(spe.kmeans_S3_Mutant2$cluster==1),] %>% colSums() %>% t()
for (j in 2:4){  
  tmp_S3_Mutant2<-result3_Mutant2[which(spe.kmeans_S3_Mutant2$cluster==j),] %>% colSums() %>% t() 
  df_S3_Mutant2<-rbind(df_S3_Mutant2,tmp_S3_Mutant2)
}

#need same number of columns to bind them.

#colnames(df2_Mutant2)<-colnames(df_Mutant2)
colnames(df_S2_Mutant2)=colnames(df_S1_Mutant2)
colnames(df_S3_Mutant2)<-colnames(df_S1_Mutant2)

#df_Mutant2<-apply(df_Mutant2,2,as.integer) %>% as.data.frame()
#df_S1_Mutant2<-apply(df_S1_Mutant2,2,as.integer) %>% as.data.frame()
#df_S2_Mutant2<-apply(df_S2_Mutant2,2,as.integer) %>% as.data.frame()
#df_S3_Mutant2<-apply(df_S3_Mutant2,2,as.integer) %>% as.data.frame()


#dfTOTAL<-rbind(df_S1_Mutant2,df_S2_Mutant2, df_S3_Mutant2, df_S1_Mutant2, df_S2_Mutant2, df_S3_WT3, df_S1_Mutant1,
#df_S2_Mutant1, df_S3_Mutant1, df_S1_Mutant2, df_S2_Mutant2, df_S3_Mutant3)
dfT=rbind(df_S1_Mutant2,df_S2_Mutant2,df_S3_Mutant2)

BrayCurtis<-vegdist(dfT,method="bray")
#print(BrayCurtis)
hc<-hclust(BrayCurtis)
#plot(hc,labels=dfT$rownames)
plot(hc)

save.image("Mutant2_Rarefaction_150.rdata")
