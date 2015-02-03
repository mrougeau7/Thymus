setwd("~/GitHub/Thymus/Datasets/Data with 11 cell types/WT2")
setwd("~/Thymus/Datasets/Data with 11 cell types/WT2") 

load("WT2.rdata")
load("WT2_100.rdata")
load("WT2_150.rdata")
library(gplots)

msz=200  #matrix size in hist2d
cluster=4  #number of clusters produced for K means
txtfiles_WT2=list.files(pattern="*.txt")  #loads in order all files within folder

mat_WT2<-matrix(data=NA,nrow=40000,ncol=33) #this will change with size of matrix and # of sections/cell types

toFill_WT2<-as.data.frame(mat_WT2)  #changes to data frame

toUse0_WT2<-seq(1,33,1) #uses all 33 columns
toUse1_WT2<-seq(1,33,3) #uses columns that correspond to Section1
toUse2_WT2<-seq(2,33,3) #Section2
toUse3_WT2<-seq(3,33,3) #Section3

for(i in toUse0_WT2){
  tmp_WT2 = read.table(txtfiles_WT2[i], sep="\t", head=T)
  X_WT2=tmp_WT2[[1]]  #[,1]
  Y_WT2=tmp_WT2[[2]]
  my.xy_WT2<-hist2d(X_WT2,Y_WT2,nbins=c(msz,msz))
  obj_WT2 <- my.xy_WT2$counts
  
  obj_WT2[1,1] = obj_WT2[1,1]-1
  obj_WT2[200,1] = obj_WT2[200,1]-1
  obj_WT2[1,200] = obj_WT2[1,200]-1
  obj_WT2[200,200]=obj_WT2[200,200]-1
  
  count_WT2<- as.vector(obj_WT2)
  toFill_WT2[,i]<-count_WT2
  
  result_WT2<-toFill_WT2[,toUse0_WT2]
  result_WT2->result0_WT2
}

for(i in toUse1_WT2){
  tmp_WT2 = read.table(txtfiles_WT2[i], sep="\t", head=T)
  X_WT2=tmp_WT2[[1]]  #[,1]
  Y_WT2=tmp_WT2[[2]]
  my.xy_WT2<-hist2d(X_WT2,Y_WT2,nbins=c(msz,msz))
  obj_WT2 <- my.xy_WT2$counts
  
  obj_WT2[1,1] = obj_WT2[1,1]-1
  obj_WT2[200,1] = obj_WT2[200,1]-1
  obj_WT2[1,200] = obj_WT2[1,200]-1
  obj_WT2[200,200]=obj_WT2[200,200]-1
  
  count_WT2<- as.vector(obj_WT2)
  toFill_WT2[,i]<-count_WT2
  
  result_WT2<-toFill_WT2[,toUse1_WT2]
  result_WT2->result1_WT2
}

for(i in toUse2_WT2){
  tmp_WT2 = read.table(txtfiles_WT2[i], sep="\t", head=T)
  X_WT2=tmp_WT2[[1]]  #[,1]
  Y_WT2=tmp_WT2[[2]]
  my.xy_WT2<-hist2d(X_WT2,Y_WT2,nbins=c(msz,msz))
  obj_WT2 <- my.xy_WT2$counts
  
  obj_WT2[1,1] = obj_WT2[1,1]-1
  obj_WT2[200,1] = obj_WT2[200,1]-1
  obj_WT2[1,200] = obj_WT2[1,200]-1
  obj_WT2[200,200]=obj_WT2[200,200]-1
  
  count_WT2<- as.vector(obj_WT2)
  toFill_WT2[,i]<-count_WT2
  
  result_WT2<-toFill_WT2[,toUse2_WT2]
  result_WT2->result2_WT2
}

for(i in toUse3_WT2){
  tmp_WT2 = read.table(txtfiles_WT2[i], sep="\t", head=T)
  X_WT2=tmp_WT2[[1]]  #[,1]
  Y_WT2=tmp_WT2[[2]]
  my.xy_WT2<-hist2d(X_WT2,Y_WT2,nbins=c(msz,msz))
  obj_WT2 <- my.xy_WT2$counts
  
  obj_WT2[1,1] = obj_WT2[1,1]-1
  obj_WT2[200,1] = obj_WT2[200,1]-1
  obj_WT2[1,200] = obj_WT2[1,200]-1
  obj_WT2[200,200]=obj_WT2[200,200]-1
  
  count_WT2<- as.vector(obj_WT2)
  toFill_WT2[,i]<-count_WT2
  
  result_WT2<-toFill_WT2[,toUse3_WT2]
  result_WT2->result3_WT2   ###Section3
}


library (vegan)
library(labdsv)

# This step normalizes your data and is optional.
#spe.std <- decostand(spe, "normalize")  #You can also use "standardize". See 'help' for details.

# Do the k-means clustering [read 'help' for the 'kmeans' function to see what the arguments "centers" (clusters or k) and "nstart" (randomizations) mean].

spe.kmeans_All_WT2 <- kmeans(result0_WT2, centers=cluster, nstart=1000)
spe.kmeans_S1_WT2 <- kmeans(result1_WT2, centers=cluster, nstart=1000)
spe.kmeans_S2_WT2 <- kmeans(result2_WT2, centers=cluster, nstart=1000)
spe.kmeans_S3_WT2 <- kmeans(result3_WT2, centers=cluster, nstart=1000)

library(gplots) ###Call 1 set at a time; produces image of section
v0_WT2=spe.kmeans_All_WT2$cluster
v1_WT2=spe.kmeans_S1_WT2$cluster  
v2_WT2=spe.kmeans_S2_WT2$cluster
v3_WT2=spe.kmeans_S3_WT2$cluster

#load("thymus.Rdata")
#image(v2) # make pic
v0.1_WT2=matrix(v0_WT2,nrow=200,ncol=200)
v1.1_WT2=matrix(v1_WT2,nrow=200,ncol=200)
v2.1_WT2=matrix(v2_WT2,nrow=200,ncol=200)
v3.1_WT2=matrix(v3_WT2,nrow=200,ncol=200)

#v2_WT2<-v2_WT2[1:40,1:200] # trim
#image(v2_WT2) # make new pic

heatmap.2( v0.1_WT2, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v0.1_WT2,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("red", "yellow", "blue", "green")) 

heatmap.2( v1.1_WT2, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v1.1_WT2,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("red", "yellow", "blue", "green")) 

heatmap.2( v2.1_WT2, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v2.1_WT2,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("red", "yellow", "blue", "green")) 

heatmap.2( v3.1_WT2, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v3.1_WT2,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("red", "yellow", "blue", "green")) 


my.order<-c(1,2,3,4) # define the order we want to plot panels
par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
for (i in 1:max(spe.kmeans_All_WT2$cluster)){ 
  barplot(colSums(result0_WT2[which(spe.kmeans_All_WT2$cluster==i),]),main=i,ylim=c(0,12000)) 
  # we pick out desired cluster and plot
}

my.order<-c(1,2,3,4) # define the order we want to plot panels
par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
for (i in 1:max(spe.kmeans_S1_WT2$cluster)){ 
  barplot(colSums(result1_WT2[which(spe.kmeans_S1_WT2$cluster==i),]),main=i,ylim=c(0,12000)) 
  #we pick out desired cluster and plot
}

my.order<-c(1,2,3,4) # define the order we want to plot panels
par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
for (i in 1:max(spe.kmeans_S2_WT2$cluster)){ 
  barplot(colSums(result2_WT2[which(spe.kmeans_S2_WT2$cluster==i),]),main=i,ylim=c(0,12000)) 
  # we pick out desired cluster and plot
}

my.order<-c(1,2,3,4) # define the order we want to plot panels
par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
for (i in 1:max(spe.kmeans_S3_WT2$cluster)){ 
  barplot(colSums(result3_WT2[which(spe.kmeans_S3_WT2$cluster==i),]),main=i,ylim=c(0,12000)) 
  # we pick out desired cluster and plot
}

b1_WT2<-colSums(result0_WT2[which(spe.kmeans_All_WT2$cluster==1),])
b2_WT2<-colSums(result0_WT2[which(spe.kmeans_All_WT2$cluster==2),])
b3_WT2<-colSums(result0_WT2[which(spe.kmeans_All_WT2$cluster==3),])
b4_WT2<-colSums(result0_WT2[which(spe.kmeans_All_WT2$cluster==4),])

b1_S1_WT2<-colSums(result1_WT2[which(spe.kmeans_S1_WT2$cluster==1),])
b2_S1_WT2<-colSums(result1_WT2[which(spe.kmeans_S1_WT2$cluster==2),])
b3_S1_WT2<-colSums(result1_WT2[which(spe.kmeans_S1_WT2$cluster==3),])
b4_S1_WT2<-colSums(result1_WT2[which(spe.kmeans_S1_WT2$cluster==4),])

b1_S2_WT2<-colSums(result2_WT2[which(spe.kmeans_S2_WT2$cluster==1),])
b2_S2_WT2<-colSums(result2_WT2[which(spe.kmeans_S2_WT2$cluster==2),])
b3_S2_WT2<-colSums(result2_WT2[which(spe.kmeans_S2_WT2$cluster==3),])
b4_S2_WT2<-colSums(result2_WT2[which(spe.kmeans_S2_WT2$cluster==4),])

b1_S3_WT2<-colSums(result3_WT2[which(spe.kmeans_S3_WT2$cluster==1),])
b2_S3_WT2<-colSums(result3_WT2[which(spe.kmeans_S3_WT2$cluster==2),])
b3_S3_WT2<-colSums(result3_WT2[which(spe.kmeans_S3_WT2$cluster==3),])
b4_S3_WT2<-colSums(result3_WT2[which(spe.kmeans_S3_WT2$cluster==4),])

###and normalize these so that each cluster type has the same total amount of cells (in other words, we're ###getting the proportion of cell types in each cluster)

b1_WT2<-b1_WT2/sum(b1_WT2)
b2_WT2<-b2_WT2/sum(b2_WT2)
b3_WT2<-b3_WT2/sum(b3_WT2)
b4_WT2<-b4_WT2/sum(b4_WT2)

b1_S1_WT2<-b1_S1_WT2/sum(b1_S1_WT2)
b2_S1_WT2<-b2_S1_WT2/sum(b2_S1_WT2)
b3_S1_WT2<-b3_S1_WT2/sum(b3_S1_WT2)
b4_S1_WT2<-b4_S1_WT2/sum(b4_S1_WT2)

b1_S2_WT2<-b1_S2_WT2/sum(b1_S2_WT2)
b2_S2_WT2<-b2_S2_WT2/sum(b2_S2_WT2)
b3_S2_WT2<-b3_S2_WT2/sum(b3_S2_WT2)
b4_S2_WT2<-b4_S2_WT2/sum(b4_S2_WT2)

b1_S3_WT2<-b1_S3_WT2/sum(b1_S3_WT2)
b2_S3_WT2<-b2_S3_WT2/sum(b2_S3_WT2)
b3_S3_WT2<-b3_S3_WT2/sum(b3_S3_WT2)
b4_S3_WT2<-b4_S3_WT2/sum(b4_S3_WT2)

#We can look at the difference between pairs of clusters to see which cell types are most different between #clusters

par(mfrow=c(2,2))  ##All
barplot(b1_WT2-b2_WT2,main="1-2")
barplot(b1_WT2-b3_WT2,main="1-3")
#barplot(b1_WT2-b4_WT2,main="1-4")
barplot(b2_WT2-b3_WT2,main="2-3")
#barplot(b2_WT2-b4_WT2,main="2-4")
#barplot(b3_WT2-b4_WT2,main="3-4")

par(mfrow=c(2,2))
barplot(b1_S1_WT2-b2_S1_WT2,main="1-2")
barplot(b1_S1_WT2-b3_S1_WT2,main="1-3")
#barplot(b1_S1_WT2-b4_S1_WT2,main="1-4")
barplot(b2_S1_WT2-b3_S1_WT2,main="2-3")
#barplot(b2_S1_WT2-b4_S1_WT2,main="2-4")
#barplot(b3_S1_WT2-b4_S1_WT2,main="3-4")

par(mfrow=c(2,2))
#barplot(b1_S2_WT2-b2_S2_WT2,main="1-2")
barplot(b1_S2_WT2-b3_S2_WT2,main="1-3")
barplot(b1_S2_WT2-b4_S2_WT2,main="1-4")
#barplot(b2_S2_WT2-b3_S2_WT2,main="2-3")
#barplot(b2_S2_WT2-b4_S2_WT2,main="2-4")
barplot(b3_S2_WT2-b4_S2_WT2,main="3-4")

par(mfrow=c(2,2))
barplot(b1_S3_WT2-b2_S3_WT2,main="1-2")
#barplot(b1_S3_WT2-b3_S3_WT2,main="1-3")
barplot(b1_S3_WT2-b4_S3_WT2,main="1-4")
#barplot(b2_S3_WT2-b3_S3_WT2,main="2-3")
barplot(b2_S3_WT2-b4_S3_WT2,main="2-4")
#barplot(b3_S3_WT2-b4_S3_WT2,main="3-4")

###Bray Curtis

library(vegan) 
library(dplyr) # used for chaining and manipulation  
# http://blog.rstudio.org/2014/01/17/introducing-dplyr/
#reads as: take only cluster 1 data from 'spe', then take column sums (i.e. total cell type count) 
#then transpose (swap rows/cols needed for bray curtis calc)

df_WT2<-result0_WT2[which(spe.kmeans_All_WT2$cluster==1),] %>% colSums() %>% t()   # example of chaining. 
for (j in 2:4){  
  tmp_WT2<-result0_WT2[which(spe.kmeans_All_WT2$cluster==j),] %>% colSums() %>% t() 
  df_WT2<-rbind(df_WT2,tmp_WT2)
}

df_S1_WT2<-result1_WT2[which(spe.kmeans_S1_WT2$cluster==1),] %>% colSums() %>% t()
for (j in 2:4){  
  tmp_S1_WT2<-result1_WT2[which(spe.kmeans_S1_WT2$cluster==j),] %>% colSums() %>% t() 
  df_S1_WT2<-rbind(df_S1_WT2,tmp_S1_WT2)
}

df_S2_WT2<-result2_WT2[which(spe.kmeans_S1_WT2$cluster==1),] %>% colSums() %>% t()
for (j in 2:4){  
  tmp_S2_WT2<-result2_WT2[which(spe.kmeans_S1_WT2$cluster==j),] %>% colSums() %>% t() 
  df_S1_WT2<-rbind(df_S1_WT2,tmp_S1_WT2)
}

df_S3_WT2<-result3_WT2[which(spe.kmeans_S3_WT2$cluster==1),] %>% colSums() %>% t()
for (j in 2:4){  
  tmp_S3_WT2<-result3_WT2[which(spe.kmeans_S3_WT2$cluster==j),] %>% colSums() %>% t() 
  df_S3_WT2<-rbind(df_S3_WT2,tmp_S3_WT2)
}

#need same number of columns to bind them.

#colnames(df2_WT2)<-colnames(df_WT2)
colnames(df_S2_WT2)=colnames(df_S1_WT2)
colnames(df_S3_WT2)<-colnames(df_S1_WT2)

#df_WT2<-apply(df_WT2,2,as.integer) %>% as.data.frame()
#df_S1_WT2<-apply(df_S1_WT2,2,as.integer) %>% as.data.frame()
#df_S2_WT2<-apply(df_S2_WT2,2,as.integer) %>% as.data.frame()
#df_S3_WT2<-apply(df_S3_WT2,2,as.integer) %>% as.data.frame()


#dfTOTAL<-rbind(df_S1_WT2,df_S2_WT2, df_S3_WT2, df_S1_WT2, df_S2_WT2, df_S3_WT3, df_S1_Mutant1,
#df_S2_Mutant1, df_S3_Mutant1, df_S1_WT2, df_S2_WT2, df_S3_Mutant3)
dfT_WT2=rbind(df_S1_WT2,df_S2_WT2,df_S3_WT2)

BrayCurtis_WT2<-vegdist(dfT_WT2,method="bray")
#print(BrayCurtis)
hc_WT2<-hclust(BrayCurtis_WT2)
#plot(hc,labels=dfT$rownames)
plot(hc_WT2)

save.image("WT2_150.rdata")
