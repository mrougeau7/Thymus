setwd("~/GitHub/Thymus/Datasets/Data with 11 cell types/WT1")

library(gplots)

msz=100  #matrix size in hist2d
cluster=4  #number of clusters produced for K means
txtfiles_WT1=list.files(pattern="*.txt")  #loads in order all files within folder

mat_WT1<-matrix(data=NA,nrow=10000,ncol=33) #this will change with size of matrix and # of sections/cell types

toFill_WT1<-as.data.frame(mat_WT1)  #changes to data frame

toUse0_WT1<-seq(1,33,1) #uses all 33 columns
toUse1_WT1<-seq(1,33,3) #uses columns that correspond to Section1
toUse2_WT1<-seq(2,33,3) #Section2
toUse3_WT1<-seq(3,33,3) #Section3

for(i in toUse0_WT1){
  tmp_WT1 = read.table(txtfiles_WT1[i], sep="\t", head=T)
  X_WT1=tmp_WT1[[1]]  #[,1]
  Y_WT1=tmp_WT1[[2]]
  my.xy_WT1<-hist2d(X_WT1,Y_WT1,nbins=c(msz,msz))
  obj_WT1 <- my.xy_WT1$counts
  
  obj_WT1[1,1] = obj_WT1[1,1]-1
  obj_WT1[100,1] = obj_WT1[100,1]-1
  obj_WT1[1,100] = obj_WT1[1,100]-1
  obj_WT1[100,100]=obj_WT1[100,100]-1
  
  count_WT1<- as.vector(obj_WT1)
  toFill_WT1[,i]<-count_WT1
  
  result_WT1<-toFill_WT1[,toUse0_WT1]
  result_WT1->result0_WT1
}

for(i in toUse1_WT1){
  tmp_WT1 = read.table(txtfiles_WT1[i], sep="\t", head=T)
  X_WT1=tmp_WT1[[1]]  #[,1]
  Y_WT1=tmp_WT1[[2]]
  my.xy_WT1<-hist2d(X_WT1,Y_WT1,nbins=c(msz,msz))
  obj_WT1 <- my.xy_WT1$counts
  
  obj_WT1[1,1] = obj_WT1[1,1]-1
  obj_WT1[100,1] = obj_WT1[100,1]-1
  obj_WT1[1,100] = obj_WT1[1,100]-1
  obj_WT1[100,100]=obj_WT1[100,100]-1
  
  count_WT1<- as.vector(obj_WT1)
  toFill_WT1[,i]<-count_WT1
  
  result_WT1<-toFill_WT1[,toUse1_WT1]
  result_WT1->result1_WT1
}

for(i in toUse2_WT1){
  tmp_WT1 = read.table(txtfiles_WT1[i], sep="\t", head=T)
  X_WT1=tmp_WT1[[1]]  #[,1]
  Y_WT1=tmp_WT1[[2]]
  my.xy_WT1<-hist2d(X_WT1,Y_WT1,nbins=c(msz,msz))
  obj_WT1 <- my.xy_WT1$counts
  
  obj_WT1[1,1] = obj_WT1[1,1]-1
  obj_WT1[100,1] = obj_WT1[100,1]-1
  obj_WT1[1,100] = obj_WT1[1,100]-1
  obj_WT1[100,100]=obj_WT1[100,100]-1
  
  count_WT1<- as.vector(obj_WT1)
  toFill_WT1[,i]<-count_WT1
  
  result_WT1<-toFill_WT1[,toUse2_WT1]
  result_WT1->result2_WT1
}

for(i in toUse3_WT1){
  tmp_WT1 = read.table(txtfiles_WT1[i], sep="\t", head=T)
  X_WT1=tmp_WT1[[1]]  #[,1]
  Y_WT1=tmp_WT1[[2]]
  my.xy_WT1<-hist2d(X_WT1,Y_WT1,nbins=c(msz,msz))
  obj_WT1 <- my.xy_WT1$counts
  
  obj_WT1[1,1] = obj_WT1[1,1]-1
  obj_WT1[100,1] = obj_WT1[100,1]-1
  obj_WT1[1,100] = obj_WT1[1,100]-1
  obj_WT1[100,100]=obj_WT1[100,100]-1
  
  count_WT1<- as.vector(obj_WT1)
  toFill_WT1[,i]<-count_WT1
  
  result_WT1<-toFill_WT1[,toUse3_WT1]
  result_WT1->result3_WT1   ###Section3
}


library (vegan)
library(labdsv)

# This step normalizes your data and is optional.
#spe.std <- decostand(spe, "normalize")  #You can also use "standardize". See 'help' for details.

# Do the k-means clustering [read 'help' for the 'kmeans' function to see what the arguments "centers" (clusters or k) and "nstart" (randomizations) mean].

spe.kmeans_All_WT1 <- kmeans(result0_WT1, centers=cluster, nstart=10)
spe.kmeans_S1_WT1 <- kmeans(result1_WT1, centers=cluster, nstart=10)
spe.kmeans_S2_WT1 <- kmeans(result2_WT1, centers=cluster, nstart=10)
spe.kmeans_S3_WT1 <- kmeans(result3_WT1, centers=cluster, nstart=10)

library(gplots) ###Call 1 set at a time; produces image of section
v0_WT1=spe.kmeans_All_WT1$cluster
v1_WT1=spe.kmeans_S1_WT1$cluster  
v2_WT1=spe.kmeans_S2_WT1$cluster
v3_WT1=spe.kmeans_S3_WT1$cluster

#load("thymus.Rdata")
#image(v2) # make pic
v0.1_WT1=matrix(v0_WT1,nrow=100,ncol=100)
v1.1_WT1=matrix(v1_WT1,nrow=100,ncol=100)
v2.1_WT1=matrix(v2_WT1,nrow=100,ncol=100)
v3.1_WT1=matrix(v3_WT1,nrow=100,ncol=100)

#v2_WT1<-v2_WT1[1:40,1:100] # trim
#image(v2_WT1) # make new pic

heatmap.2( v0.1_WT1, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v0.1_WT1,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("red", "yellow", "blue", "green")) 

heatmap.2( v1.1_WT1, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v1.1_WT1,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("red", "yellow", "blue", "green")) 

heatmap.2( v2.1_WT1, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v2.1_WT1,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("red", "yellow", "blue", "green")) 

heatmap.2( v3.1_WT1, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v3.1_WT1,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("red", "yellow", "blue", "green")) 

clust.ext<-v2_WT1[1,1] # get cluster number of corner pixel (should always represent 'exterior')
idx_Wt1<-which(v2_WT1==clust.ext_WT1,arr.ind=T) # get addresses of all pixels of that number
#mat[mat < 0.1] <- NA
idx2_WT1<-which(idx_Wt1[,1]<40) # find such pixels for which there is another pixel immediately to right
idx_WT1<-idx_Wt1[idx2_Wt1,] # use only those pixels
idx_Wt1[,1]<-idx_Wt1[,1]+1 # find pixels immediately to right
z_Wt1<-table(v2_WT1[idx_Wt1]) # make sorted table of cluster numbers of these pixels
clust.int_WT1<-as.integer(which(z_Wt1==max(z_WT1[z_Wt1!=max(z_WT1)]))) # get the cluster number that is 2nd most common (this is the thing most often 'touching' exterior, i.e. interior but not medular)
#mat[mat < 0.1] <- NA

clust.med1_Wt1<-as.integer(which(z_Wt1==min(z_WT1[z_WT1!=min(z_WT1)])))
clust.med2_Wt1<-as.integer(which(z_Wt1==min(z_WT1)))
idx.med1_Wt1<-which(v2_WT1==clust.med1_WT1)
idx.med2_WT1<-which(v2_WT1==clust.med2_Wt1)

####ADD s

mean(toFill[idx.med1,1],na.rm=T) #looking at certain column numbers
mean(toFill[idx.med2,1],na.rm=T)
mean(toFill[idx.med1,22],na.rm=T)
mean(toFill[idx.med2,22],na.rm=T)
mean(toFill[idx.med1,25],na.rm=T)
mean(toFill[idx.med2,25],na.rm=T)

col1.01<-1*(toFill[,1]>0)
col22.01<-1*(toFill[,22]>0)
col25.01<-1*(toFill[,25]>0)

mean(col1.01[idx.med1],na.rm=T)
mean(col1.01[idx.med2],na.rm=T)
mean(col22.01[idx.med1],na.rm=T)
mean(col22.01[idx.med2],na.rm=T)
mean(col25.01[idx.med1],na.rm=T)
mean(col25.01[idx.med2],na.rm=T)

heatmap.2( v2_WT1, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v2_WT1,
notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
margins = c(0,0),col=c("green", "yellow", "blue", "red")) 


my.order<-c(1,2,3,4) # define the order we want to plot panels
par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
for (i in 1:max(spe.kmeans_All_WT1$cluster)){ 
barplot(colSums(result0_WT1[which(spe.kmeans_All_WT1$cluster==i),]),main=i,ylim=c(0,12000)) 
# we pick out desired cluster and plot
}

my.order<-c(1,2,3,4) # define the order we want to plot panels
par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
for (i in 1:max(spe.kmeans_S1_WT1$cluster)){ 
  barplot(colSums(result1_WT1[which(spe.kmeans_S1_WT1$cluster==i),]),main=i,ylim=c(0,12000)) 
 #we pick out desired cluster and plot
}

my.order<-c(1,2,3,4) # define the order we want to plot panels
par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
for (i in 1:max(spe.kmeans_S2_WT1$cluster)){ 
  barplot(colSums(result2_WT1[which(spe.kmeans_S2_WT1$cluster==i),]),main=i,ylim=c(0,12000)) 
# we pick out desired cluster and plot
}

my.order<-c(1,2,3,4) # define the order we want to plot panels
par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
for (i in 1:max(spe.kmeans_S3_WT1$cluster)){ 
  barplot(colSums(result3_WT1[which(spe.kmeans_S3_WT1$cluster==i),]),main=i,ylim=c(0,12000)) 
# we pick out desired cluster and plot
}

b2_WT1<-colSums(result0_WT1[which(spe.kmeans_All_WT1$cluster==2),])
b3_WT1<-colSums(result0_WT1[which(spe.kmeans_All_WT1$cluster==3),])
b4_WT1<-colSums(result0_WT1[which(spe.kmeans_All_WT1$cluster==4),])

b2_S1_WT1<-colSums(result1_WT1[which(spe.kmeans_S1_WT1$cluster==2),])
b3_S1_WT1<-colSums(result1_WT1[which(spe.kmeans_S1_WT1$cluster==3),])
b4_S1_WT1<-colSums(result1_WT1[which(spe.kmeans_S1_WT1$cluster==4),])

b2_S2_WT1<-colSums(result2_WT1[which(spe.kmeans_S2_WT1$cluster==2),])
b3_S2_WT1<-colSums(result2_WT1[which(spe.kmeans_S2_WT1$cluster==3),])
b4_S2_WT1<-colSums(result2_WT1[which(spe.kmeans_S2_WT1$cluster==4),])

b2_S3_WT1<-colSums(result3_WT1[which(spe.kmeans_S3_WT1$cluster==2),])
b3_S3_WT1<-colSums(result3_WT1[which(spe.kmeans_S3_WT1$cluster==3),])
b4_S3_WT1<-colSums(result3_WT1[which(spe.kmeans_S3_WT1$cluster==4),])

###and normalize these so that each cluster type has the same total amount of cells (in other words, we're ###getting the proportion of cell types in each cluster)


b2_WT1<-b2_WT1/sum(b2_WT1)
b3_WT1<-b3_WT1/sum(b3_WT1)
b4_WT1<-b4_WT1/sum(b4_WT1)

b2_S1_WT1<-b2_S1_WT1/sum(b2_S1_WT1)
b3_S1_WT1<-b3_S1_WT1/sum(b3_S1_WT1)
b4_S1_WT1<-b4_S1_WT1/sum(b4_S1_WT1)

b2_S2_WT1<-b2_S2_WT1/sum(b2_S2_WT1)
b3_S2_WT1<-b3_S2_WT1/sum(b3_S2_WT1)
b4_S2_WT1<-b4_S2_WT1/sum(b4_S2_WT1)

b2_S3_WT1<-b2_S3_WT1/sum(b2_S3_WT1)
b3_S3_WT1<-b3_S3_WT1/sum(b3_S3_WT1)
b4_S3_WT1<-b4_S3_WT1/sum(b4_S3_WT1)

#We can look at the difference between pairs of clusters to see which cell types are most different between #clusters

par(mfrow=c(2,2))  ##All
barplot(b2_WT1-b3_WT1,main="2-3")
barplot(b2_WT1-b4_WT1,main="2-4")
barplot(b3_WT1-b4_WT1,main="3-4")

par(mfrow=c(2,2))
barplot(b2_S1_WT1-b3_S1_WT1,main="2-3")
barplot(b2_S1_WT1-b4_S1_WT1,main="2-4")
barplot(b3_S1_WT1-b4_S1_WT1,main="3-4")

par(mfrow=c(2,2))
barplot(b2_S2_WT1-b3_S2_WT1,main="2-3")
barplot(b2_S2_WT1-b4_S2_WT1,main="2-4")
barplot(b3_S2_WT1-b4_S2_WT1,main="3-4")

par(mfrow=c(2,2))
barplot(b2_S3_WT1-b3_S3_WT1,main="2-3")
barplot(b2_S3_WT1-b4_S3_WT1,main="2-4")
barplot(b3_S3_WT1-b4_S3_WT1,main="3-4")

###Bray Curtis

library(vegan) 
library(dplyr) # used for chaining and manipulation  
# http://blog.rstudio.org/2014/01/17/introducing-dplyr/
#reads as: take only cluster 1 data from 'spe', then take column sums (i.e. total cell type count) 
#then transpose (swap rows/cols needed for bray curtis calc)

df_WT1<-result0_WT1[which(spe.kmeans_All_WT1$cluster==1),] %>% colSums() %>% t()   # example of chaining. 
for (j in 2:4){  
  tmp_WT1<-result0_WT1[which(spe.kmeans_All_WT1$cluster==j),] %>% colSums() %>% t() 
  df_WT1<-rbind(df_WT1,tmp_WT1)
}

df_S1_WT1<-result1_WT1[which(spe.kmeans_S1_WT1$cluster==1),] %>% colSums() %>% t()
for (j in 2:4){  
  tmp_S1_WT1<-result1_WT1[which(spe.kmeans_S1_WT1$cluster==j),] %>% colSums() %>% t() 
  df_S1_WT1<-rbind(df_S1_WT1,tmp_S1_WT1)
}

df_S2_WT1<-result2_WT1[which(spe.kmeans_S1_WT1$cluster==1),] %>% colSums() %>% t()
for (j in 2:4){  
  tmp_S2_WT1<-result2_WT1[which(spe.kmeans_S1_WT1$cluster==j),] %>% colSums() %>% t() 
  df_S1_WT1<-rbind(df_S1_WT1,tmp_S1_WT1)
}

df_S3_WT1<-result3_WT1[which(spe.kmeans_S3_WT1$cluster==1),] %>% colSums() %>% t()
for (j in 2:4){  
  tmp_S3_WT1<-result3_WT1[which(spe.kmeans_S3_WT1$cluster==j),] %>% colSums() %>% t() 
  df_S3_WT1<-rbind(df_S3_WT1,tmp_S3_WT1)
}

#need same number of columns to bind them.

#colnames(df2_WT1)<-colnames(df_WT1)
colnames(df_S1_WT1)=colnames(df_S2_WT1)
colnames(df_S3_WT1)<-colnames(df_S2_WT1)

#df_WT1<-apply(df_WT1,2,as.integer) %>% as.data.frame()
#df_S1_WT1<-apply(df_S1_WT1,2,as.integer) %>% as.data.frame()
#df_S2_WT1<-apply(df_S2_WT1,2,as.integer) %>% as.data.frame()
#df_S3_WT1<-apply(df_S3_WT1,2,as.integer) %>% as.data.frame()


#dfTOTAL<-rbind(df_S1_WT1,df_S2_WT1, df_S3_WT1, df_S1_WT2, df_S2_WT2, df_S3_WT3, df_S1_Mutant1,
#df_S2_Mutant1, df_S3_Mutant1, df_S1_Mutant2, df_S2_Mutant2, df_S3_Mutant3)
dfT=rbind(df_S1_WT1,df_S2_WT1,df_S3_WT1)

BrayCurtis<-vegdist(dfT,method="bray")
#print(BrayCurtis)
hc<-hclust(BrayCurtis)
#plot(hc,labels=dfT$rownames)
plot(hc)

save.image("WT1.rdata")