setwd("~/GitHub/Thymus/Datasets/Data with 11 cell types/Mutant2")

library(gplots)

msz=100  #matrix size in hist2d
cluster=4  #number of clusters produced for K means
txtfiles_Mutant2=list.files(pattern="*.txt")  #loads in order all files within folder

mat_Mutant2<-matrix(data=NA,nrow=10000,ncol=33) #this will change with size of matrix and # of sections/cell types

toFill_Mutant2<-as.data.frame(mat_Mutant2)  #changes to data frame

toUse0_Mutant2<-seq(1,33,1) #uses all 33 columns
toUse1_Mutant2<-seq(1,33,3) #uses columns that correspond to Section1
toUse2_Mutant2<-seq(2,33,3) #Section2
toUse3_Mutant2<-seq(3,33,3) #Section3

for(i in toUse0_Mutant2){
  tmp_Mutant2 = read.table(txtfiles_Mutant2[i], sep="\t", head=T)
  X_Mutant2=tmp_Mutant2[[1]]  #[,1]
  Y_Mutant2=tmp_Mutant2[[2]]
  my.xy_Mutant2<-hist2d(X_Mutant2,Y_Mutant2,nbins=c(msz,msz))
  obj_Mutant2 <- my.xy_Mutant2$counts
  
  obj_Mutant2[1,1] = obj_Mutant2[1,1]-1
  obj_Mutant2[100,1] = obj_Mutant2[100,1]-1
  obj_Mutant2[1,100] = obj_Mutant2[1,100]-1
  obj_Mutant2[100,100]=obj_Mutant2[100,100]-1
  
  count_Mutant2<- as.vector(obj_Mutant2)
  toFill_Mutant2[,i]<-count_Mutant2
  
  result_Mutant2<-toFill_Mutant2[,toUse0_Mutant2]
  result_Mutant2->result0_Mutant2
}

for(i in toUse1_Mutant2){
  tmp_Mutant2 = read.table(txtfiles_Mutant2[i], sep="\t", head=T)
  X_Mutant2=tmp_Mutant2[[1]]  #[,1]
  Y_Mutant2=tmp_Mutant2[[2]]
  my.xy_Mutant2<-hist2d(X_Mutant2,Y_Mutant2,nbins=c(msz,msz))
  obj_Mutant2 <- my.xy_Mutant2$counts
  
  obj_Mutant2[1,1] = obj_Mutant2[1,1]-1
  obj_Mutant2[100,1] = obj_Mutant2[100,1]-1
  obj_Mutant2[1,100] = obj_Mutant2[1,100]-1
  obj_Mutant2[100,100]=obj_Mutant2[100,100]-1
  
  count_Mutant2<- as.vector(obj_Mutant2)
  toFill_Mutant2[,i]<-count_Mutant2
  
  result_Mutant2<-toFill_Mutant2[,toUse1_Mutant2]
  result_Mutant2->result1_Mutant2
}

for(i in toUse2_Mutant2){
  tmp_Mutant2 = read.table(txtfiles_Mutant2[i], sep="\t", head=T)
  X_Mutant2=tmp_Mutant2[[1]]  #[,1]
  Y_Mutant2=tmp_Mutant2[[2]]
  my.xy_Mutant2<-hist2d(X_Mutant2,Y_Mutant2,nbins=c(msz,msz))
  obj_Mutant2 <- my.xy_Mutant2$counts
  
  obj_Mutant2[1,1] = obj_Mutant2[1,1]-1
  obj_Mutant2[100,1] = obj_Mutant2[100,1]-1
  obj_Mutant2[1,100] = obj_Mutant2[1,100]-1
  obj_Mutant2[100,100]=obj_Mutant2[100,100]-1
  
  count_Mutant2<- as.vector(obj_Mutant2)
  toFill_Mutant2[,i]<-count_Mutant2
  
  result_Mutant2<-toFill_Mutant2[,toUse2_Mutant2]
  result_Mutant2->result2_Mutant2
}

for(i in toUse3_Mutant2){
  tmp_Mutant2 = read.table(txtfiles_Mutant2[i], sep="\t", head=T)
  X_Mutant2=tmp_Mutant2[[1]]  #[,1]
  Y_Mutant2=tmp_Mutant2[[2]]
  my.xy_Mutant2<-hist2d(X_Mutant2,Y_Mutant2,nbins=c(msz,msz))
  obj_Mutant2 <- my.xy_Mutant2$counts
  
  obj_Mutant2[1,1] = obj_Mutant2[1,1]-1
  obj_Mutant2[100,1] = obj_Mutant2[100,1]-1
  obj_Mutant2[1,100] = obj_Mutant2[1,100]-1
  obj_Mutant2[100,100]=obj_Mutant2[100,100]-1
  
  count_Mutant2<- as.vector(obj_Mutant2)
  toFill_Mutant2[,i]<-count_Mutant2
  
  result_Mutant2<-toFill[,toUse3_Mutant2]
  result_Mutant2->result3_Mutant2   ###Section3
}


library (vegan)
library(labdsv)

# This step normalizes your data and is optional.
#spe.std <- decostand(spe, "normalize")  #You can also use "standardize". See 'help' for details.

# Do the k-means clustering [read 'help' for the 'kmeans' function to see what the arguments "centers" (clusters or k) and "nstart" (randomizations) mean].

spe.kmeans_Mutant1_All_Mutant2 <- kmeans(result0_Mutant2, centers=cluster, nstart=10)
spe.kmeans_Mutant1_S1_Mutant2 <- kmeans(result1_Mutant2, centers=cluster, nstart=10)
spe.kmeans_Mutant1_S2_Mutant2 <- kmeans(result2_Mutant2, centers=cluster, nstart=10)
spe.kmeans_Mutant1_S3_Mutant2 <- kmeans(result3_Mutant2, centers=cluster, nstart=10)

library(gplots) ###Call 1 set at a time; produces image of section

v1_Mutant2 <- spe.kmeans_Mutant1_S1_Mutant2$cluster  #change this to section you are interested in looking at

#load("thymus.Rdata")
#image(v2) # make pic
v2_Mutant2 = matrix(v1_Mutant2, nrow = 100, ncol=100)
v2_Mutant2<-v2_Mutant2[1:40,1:100] # trim
image(v2_Mutant2) # make new pic

heatmap.2( v2_Mutant2, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v2_Mutant2,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("green", "yellow", "blue", "red")) 

clust.ext<-v2_Mutant2[1,1] # get cluster number of corner pixel (should always represent 'exterior')
idx_Mutant2<-which(v2_Mutant2==clust.ext_Mutant2,arr.ind=T) # get addresses of all pixels of that number
#mat[mat < 0.1] <- NA
idx2_Mutant2<-which(idx_Mutant2[,1]<40) # find such pixels for which there is another pixel immediately to right
idx_Mutant2<-idx_Mutant2[idx2_Mutant2,] # use only those pixels
idx_Mutant2[,1]<-idx_Mutant2[,1]+1 # find pixels immediately to right
z_Mutant2<-table(v2_Mutant2[idx_Mutant2]) # make sorted table of cluster numbers of these pixels
clust.int_Mutant2<-as.integer(which(z_Mutant2==max(z_Mutant2[z_Mutant2!=max(z_Mutant2)]))) # get the cluster number that is 2nd most common (this is the thing most often 'touching' exterior, i.e. interior but not medular)
#mat[mat < 0.1] <- NA

clust.med1_Mutant2<-as.integer(which(z_Mutant2==min(z_Mutant2[z_Mutant2!=min(z_Mutant2)])))
clust.med2_Mutant2<-as.integer(which(z_Mutant2==min(z_Mutant2)))
idx.med1_Mutant2<-which(v2_Mutant2==clust.med1_Mutant2)
idx.med2_Mutant2<-which(v2_Mutant2==clust.med2_Mutant2)

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



heatmap.2( v2_Mutant2, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=v2_Mutant2,
           notecol="black", trace='none', key=FALSE,lwid = c(.01,0.99),lhei = c(.01,.99),
           margins = c(0,0),col=c("green", "yellow", "blue", "red")) 


my.order<-c(1,2,3,4) # define the order we want to plot panels
par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
for (i in 1:max(spe.kmeans_Mutant1_All$cluster)){ 
  barplot(colSums(result0[which(spe.kmeans_Mutant1_All1$cluster==i),]),main=i,ylim=c(0,12000)) 
  # we pick out desired cluster and plot
}

#my.order<-c(1,2,3,4) # define the order we want to plot panels
#par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
#for (i in 1:max(spe.kmeans_Mutant1_S1$cluster)){ 
#  barplot(colSums(result1[which(spe.kmeans_Mutant1_S1$cluster==i),]),main=i,ylim=c(0,12000)) 
# we pick out desired cluster and plot
#}

#my.order<-c(1,2,3,4) # define the order we want to plot panels
#par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
#for (i in 1:max(spe.kmeans_Mutant1_S2$cluster)){ 
#  barplot(colSums(result2[which(spe.kmeans_Mutant1_S2$cluster==i),]),main=i,ylim=c(0,12000)) 
# we pick out desired cluster and plot
#}

#my.order<-c(1,2,3,4) # define the order we want to plot panels
#par(mfrow=c(2,2)) # make 4 subplots in 2x2 style
#for (i in 1:max(spe.kmeans_Mutant1_S3$cluster)){ 
#  barplot(colSums(result3[which(spe.kmeans_Mutant1_S3$cluster==i),]),main=i,ylim=c(0,12000)) 
# we pick out desired cluster and plot
#}

b2<-colSums(result0[which(spe.kmeans_Mutant1_All$cluster==2),])
b3<-colSums(result0[which(spe.kmeans_Mutant1_All$cluster==3),])
b4<-colSums(result0[which(spe.kmeans_Mutant1_All$cluster==4),])

#b2_S1<-colSums(result1[which(spe.kmeans_Mutant1_S1$cluster==2),])
#b3_S1<-colSums(result1[which(spe.kmeans_Mutant1_S1$cluster==3),])
#b4_S1<-colSums(result1[which(spe.kmeans_Mutant1_S1$cluster==4),])

#b2_S2<-colSums(result2[which(spe.kmeans_Mutant1_S2$cluster==2),])
#b3_S2<-colSums(result2[which(spe.kmeans_Mutant1_S2$cluster==3),])
#b4_S2<-colSums(result2[which(spe.kmeans_Mutant1_S2$cluster==4),])

#b2_S3<-colSums(result3[which(spe.kmeans_Mutant1_S3$cluster==2),])
#b3_S3<-colSums(result3[which(spe.kmeans_Mutant1_S3$cluster==3),])
#b4_S3<-colSums(result3[which(spe.kmeans_Mutant1_S3$cluster==4),])

###and normalize these so that each cluster type has the same total amount of cells (in other words, we're ###getting the proportion of cell types in each cluster)


b2<-b2/sum(b2)
b3<-b3/sum(b3)
b4<-b4/sum(b4)

#b2_S1<-b2/sum(b2_S1)
#b3_S1<-b3/sum(b3_S1)
#b4_S1<-b4/sum(b4_S1)

#b2_S2<-b2/sum(b2_S2)
#b3_S2<-b3/sum(b3_S2)
#b4_S2<-b4/sum(b4_S2)

#b2_S3<-b2/sum(b2_S3)
#b3_S3<-b3/sum(b3_S3)
#b4_S3<-b4/sum(b4_S3)

#We can look at the difference between pairs of clusters to see which cell types are most different between #clusters

par(mfrow=c(2,2))  ##All
barplot(b2-b3,main="2-3")
barplot(b2-b4,main="2-4")
barplot(b3-b4,main="3-4")

#par(mfrow=c(2,2))
#barplot(b2-b3,main="2-3")
#barplot(b2-b4,main="2-4")
#barplot(b3-b4,main="3-4")

#par(mfrow=c(2,2))
#barplot(b2-b3,main="2-3")
#barplot(b2-b4,main="2-4")
#barplot(b3-b4,main="3-4")

#par(mfrow=c(2,2))
#barplot(b2-b3,main="2-3")
#barplot(b2-b4,main="2-4")
#barplot(b3-b4,main="3-4")

###Bray Curtis

library(vegan) 
library(dplyr) # used for chaining and manipulation  
# http://blog.rstudio.org/2014/01/17/introducing-dplyr/
#reads as: take only cluster 1 data from 'spe', then take column sums (i.e. total cell type count) 
#then transpose (swap rows/cols needed for bray curtis calc)

df<-result0[which(spe.kmeans_Mutant1_All$cluster==1),] %>% colSums() %>% t()   # example of chaining. 
for (j in 2:4){  
  tmp<-result0[which(spe.kmeans_Mutant1_All$cluster==j),] %>% colSums() %>% t() 
  df<-rbind(df,tmp)
}

#df1.1<-spe1.1[which(spe.kmeans1.1$cluster==1),] %>% colSums() %>% t()
#for (j in 2:4){  
#  tmp1.1<-spe1.1[which(spe.kmeans1.1$cluster==j),] %>% colSums() %>% t() 
#  df1.1<-rbind(df1.1,tmp1.1)
#}

#df1.2<-spe1.2[which(spe.kmeans1.2$cluster==1),] %>% colSums() %>% t()
#for (j in 2:4){  
#  tmp1.2<-spe1.2[which(spe.kmeans1.2$cluster==j),] %>% colSums() %>% t() 
#  df1.2<-rbind(df1.2,tmp1.2)
#}

#df1.3<-spe1.3[which(spe.kmeans1.3$cluster==1),] %>% colSums() %>% t()
#for (j in 2:4){  
#  tmp1.2<-spe1.3[which(spe.kmeans1.3$cluster==j),] %>% colSums() %>% t() 
#  df1.3<-rbind(df1.3,tmp1.3)
#}


df<-apply(df,2,as.integer) %>% as.data.frame()
#df1.1<-apply(df1.1,2,as.integer) %>% as.data.frame()
#df1.2<-apply(df1.2,2,as.integer) %>% as.data.frame()
#df1.3<-apply(df1.3,2,as.integer) %>% as.data.frame()

#need same number of columns to bind them.
colnames(df2)<-colnames(df)
#colnames(df1.1)=colnames(df2.1)
#colnames(df1.2)<-colnames(df2.1)
#colnames(df1.2)<-colnames(df2.1)

#dfTOTAL<-rbind(df1.2,df2.2, df1.2, df1.2, df1.1, df1.1)

BrayCurtis<-vegdist(dfTOTAL,method="bray")
#print(BrayCurtis)
hc<-hclust(BrayCurtis)
plot(hc,labels=dfTOTAL$rownames)

save.image("Mutant2.rdata")