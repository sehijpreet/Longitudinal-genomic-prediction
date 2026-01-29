rm(list=ls())

library(dplyr)

setwd('')

list.files()

geno <- read.csv('SNP_fil.csv')
head(geno)
dim(geno)

duplicates <- geno %>% 
  group_by(Genotype) %>%
  filter(n() > 1) %>%
  ungroup()

duplicates
dim(duplicates)

SNPs_1 <- geno %>% distinct(Genotype, .keep_all = TRUE)
dim(SNPs_1)
SNPs_1[1:5, 1:5]

data <- read.csv('leg_coef.csv')
head(data)
dim(data)
length(unique(data$Genotype))

X0 <- SNPs_1
X0[1:5,9846:9848]
dim(X0)

X2 <- X0
dim(X2)

rownames(X2) <- X2[,1]
X2 <- X2[,-1]
X2[1:5,1:5]


NaN.freq.i <- 2864*0.2
prop.MAF.j <-  0.05

X2<- X2+1
SNPs <- X2
dim(SNPs)
SNPs[1:5, 1:5]


NaNs <- matrix(NA,nrow=1, ncol=dim(SNPs)[2])

for(i in 1:dim(SNPs)[2])
{
  NaNs[i] <- length(which(is.na(SNPs[,i])))
}

index.1 <- order(NaNs) 
NaNs2 <- NaNs[index.1]

NaN.freq <- unique(NaNs2)   
percentage.NaNs <- NaN.freq * 100 / dim(SNPs)[1]

write.table(NaN.freq, file='NaNs.freq.csv', sep=',',row.names=F, col.names=F)
write.table(percentage.NaNs, file='percentage.NaNs.csv', sep=',',row.names=F, col.names=F)

index.2 <- which(NaNs <= NaN.freq.i)
X <-SNPs[,index.2]
X[1:5, 1:5]

p <- colMeans(X,na.rm=T)/2    
p <- ifelse(p <= 0.5,p,1-p)

print( length( which( p >= prop.MAF.j ) ) )
index.3 <- which( p >= prop.MAF.j ) 
write.csv(X[,index.3], file='X.csv')
print( paste('NaNs_',NaN.freq.i,'_MAF_',prop.MAF.j,sep='') )

write.csv(X, file='X.csv')

dim(X)

X
sds <- apply(X,2,sd)
sds
sum(sds)

Y <- data
head(Y)
Y <- Y[,-1]
colIDy <- 6
IDs<-Y[,colIDy]

dim(X)

X3 <- scale(X)

sds <- apply(X3,2,sd)
td <- which(is.na(sds))
X4 <- X3[,-td]
dim(X4)

X5 <- X4 

G <- tcrossprod(X5)/dim(X5)[2]
 G[1:5, 1:5]
dim(G)
image(G)
diag(G)

order(diag(G))
rownames(G)<- SNPs_1$Genotype
#save(G,file='G.rda')

if(!is.null(colIDy)){ stopifnot(all(IDs%in%rownames(G))) }

if(!is.null(colIDy)){
  IDs<-factor(IDs,levels=rownames(G))
  Z<-as.matrix(model.matrix(~IDs-1))
  G<-tcrossprod(tcrossprod(Z,G),Z)
}

dim(G)
save(G,file='G.rda')



###### E #################



head(Y)
# We assign the column that links phenotypes with environments
colIDy <- 5
IDs <- Y[,colIDy]
Y[,colIDy]<-factor(Y[,colIDy])  

# Must create a matrix  
Z<-as.matrix(model.matrix(~Y[,colIDy]-1)) 
d<-colSums(Z)
V<-Z
n<-length(levels(Y[,colIDy]))
for(i in 1:ncol(Z)){ V[,i]<-V[,i]/sqrt(d[i]) } 
EVD <- list(vectors=V,values=d)
E <- tcrossprod(Z)
dim(E)
save(E,file='E.rda')
##save(EVD,file='EVD.E.rda')

image(E)


GxE <- G*E
EVD <- eigen(GxE)
image(GxE)
#order(diag(GxE))
#G[929,929]
save(GxE,file='GxE.rda')
dim(GxE)


