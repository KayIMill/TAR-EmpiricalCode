
library(vegan)
library(gmt)
setwd(choose.dir())

## Bringing in matrix from mothur ##
## Code from Mario Muscarella ##

source("DiversityFunctions.r")  

shared = "AGTAR.sil.shared"
design = "AGTAR.design"

# Import Site by OTU Matrix
AGTAR <- t(read.otu(shared, "0.03"))
design <- read.delim(design, header=T, row.names=1)
  
# Remove Zero Sum OTUs
AGTAR <- AGTAR[,!(colSums(abs(AGTAR)) ==0)] 

# Select only OTUs with > 1 observations in >1 sample
AGTAR <- AGTAR[rowSums(AGTAR) > 1,]
AGTAR <- AGTAR[rowSums((AGTAR > 0)*1) > 1, ]
head(AGTAR)

# Calculate Presence Absence
dataPA <- (AGTAR > 0)*1 

tdataPA <- decostand(t(dataPA),method="log")
head(tdataPA)

## Geographic Distance ##
## From Zinger## 

env2 <- read.table("coords_ag.txt", header=TRUE) 
head(env2)
env2
rownames = list(env2$cid)
rownames

dist.geo2 = matrix(0, nrow(env2), nrow(env2), dimnames=list(rownames(rownames),rownames(rownames))
for(j in 1:(nrow(env2)-1)) for(i in (j+1):nrow(env2)){ #for each sample pairwise comparison:
  dist.geo2[i,j] = 100 * (geodist(as.numeric(as.vector(env2[i,"LAT"])), as.numeric(as.vector(env2[i,"LONG"])), as.numeric(as.vector(env2[j,"LAT"])), as.numeric(as.vector(env2[j,"LONG"])), units="km"))
  }

dist.geo2


## distances for 33 samples aka only DNA or RNA ##

halfcoords <- read.table("coords_half.txt", header=TRUE) 
head(halfcoords)
halfcoords
rownames = list(halfcoords$cid)
rownames

dists = matrix(0, nrow(halfcoords), nrow(halfcoords), dimnames=list(rownames(rownames),rownames(rownames)))

for(j in 1:(nrow(halfcoords)-1)) for(i in (j+1):nrow(halfcoords)){ #for each sample pairwise comparison:
  dists[i,j] = 100 * (geodist(as.numeric(as.vector(halfcoords[i,"LAT"])), as.numeric(as.vector(halfcoords[i,"LONG"])), as.numeric(as.vector(halfcoords[j,"LAT"])), as.numeric(as.vector(halfcoords[j,"LONG"])), units="km"))
  }

dists

write.csv(dists, "dists33.csv")
write.csv(dist.geo2, "distances.csv")


## subsetting DNA/RNA ## 
# datafram is transformed presence absence ##
# tdataPA ##

tdataPA <- decostand(t(dataPA),method="log")
head(tdataPA)

DNA.PA <- tdataPA[c(1:11, 23:33, 45:55),c(1:11, 23:33, 45:55)]
DNA.PA


RNA.PA <- tdataPA[c(12:22, 34:44, 56:66),c(12:22, 34:44, 56:66)]
RNA.PA

## Community Distance ##

##Subsets 

DNAPA.BC.dist <- as.matrix(vegdist(DNA.PA, method="bray"))
DNAPA.BC.dist
write.csv(DNAPA.BC.dist, "DNABCdist.csv")

DNAPA.J.dist <- as.matrix(vegdist(DNA.PA, method="jaccard"))
DNAPA.J.dist

RNAPA.BC.dist <- as.matrix(vegdist(RNA.PA, method="bray"))
RNAPA.BC.dist
write.csv(RNAPA.BC.dist, "RNABCdist.csv")


## All data

bray.PA.dist <- vegdist(decostand(t(dataPA),method="log"),method="bray")
bray.PA.dist

bray.PA.matrix = as.matrix(vegdist(decostand(t(dataPA),method="log"),method="bray"))
bray.PA.matrix


#exporting data
write.csv(bray.PA.matrix, "braydistances.csv")
write.csv(jaccards.PA.matrix, "jaccarddistances.csv")



## Trying to read Matricies ##
## Code from Zimmer ##

eco = list(ELAG=which(env$cid==2:12),
             ELFOR=which(env$cid==17:27),
             KBSAG=which(env$cid==32:42),
             KBSFOR=which(env$cid==47:57),
             RUAG=which(env$cid==62:72),
             RUFOR=which(env$cid==77:87))

geo.sor.raw = list()
for (i in 1:length(eco)){
	geo.sor.raw[[i]] = lm(log(as.vector(as.dist(MAT.sor.raw[eco[[i]],eco[[i]]], diag=F, upper=F))+0.01)~log(as.vector(as.dist(dist.geo[eco[[i]],eco[[i]]], diag=F, upper=F))+0.01))
    print(i)
}


## Trying to read Matricies ## 
## The for-loops ##
## matrix_name[row#, col#]
## geograpic distances in dists

dists.vector.test <- NULL
for (i in 1:33) dists.vector.test <- append(dists.vector.test, dists[y:33, x])
dists.vector.test


## This just needs to iterate through 33 times... ## 

x = 1
y = 2

dists.vector2 <- append(dists.vector, dists[y:33, x])
dists.vector2
x
y

x = x + 1
y = y + 1



## Graphing Distance Decay Relationships ## 

dds <- read.table("finaldataset.txt", header= TRUE)
head(dds)

lmDNA <- lm(dds$DNAdisim~dds$lndist)
lmRNA <- lm(dds$RNAdisim~dds$lndist)


plot(dds$lndist,dds$DNAdisim, xlab="Geographic Distance (ln(m))", ylab="Bray-Curtis dissimilarity", pch=16, col="blue")
points(dds$lndist, dds$RNAdisim, pch=18, col="red")
abline(lmRNA, col="red")
abline(lmDNA, col="blue")

title(main="Distance Decay")

lmRNA
lmDNA