suppressPackageStartupMessages({
  # Bioconductor
  library(BiocParallel)
  library(SingleCellExperiment)
  library(clusterExperiment)
  library(scone)
  library(zinbwave)
  library(slingshot)
  # CRAN
  library(gam)
  library(scater)
  library(doParallel)
  library(RColorBrewer)
})

load("/mnt/cbis/home/ahjung/projects/bivalent/modalv/data/GSE95601_oeHBCdiff_RSEM_eSet.Rda")
source("~/projects/bivalent/modalv/modalv_function.R")

sce <- SingleCellExperiment(assays=list(tpm=assayData(RSEM_eSet)$tpm_table))
filter <- rowSums(is.nan(assay(sce)))==ncol(assay(sce)) # remove NaN
sce <- sce[!filter,]
filter <- rowSums(assay(sce),na.rm=TRUE)==0 # remove cells with all zero
sce <- sce[!filter,]

messup <- function(sc_cell) {
    # sc_cell must be in log scale
    nnoise <- matrix(rnorm(nrow(sc_cell)*ncol(sc_cell),
                           mean=0, sd=0.1),ncol=ncol(sc_cell))
    sc_cell_messup <- sc_cell+nnoise
    return(sc_cell_messup)
}

logtpm <- log(assay(sce)+1)
logtpm_messup <- messup(logtpm)

logtpm_messup_bimodal <- apply(logtpm_messup, 1, test_bimodal)

cl <- makeCluster(detectCores()-2) # create a cluster with max-2 cores
registerDoParallel(cl) # register the cluster

res_tpm = foreach(i = 1:nrow(logtpm_messup),
    .combine = "rbind",
              .packages="mixtools")     %dopar% {
  # generate a bootstrap sample              
        fit_bimodal(logtpm_messup[i,],rownames(logtpm_messup)[i])
}


stopCluster(cl) # shut down the cluster

logtpm_messup_bimodal <- logtpm_messup_bimodal[!is.na(res_tpm$gene)]
res_tpm <- res_tpm[!is.na(res_tpm$gene),]

condition <- (logtpm_messup_bimodal == 0) &
    sapply(1:nrow(res_tpm), function(i) min(res_tpm[i,]$mu1, res_tpm[i,]$mu2) <= 1) &
        sapply(1:nrow(res_tpm), function(i) max(res_tpm[i,]$mu1, res_tpm[i,]$mu2) > 1)  &
            sapply(1:nrow(res_tpm), function(i) res_tpm[i,]$lambda1 <= 0.6)

bimo <- names(condition)[condition]
bimo <- sapply(bimo, function(x) strsplit(x,"[.]")[[1]][1])
set.seed(20)
register(SerialParam())

sapply(1:50, function(x) plottrough(logtpm_messup[bimo,][x,]))

findtrough <- function(densityy){
    trough <- which(diff(sign(diff(densityy)))==2)
    trough <- ifelse(length(trough)>1,trough[1],trough)
    return(trough)
}

findcutoff <- function(exp) {
    dens <- density(exp)
    return(dens$x[findtrough(dens$y)])
}

binarizeexp <- function(exp) {
    cutoff <- findcutoff(exp)
    return(exp > cutoff)
}

plottrough <- function(exp){
    dens <- density(exp)
    trough <- findtrough(dens$y)
    plot(dens,main=rownames(exp))
    abline(v=dens$x[trough],col="red")
}


logtpm_messup_bin <- t(apply(logtpm_messup[bimo,],1,binarizeexp))



tpmvalues <- assay(sce)

plot(log(apply(tpmvalues,1,mean),10),
     log(apply(tpmvalues,1,sd)^2/apply(tpmvalues,1,mean),10))

points(log(apply(tpmvalues[rownames(tpmvalues) %in% bimo,],1,mean),10),
     log(apply(tpmvalues[rownames(tpmvalues) %in% bimo,],1,sd)^2/apply(tpmvalues[rownames(tpmvalues) %in% bimo,],1,mean),10), col="red")

points(log(apply(tpmvalues[rownames(tpmvalues) %in% rownames(filtered),],1,mean),10),
     log(apply(tpmvalues[rownames(tpmvalues) %in% rownames(filtered),],1,sd)^2/apply(tpmvalues[rownames(tpmvalues) %in% rownames(filtered),],1,mean),10), col="blue")



load("/mnt/gtklab01/ahjung/bivalent/data/GSE95601/fletcher.rda")
load("/mnt/gtklab01/ahjung/bivalent/data/GSE95601/clustered.rda")


### 125 housekeeping genes are found in the bimodal genes.
data("housekeeping")
hk = rownames(fletcher)[toupper(rownames(fletcher)) %in% housekeeping$V1]

mfilt <- metric_sample_filter(counts(fletcher), 
                              nreads = colData(fletcher)$NREADS,
                              ralign = colData(fletcher)$RALIGN,
                              pos_controls = rownames(fletcher) %in% hk,
                              zcut = 3, mixture = FALSE,
                              plot = TRUE)
mfilt <- !apply(simplify2array(mfilt[!is.na(mfilt)]), 1, any)
filtered <- fletcher[, mfilt]
filtered <- filtered[!apply(counts(filtered),1,function(x) all(x==0)),]


filtered <- makeFilterStats(filtered, filterStats="var", transFun = log1p)
filtered <- filterData(filtered, percentile=1000, filterStats="var")
#filtered <- fletcher[rownames(fletcher) %in% bimo,]
filtered

publishedClusters <- colData(fletcher)[, "publishedClusters"]

col_clus <- c("black", "#1B9E77", "antiquewhite2", "cyan", "#E7298A", 
              "#A6CEE3", "#666666", "#E6AB02", "#FFED6F", "darkorchid2", 
              "#B3DE69", "#FF7F00", "#A6761D", "#1F78B4")
names(col_clus) <- sort(unique(publishedClusters))
table(publishedClusters)

clustered_bimo <- zinbwave(filtered, K = 50, X = "~ Batch", residuals = TRUE, normalizedValues = TRUE, BPPARAM=MulticoreParam(50))

assayNames(clustered)

###################### 
norm <- assay(clustered, "normalizedValues")
norm[1:3,1:3]

norm <- logtpm[bimo,]
norm <- logtpm_messup_bin

d <- as.dist(1 - cor(norm, method = 'sp') ^ 2)
PC <- prcomp(d)
rd1 <- PC$x[,1:2]
plot(rd1, pch=16, asp = 1)

library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1,10)$classification
points(rd1, col = col_clus[cl1], pch=16, asp = 1)

hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(20)


heatmap.2(matrix(as.numeric(norm),nrow=972), #rownames(sc_cell_messup) %in% bigenes
          trace="none",
          ColSideColors=col_clus[cl1],
#          distfun = function(x) as.dist(1 - cor(t(x), method = 'sp') ^ 2),
          col=hmcol#,
#                                        Colv=F
          )

???????????????????????????



library(ade4)

my.mst <- mstree(dist(rd1),1)
s.label(rd1, clab = 0, cpoi = 2, neig = my.mst, cnei = 1, label = cl1)




M <- mst(d)
#opar <- par()
par(mfcol = c(2, 2))
plot(M)
plot(M, graph = "nsca")
E(mst)$width  <-  3
plot(M, x1 = PC$x[, 1], x2 = PC$x[, 2])#,col=brewer.pal(9,"Set1")[cl1])
par(opar)


library(igraph)
 
mst <- minimum.spanning.tree(g)
E(mst)$color <- 'red'
 
layout  <- layout.kamada.kawai(mst)
 
plot(g, layout=layout)
par(new=T)
E(mst)$width  <-  3
plot(mst, layout=layout)
par(new=F)


d <- dist(rd1)
G <- graph.adjacency(d)
mmst <- minimum.spanning.tree(G)
plot(mmst, main = "MST")
plot(G)

plot(G, edge.color=c("red","green")[sign(E(G)$weight)/2 + 1.5], 
     edge.width = 3 *abs(E(G)$weight))



g <- erdos.renyi.game(10, 3/10)
mst <- minimum.spanning.tree(g)



cl2 <- kmeans(rd1, centers = 4)$cluster
colData(clustered)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

sce_sling <- slingshot(dist(norm))

summary(sce$slingPseudotime_1)


pseudoCe <- clustered[,!primaryClusterNamed(clustered) %in% c("-1")]
X <- reducedDim(pseudoCe,type="zinbwave")
mds <- cmdscale(dist(X), eig = TRUE, k = 4)
lineages <- slingshot(mds$points, clusterLabels = primaryClusterNamed(pseudoCe), start.clus = "c1")

colorCl<-convertClusterLegend(pseudoCe,whichCluster="primary",output="matrixColors")[,1]
pairs(lineages, type="lineages", col = colorCl)

par(mfrow=c(5,10))
sapply(1:50, function(i) 
    plot(density(norm[bimodal[bimodal %in% rownames(clustered)],][i,])))
