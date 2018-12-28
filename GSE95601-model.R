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
source("~/projects/bivalent/modalv/GSE75748_data_function.R")

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


set.seed(20)
register(SerialParam())


load("/mnt/gtklab01/ahjung/bivalent/data/GSE95601/fletcher.rda")
load("/mnt/gtklab01/ahjung/bivalent/data/GSE95601/clustered.rda")


### 123 housekeeping genes are found in the bimodal genes.
## data("housekeeping")
## hk = rownames(fletcher)[toupper(rownames(fletcher)) %in% housekeeping$V1]

## mfilt <- metric_sample_filter(counts(fletcher), 
##                               nreads = colData(fletcher)$NREADS,
##                               ralign = colData(fletcher)$RALIGN,
##                               pos_controls = rownames(fletcher) %in% hk,
##                               zcut = 3, mixture = FALSE,
##                               plot = TRUE)

## mfilt <- !apply(simplify2array(mfilt[!is.na(mfilt)]), 1, any)
## filtered <- fletcher[, mfilt]
#filtered <- filtered[!apply(counts(filtered),1,function(x) all(x==0)),]


## filtered <- makeFilterStats(filtered, filterStats="var", transFun = log1p)
## filtered <- filterData(filtered, percentile=1000, filterStats="var")
filtered <- fletcher[rownames(fletcher) %in% bimo,]
filtered

publishedClusters <- colData(filtered)[, "publishedClusters"]
col_clus <- c("transparent", "#1B9E77", "antiquewhite2", "cyan", "#E7298A", 
              "#A6CEE3", "#666666", "#E6AB02", "#FFED6F", "darkorchid2", 
              "#B3DE69", "#FF7F00", "#A6761D", "#1F78B4")
names(col_clus) <- sort(unique(publishedClusters))
table(publishedClusters)

clustered <- zinbwave(filtered, K = 50, X = "~ Batch", residuals = TRUE, normalizedValues = TRUE, BPPARAM=MulticoreParam(50))

assayNames(clustered)

###################### 
norm <- assay(clustered, "normalizedValues")
norm[1:3,1:3]

reducedDimNames(clustered)
#> [1] "zinbwave"
W <- reducedDim(clustered, "zinbwave")
dim(W)
#> [1] 747  50
W[1:3, 1:3]

W <- reducedDim(clustered)
d <- dist(W)
fit <- cmdscale(d, eig = TRUE, k = 2)
plot(fit$points, col = col_clus[as.character(publishedClusters)], main = "",
     pch = 20, xlab = "Component 1", ylab = "Component 2")
legend(x = "topleft", legend = unique(names(col_clus)), cex = .5, fill = unique(col_clus), title = "Sample")

plotClusters(clustered)

plotCoClustering(clustered)

pseudoCe <- clustered[,!primaryClusterNamed(clustered) %in% c("-1")]
X <- reducedDim(pseudoCe,type="zinbwave")
mds <- cmdscale(dist(X), eig = TRUE, k = 4)
lineages <- slingshot(mds$points, clusterLabels = primaryClusterNamed(pseudoCe), start.clus = "c1")

colorCl<-convertClusterLegend(pseudoCe,whichCluster="primary",output="matrixColors")[,1]
pairs(lineages, type="lineages", col = colorCl)

par(mfrow=c(5,10))
sapply(1:50, function(i) 
    plot(density(norm[bimodal[bimodal %in% rownames(clustered)],][i,])))
