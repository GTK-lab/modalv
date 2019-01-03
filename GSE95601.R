df_protocol <- read.table("/mnt/gtklab01/ahjung/bivalent/data/GSE95601/GSE95601_oeHBCdiff_RSEM_eSet_protocolData.txt")

df <- read.table("/mnt/gtklab01/ahjung/bivalent/data/GSE95601/GSE95601_oeHBCdiff_RSEM_eSet_expectedCounts_table.txt")

df2 <- df[!apply(df,1,sum,na.rm=TRUE) < 10,]

test <- t2gm[t2gm$ens_gene %in% rownames(df2),]
test <- unique(test[,c("ens_gene","length")])
test <- test[!duplicated(test$ens_gene),]

df2_test <- df2[unique(test$ens_gene),]

df2_rpk <- invisible(t(sapply(1:nrow(df2_test),
                                 function(x) df2_test[x,]/(test$length[x]/1000))))
rownames(df2_rpk) <- rownames(df2_test)

df2_rpk2 <- apply(df2_rpk,2,function(x) as.numeric(x))
rownames(df2_rpk2) <- rownames(df2_rpk)

logdf2_rpk <- log(df2_rpk2+1)

df2_tpmfactor <- apply(df2_rpk2,2,function(x) sum(as.numeric(x)))

df2_tpm <- invisible(sapply(1:ncol(df2_rpk2),function(x) df2_rpk2[,x]/df2_tpmfactor[x]))
colnames(df2_tpm) <- colnames(df2_rpk)
rownames(df2_tpm) <- rownames(df2_rpk)


df2_tpm <- df2_tpm[!apply(df2_tpm,1,sum)==0,]
df2_tpm <- df2_tpm*1000000
logdf2_tpm <- log(df2_tpm+1)

logdf2_tpm_messup <- messup(logdf2_tpm)

logdf2_tpm_messup_bimodal <- apply(logdf2_tpm_messup, 1, test_bimodal)

cl <- makeCluster(detectCores()-2) # create a cluster with max-2 cores
registerDoParallel(cl) # register the cluster

res_tpm = foreach(i = 1:nrow(logdf2_tpm_messup),
    .combine = "rbind",
              .packages="mixtools")     %dopar% {
  # generate a bootstrap sample              
        fit_bimodal(logdf2_tpm_messup[i,],rownames(logdf2_tpm_messup)[i])
}


stopCluster(cl) # shut down the cluster

res_tpm <- res_tpm[!is.na(res_tpm$gene),]

condition <- (logdf2_tpm_messup_bimodal == 0) &
    (min(res_tpm$mu1, res_tpm$mu2) < 1)  &
         (res_tpm$lambda1 <= 0.8)

par(mfrow=c(5,10))

sapply(1:50, function(x)
    #plot(density(logdf2_tpm[condition,][x,]),main=rownames(logdf2_tpm)[condition][x]))

plot(density(logdf2_tpm[condition,][x,]),main=rownames(logdf2_tpm)[condition][x], pch=20
))

set.seed(1)
tsne_out <- Rtsne(t(unique(logdf2_tpm[condition,])),pca=FALSE,perplexity=30,theta=0.0)

## diffexp <- c("Krt5","Aqp3","Krt14","Krt17","Cbr2","Rtn1","Car12","Gstm2","Apoe","Tubb3","Trp63","Gstm1","Scgb1c1","Hes1","Cyr61","Zfp36","Fstl")





ggplot(data=data.frame("C1"=tsne_out$Y[,1],"C2"=tsne_out$Y[,2],"lib"=d_tsne_1$cl_hierarchical), aes(C1,C2,col=lib)) +
    geom_point()


hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(50)
cols_all <- palette(brewer.pal(8, "Dark2"))[as.numeric(mcluster)]

heatmap.2(logdf2_tpm[condition,][,order(mcluster)], #  & rownames(sc_cell_messup) %in% bigenes
trace="none",ColSideColors=cols_all,
#RowSideColors=cols2,
col=hmcol
#                                        Colv=F
          )


pca <- prcomp(t(logdf2_tpm[condition,]), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, pch=16, asp = 1)

library(destiny)

dm <- DiffusionMap(t(logdf2_tpm[condition,]))
rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)

plot(rd2, pch=16, asp = 1)

reducedDims(sim) <- SimpleList(PCA = rd1, DiffMap = rd2)

library(mclust)

cl1 <- Mclust(rd1)$classification
#colData(sim)$GMM <- cl1

plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 4)$cluster
#colData(sim)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

sce <- slingshot(rd1, clusterLabels = 'GMM', reducedDim = 'PCA')

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(rd1, col = colors[cut(sce$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)

d_tsne_1 <- as.data.frame(tsne_out$Y)
fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))

d_tsne_1$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=10))  

sub_tpm <- logdf2_tpm[condition,]
mcluster <- d_tsne_1$cl_hierarchical

                                        # per cluster, get ratio between two? (zero and non zero)
ncluster <- 1

cluster_ratio <- function(ncluster){
x <- (apply(logdf2_tpm[condition,mcluster==ncluster],1,function(x) sum(x>=1)))/ncol(logdf2_tpm[condition,mcluster==ncluster])*100
return(x) }

cratio <- data.frame("cluster1"=cluster_ratio(1))
invisible(sapply(2:10, function(x) cratio <<- cbind(cratio, cluster_ratio(x))))

colnames(cratio) <- paste0("cluster",1:10)


heatmap.2(as.matrix(cratio),
trace="none",#ColSideColors=cols_all,
#RowSideColors=cols2,
col=hmcol,distfun = function(x) as.dist(1 - cor(t(x), method = 'sp') ^ 2)
          )

for (i in sample(1:nrow(sub_tpm),50)) {
plot(density(sub_tpm[i,mcluster==8]),ylim=c(0,0.8),main=rownames(sub_tpm)[i],col="red",xlim=c(-2,10))
lines(density(sub_tpm[i,mcluster==9]),col="orange")
lines(density(sub_tpm[i,mcluster==7]),col="green")
## lines(density(sub_tpm[i,mcluster==5]),col="blue")
## lines(density(sub_tpm[i,mcluster==10]),col="purple")
## lines(density(sub_tpm[i,mcluster==1]),col="black")
}

