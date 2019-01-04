setwd("~/projects/bivalent/modalv")

source("GSE75748_data.R")
source("GSE75748_functions.R")

## sc_cell_original <- sc_cell
sc_cell <- sc_cell[,sc_cell_coldata$cell != c("H9")]
sc_cell <- sc_cell[!apply(sc_cell,1,sum)==0,]
sc_cell_coldata_H1 <- sc_cell_coldata[sc_cell_coldata$cell != c("H9"),]
sc_cell_coldata_H1$cell <- as.factor(as.character(sc_cell_coldata_H1$cell))

sc_time_tpm <- get_tpm(sc_time)
log_sc_time_tpm <- log(sc_time_tpm+1)
log_sc_time_tpm_messup <- messup(log_sc_time_tpm)

log_sc_time_tpm_messup_bimodal <- apply(log_sc_time_tpm_messup, 1, test_bimodal)

bimofit <- fit_bimodal_multi(log_sc_time_tpm_messup)

bimocondition <- filter_condition(log_sc_time_tpm_messup_bimodal, bimofit)

nclass <- table(sc_time_coldata$exp)
cols <- get_colors(nclass)

## par(mfrow=c(5,10))
## color_check(cols,sc_time_coldata$exp)
## color_check(cols,cl1)
## sapply(1:49, function(x)
##     plot_cluster(nclass, log_sc_time_tpm[bimocondition,][x,], rownames(log_sc_time_tpm)[bimocondition][x],
##                  cols))


## TSNE
set.seed(2)

bimocondition_sub <- (rownames(log_sc_time_tpm) %in% bigenes) & bimocondition

tsne_out <- Rtsne(t(unique(log_sc_time_tpm[bimocondition_sub,])),
                  pca=FALSE,
                  perplexity=30,
                  theta=0.0)

ggplot(data=data.frame("C1"=tsne_out$Y[,1],"C2"=tsne_out$Y[,2],
           "type"=sc_time_coldata$exp),
       aes(C1,C2,col=type)) +
    geom_point()

rd1 <- tsne_out$Y[,1:2]
rownames(rd1) <- sc_time_coldata$exp

## pseudotime

log_sc_time_tpm_messup_bin <- t(apply(log_sc_time_tpm_messup[bimocondition_sub,],1,binarizeexp))

norm <- log_sc_time_tpm_messup_bin

d <- as.dist(1 - cor(norm, method = 'sp') ^ 2)
PC <- prcomp(d)
rd1 <- PC$x[,1:2]
plot(rd1, pch=16, asp = 1,col=cols[sc_time_coldata$exp])
plot(rd1, pch=16, asp = 1,col=cols[cl1])

scatterplot3d(PC$x[,1],PC$x[,2],PC$x[,3],color=cols[sc_time_coldata$exp])

cl1 <- Mclust(rd1,6)$classification


# MST
my.mst <- mstree(dist(rd1),1)

s.label(rd1, clab = 0, cpoi = 2, neig = my.mst, cnei = 1, label = cols[sc_time_coldata_H1$cell])

points(rd1, col = cols[sc_time_coldata$exp], pch=16, asp = 1)
points(rd1, col = cols[cl1], pch=16, asp = 1)

# finding shortest path to order cells

mygraph <- graph_from_adjacency_matrix(neig2mat(my.mst))
a <- shortest.paths(mygraph, 1) # starting from any cell from time 0

# heatmap of cells ordered according to shortest path (pseudotime)

hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(20)
hmcol <- colorRampPalette(brewer.pal(9, "YlGnBu"))(20)
heatmap.2(matrix(as.numeric(log_sc_time_tpm[bimocondition_sub,]),
                 nrow=nrow(norm))[,order(as.numeric(a))], #rownames(sc_time_messup) %in% bigenes
          trace="none",
          ColSideColors=cols[cl1][order(as.numeric(a))],
          distfun = function(x) as.dist(1 - cor(t(x), method = 'sp') ^ 2),
          col=hmcol,
          Colv=F
          )

norm_bin <- (norm*1)[,order(as.numeric(a))]

ordered_cell <- matrix(as.numeric(log_sc_time_tpm[bimocondition_sub,]),
                       nrow=nrow(norm))[,order(as.numeric(a))]
rownames(ordered_cell) <- rownames(log_sc_time_tpm[bimocondition_sub,])
nclass <- table(cl1[order(as.numeric(a))])

ordered_cell_ratio <- t(sapply(1:nrow(ordered_cell),
                               function(i)
                                   clusterratio(ordered_cell[i,],nclass)[,"ON"]))

ordered_cell_ratio_wd <- t(sapply(1:nrow(ordered_cell),
                                  function(i)
                                   windowratio(ordered_cell[i,],20)[,"ON"]))

## hr <- hclust(as.dist(1-cor(t(ordered_cell_ratio), method="pearson")),
##              method="complete")

hr <- hclust(dist(ordered_cell_ratio_wd))
## hr <- hclust(dist(muin_bin, method="euclidean"),method="ward.D2")

heatmap.2(ordered_cell[hr$order,][names(cl_ratio)[cl_ratio==6],],
          trace="none",
#ColSideColors=cols[as.numeric(sc_time_coldata$exp)][order(as.numeric(a))],
#          distfun = function(x) as.dist(1 - cor(t(x), method = 'sp') ^ 2),
          col=hmcol,
          ## dist=dist,
          ## hclust=hclust,
          ColSideColors=cols[cl1][order(as.numeric(a))],
          Colv=FALSE,
          ## Rowv=FALSE,
          dendrogram = "none"
          )



## heatmap.2(muin_bin,
##           trace="none", col=hmcol,
##           dendrogram="none")


## idmat <- expand_bit(nclass)
## idmat_c <- t(get_bit(1:6))

## muin_bin <- muin_all(idmat, norm_bin)
## muin_bin <- muin_all(idmat_c, ordered_cell_ratio)

tsne_out <- Rtsne(ordered_cell_ratio_wd,
                  pca=FALSE,
                  perplexity=30,
                  theta=0.0)

rd1 <- tsne_out$Y[,1:2]
rownames(rd1) <- rownames(ordered_cell)
cl_ratio <- Mclust(rd1,6)$classification

heatmap.2(ordered_cell[order(cl_ratio),],
          trace="none",
#ColSideColors=cols[as.numeric(sc_time_coldata$exp)][order(as.numeric(a))],
#          distfun = function(x) as.dist(1 - cor(t(x), method = 'sp') ^ 2),
          col=hmcol,
          ## dist=dist,
          ## hclust=hclust,
          ColSideColors=cols[cl1][order(as.numeric(a))],
          Colv=FALSE,
          Rowv=FALSE,
          dendrogram = "none"
          )

par(mfrow=c(4,2))

plot_cluster(table(cl1), log_sc_time_tpm["SOX2",],"SOX2",cols,ylim=c(0,2),xlim=c(-2,8))

plot_cluster(table(cl1), log_sc_time_tpm["CER1",],"CER1",cols,ylim=c(0,1),xlim=c(-2,12))

plot_cluster(table(cl1), log_sc_time_tpm["DNMT3B",],"DNMT3B",cols,ylim=c(0,1),xlim=c(-2,8))

plot_cluster(table(cl1), log_sc_time_tpm["COCH",],"COCH",cols,ylim=c(0,0.5),xlim=c(-2,8))

plot_cluster(table(cl1), log_sc_time_tpm["BAG3",],"BAG3",cols,ylim=c(0,0.8),xlim=c(-2,8))

plot_cluster(table(cl1), log_sc_time_tpm["PRDM1",],"PRDM1",cols,ylim=c(0,0.8),xlim=c(-2,8))

plot_cluster(table(cl1), log_sc_time_tpm["CYP26A1",],"CYP26A1",cols,ylim=c(0,0.8),xlim=c(-2,11))



par(mfrow=c(1,6))
#color_check(cols,sc_time_coldata$exp)
color_check(cols,1:6)

sapply(100:149, function(x)
    plot_cluster(table(cl1), log_sc_time_tpm[bimocondition_sub,][x,], rownames(log_sc_time_tpm)[bimocondition_sub][x],
                 cols))


#rownames(ordered_cell_ratio_wd) <- rownames(ordered_cell)

sapply(c("BAG3","BID","GATA6","NBPF15","RARRES2","MRPS26"),
       function(x) scatter.smooth(ordered_cell_ratio_wd[x,], lpars = list(col="red", lwd=3), main=x, ylab="Fraction of ON", xlab="windows"))

#plot.new()

cumsum(table(cl1))

invisible(sapply(c("BAG3","BID","GATA6","NBPF15","RARRES2","MRPS26"),
       function(x) plot(density(as.numeric(log_sc_time_tpm[x,712:758])),main="")))

table(cl1),log_sc_time_tpm[x,],
                                x,
                                cols))
color_check(cols,1:6)

library(gridExtra)
p6 <- plot_time(nclass,
               as.numeric(log_sc_time_tpm["MRPS26",]),
               "MRPS26")
grid.arrange(p,p2,p3,p4,p5,p6, ncol=6)



#####TODO ggplot, facet_wrap, plot_time function
