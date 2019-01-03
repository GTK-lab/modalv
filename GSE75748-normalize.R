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
heatmap.2(matrix(as.numeric(log_sc_time_tpm[bimocondition_sub,]),
                 nrow=nrow(norm))[,order(as.numeric(a))], #rownames(sc_time_messup) %in% bigenes
          trace="none",
          ColSideColors=cols[sc_time_coldata$exp][order(as.numeric(a))],
          distfun = function(x) as.dist(1 - cor(t(x), method = 'sp') ^ 2),
          col=hmcol,
          Colv=F
          )
