setwd("~/projects/bivalent/modalv")

source("GSE75748_data.R")
source("GSE75748_function.R")

## sc_cell_original <- sc_cell
## sc_cell <- sc_cell[,sc_cell_coldata$cell != c("H9")]
## sc_cell <- sc_cell[!apply(sc_cell,1,sum)==0,]
## sc_cell_coldata_H1 <- sc_cell_coldata[sc_cell_coldata$cell != c("H9"),]
## sc_cell_coldata_H1$cell <- as.factor(as.character(sc_cell_coldata_H1$cell))
sc_cell_coldata_H1$cell = factor(sc_cell_coldata_H1$cell,
    levels(sc_cell_coldata_H1$cell)[c(3,1,2,4,5,6)])


sc_cell_tpm <- get_tpm(sc_cell)
log_sc_cell_tpm <- log(sc_cell_tpm+1)
log_sc_cell_tpm_messup <- messup(log_sc_cell_tpm, 0.1)

log_sc_cell_tpm_messup_bimodal <- apply(log_sc_cell_tpm_messup, 1, test_bimodal)

bimofit <- fit_bimodal_multi(log_sc_time_tpm_messup)

bimocondition <- filter_condition(log_sc_time_tpm_messup_bimodal, bimofit)

mclass <- sc_time_coldata$exp
nclass <- table(mclass)
cols <- get_colors(nclass)
## par(mfrow=c(5,10))
## color_check(cols,sc_cell_coldata$exp)
## color_check(cols,cl1)
## sapply(1:49, function(x)
##     plot_cluster(nclass, log_sc_cell_tpm[bimocondition,][x,], rownames(log_sc_cell_tpm)[bimocondition][x],
##                  cols))


## TSNE
set.seed(1)

bimocondition_sub <- (rownames(log_sc_cell_tpm) %in% bigenes) & bimocondition

tsne_out <- Rtsne(t(unique(log_sc_cell_tpm_messup[bimocondition_sub,])),
                  pca=FALSE,
                  perplexity=30,
                  theta=0.0)

tsne1 <- ggplot(data=data.frame("C1"=tsne_out$Y[,1],"C2"=tsne_out$Y[,2],
           "type"=mclass),
       aes(C1,C2,col=type)) +
    geom_point()

rd1 <- tsne_out$Y[,1:2]
rownames(rd1) <- mclass

cl1 <- Mclust(rd1,length(nclass))$classification
#plot(rd1_pc, pch=16, asp = 1,col=cols[cl1])

tsne2 <- ggplot(data=data.frame("C1"=tsne_out$Y[,1],"C2"=tsne_out$Y[,2],
           "type"=as.factor(cl1)),
       aes(C1,C2,col=type)) +
    geom_point()

grid.arrange(tsne1,tsne2,ncol=2)

## pseudocell

# binarize 
log_sc_cell_tpm_messup_bin <- t(apply(log_sc_cell_tpm_messup[bimocondition_sub,],1,binarizeexp))

norm <- log_sc_cell_tpm_messup_bin#[,mclass %in% c("H1","DEC","EC")]

## d <- as.dist(1 - cor(norm, method = 'sp') ^ 2)
## PC <- prcomp(d)
## rd1_pc <- PC$x[,1:2]

## par(mfrow=c(1,2))
## plot(rd1_pc, pch=16, asp = 1,col=cols[mclass])

#scatterplot3d(PC$x[,1],PC$x[,2],PC$x[,3],color=cols[sc_cell_coldata$exp])


# MST
my.mst <- mstree(dist(rd1),1)
mygraph <- graph_from_adjacency_matrix(neig2mat(my.mst))
a <- shortest.paths(mygraph, 1) # starting from any cell from cell 0
a2 <- get.shortest.paths(mygraph,1) # need to keep the option to specify end point

## using igraph
longest <- which(unlist(lapply(a2$vpath,length))==max(unlist(lapply(a2$vpath,length))))

V(mygraph)$color <- cols[cl1]
E(mygraph, path=a2$vpath[[longest]])$color <- "red"
V(mygraph)[as.numeric(a2$vpath[[longest]])]$color <- "black"

plot(mygraph,vertex.size=2, layout=rd1,vertex.label=NA,edge.arrow.size=0.5)

## s.label(rd1_pc, clab = 0, cpoi = 1, neig = my.mst, cnei = 1)#, label = cols[sc_cell_coldata_H1$cell])

## points(rd1_pc, col = cols[mclass], pch=16, asp = 1, cex=0.5)
## points(rd1_pc, col = cols[cl1], pch=16, asp = 1, cex=0.5)

# finding shortest path to order cells




p1 <- ggplot(data.frame("x"=1:ncol(norm),"y"=rep(1,ncol(norm)),"group"=mclass[order(as.numeric(a))]),
            aes(x=x, y=y, fill=group)) +
  geom_bar(stat="identity")+theme_minimal() + scale_fill_manual(values=cols)

p2 <- ggplot(data.frame("x"=1:ncol(norm),"y"=rep(1,ncol(norm)),"group"=as.factor(cl1[order(as.numeric(a))])),
            aes(x=x, y=y, fill=group)) +
  geom_bar(stat="identity")+theme_minimal() + scale_fill_manual(values=cols)

grid.arrange(p1,p2,ncol=1)

# heatmap of cells ordered according to shortest path (pseudocell)

hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(20)
#hmcol <- colorRampPalette(brewer.pal(9, "YlGnBu"))(20)
heatmap.2(matrix(as.numeric(log_sc_cell_tpm[bimocondition_sub,]),
                 nrow=nrow(norm))[,mclass %in% c("H1","DEC","EC")][,order(as.numeric(a))], #rownames(sc_cell_messup) %in% bigenes
          trace="none",
          ColSideColors=cols[mclass][mclass %in% c("H1","DEC","EC")][order(as.numeric(a))],
          distfun = function(x) as.dist(1 - cor(t(x), method = 'sp') ^ 2),
          col=hmcol,
#          Colv=F,
#          Rowv=F,
#          dendrogram = "none"
          )

norder <- order(as.numeric(cl1))  # order(as.numeric(a))

norm_bin <- (norm*1)[,norder]
ordered_cell <- matrix(as.numeric(log_sc_cell_tpm[bimocondition_sub,]),
                       nrow=nrow(norm))[,norder]
rownames(ordered_cell) <- rownames(log_sc_cell_tpm[bimocondition_sub,])
nclass <- table(cl1[norder])

## ordered_cell_ratio <- t(sapply(1:nrow(ordered_cell),
##                                function(i)
##                                    clusterratio(ordered_cell[i,],nclass)[,"ON"]))

ncluster <- 30
ordered_cell_ratio_wd <- t(sapply(1:nrow(ordered_cell),
                                  function(i)
                                   windowratio(ordered_cell[i,],ncluster)[,"ON"]))

## hr <- hclust(as.dist(1-cor(t(ordered_cell_ratio), method="pearson")),
##              method="complete")

hr <- hclust(dist(ordered_cell_ratio_wd))
## hr <- hclust(dist(muin_bin, method="euclidean"),method="ward.D2")

heatmap.2(ordered_cell_ratio_wd[subcondition,],#hr$order,],
          trace="none",
#ColSideColors=cols[as.numeric(sc_cell_coldata$exp)][order(as.numeric(a))],
#          distfun = function(x) as.dist(1 - cor(t(x), method = 'sp') ^ 2),
          col=hmcol,
          ## dist=dist,
          ## hclust=hclust,
#          ColSideColors=cols[cl1][order(as.numeric(a))],
          Colv=FALSE,
          Rowv=FALSE,
          dendrogram = "none",
          RowSideColors=as.character(cl_ratio[order(cl_ratio)]),
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

rd1_wd <- tsne_out$Y[,1:2]
rownames(rd1_wd) <- rownames(ordered_cell)

cl_ratio <- Mclust(rd1_wd,8)$classification

heatmap.2(ordered_cell[order(cl_ratio),],
          trace="none",
#ColSideColors=cols[as.numeric(sc_cell_coldata$exp)][order(as.numeric(a))],
#          distfun = function(x) as.dist(1 - cor(t(x), method = 'sp') ^ 2),
          col=hmcol,
          ## dist=dist,
          ## hclust=hclust,
          ColSideColors=cols[cl1][norder],
          Colv=FALSE,
          Rowv=FALSE,
          RowSideColors=as.character(cl_ratio[order(cl_ratio)]),
          dendrogram = "none"
          )



plot_cluster(table(cl1), log_sc_cell_tpm["SOX2",],"SOX2",cols,ylim=c(0,2),xlim=c(-2,8))

plot_cluster(table(cl1), log_sc_cell_tpm["CER1",],"CER1",cols,ylim=c(0,1),xlim=c(-2,12))

plot_cluster(table(cl1), log_sc_cell_tpm["DNMT3B",],"DNMT3B",cols,ylim=c(0,1),xlim=c(-2,8))

plot_cluster(table(cl1), log_sc_cell_tpm["COCH",],"COCH",cols,ylim=c(0,0.5),xlim=c(-2,8))

plot_cluster(table(cl1), log_sc_cell_tpm["BAG3",],"BAG3",cols,ylim=c(0,0.8),xlim=c(-2,8))

plot_cluster(table(cl1), log_sc_cell_tpm["PRDM1",],"PRDM1",cols,ylim=c(0,0.8),xlim=c(-2,8))

plot_cluster(table(cl1), log_sc_cell_tpm["CYP26A1",],"CYP26A1",cols,ylim=c(0,0.8),xlim=c(-2,11))



par(mfrow=c(1,6))
#color_check(cols,sc_cell_coldata$exp)
color_check(cols,1:6)

sapply(100:149, function(x)
    plot_cluster(table(cl1), log_sc_cell_tpm[bimocondition_sub,][x,], rownames(log_sc_cell_tpm)[bimocondition_sub][x],
                 cols))


#rownames(ordered_cell_ratio_wd) <- rownames(ordered_cell)
par(mfrow=c(1,6))

sapply(c("BAG3","GATA6","NBPF15","BID","RARRES2","MRPS26"),
       function(x) scatter.smooth(ordered_cell_ratio_wd[x,],
                                  lpars = list(col="red", lwd=3),
                                  main=x,
                                  ylab="Fraction of ON",
                                  xlab="windows"))

#plot.new()

cumsum(table(cl1))

invisible(sapply(c("BAG3","BID","GATA6","NBPF15","RARRES2","MRPS26"),
       function(x) plot(density(as.numeric(log_sc_cell_tpm[x,712:758])),main="")))

table(cl1),log_sc_cell_tpm[x,],
                                x,
                                cols))
color_check(cols,1:6)

library(gridExtra)

genes <- c("BAG3","GATA6","NBPF15","BID","RARRES2","MRPS26")

gene <- genes[6]
p6 <- plot_cell(nclass,
               as.numeric(log_sc_cell_tpm[gene,]),
               gene)

grid.arrange(p,p2,p3,p4,p5,p6, ncol=6)

library(metrumrg)

get_cv <- function(expmat, logstatus) {
    m_mean <- apply(expmat,1,mean)
m_sd <- apply(expmat,1,sd)

if (logstatus) {
m_cv <- sapply(1:nrow(expmat),
               function(i) cvLognormal(m_sd[i], m_mean[i])) } else {
m_cv <- sapply(1:nrow(expmat),
      function(i) cvNormal(m_sd[i], m_mean[i])) }

    return(m_cv)
}

plot_cv <- function(expmat,logstatus, plotstatus) {
                                        # logstatus shows whether the data is in log-scale or not
    # plotstatus shows whether to call for plot or points
    m_mean <- apply(expmat,1,mean)
m_cv <- get_cv(expmat,logstatus)
if (plotstatus) {
plot(m_mean, m_cv) } else {
points(m_mean, m_cv, col="red") }

}

# comparing variance
par(mfrow=c(1,2))
plot_cv(log_sc_time_tpm, TRUE, TRUE)
plot_cv(ordered_cell, TRUE, FALSE)

plot_cv(log_sc_time_tpm[vargenes,],TRUE,FALSE)

m_cv <- get_cv(log_sc_time_tpm,TRUE)
vargenes <- names(m_cv)[order(m_cv,decreasing=TRUE)][1:length(bimogenes)]



## boxplot(apply(log_sc_cell_tpm,1,mean))
## boxplot(apply(ordered_cell,1,mean))
## boxplot(apply(log_sc_cell_tpm[order(apply(log_sc_cell_tpm,1,sd), decreasing=TRUE),][1:1000,],1,mean))

### only bivalent genes

log_sc_cell_tpm_bi <- log_sc_cell_tpm[rownames(log_sc_cell_tpm) %in% bigenes,]

tsne_out <- Rtsne(t(unique(log_sc_time_tpm[bimogenes[bimogenes %in% bigenes],])),
                  pca=FALSE,
                  perplexity=30,
                  theta=0.0)

rd1_cell <- tsne_out$Y[,1:2]

tsne_bimo <- ggplot(data=data.frame("C1"=tsne_out$Y[,1],"C2"=tsne_out$Y[,2],
           "type"=mclass),
       aes(C1,C2,col=type)) +
    geom_point()

grid.arrange(tsne_vr, tsne_bimo,ncol=2)

cl_cell <- Mclust(rd1_cell,2)$classification


heatmap.2(log_sc_cell_tpm[vargenes,],
          trace="none",
#ColSideColors=cols[as.numeric(sc_cell_coldata$exp)][order(as.numeric(a))],
#          distfun = function(x) as.dist(1 - cor(t(x), method = 'sp') ^ 2),
          col=hmcol,
          ## dist=dist,
          ## hclust=hclust,
          ColSideColors=cols[mclass],
#          Colv=FALSE,
          ## Rowv=FALSE,
#          RowSideColors=as.character(cl_bi[order(cl_bi)]),
          dendrogram = "none"
          )


bimo_mixmdl <- get_mixmdl(log_sc_cell_tpm_messup[bimocondition,])

bimo_bin <- binarizeexp(log_sc_cell_tpm[bimocondition,])

par(mfrow=c(1,4))

plot((apply(bimo_bin,1,sum)/856*100),main="all bimodal",ylim=c(0,100))

points(apply(bimo_bin[rownames(bimo_bin) %in% k4genes,],1,sum)/856*100,col="red")




boxplot(apply(bimo_bin[rownames(bimo_bin) %in% k4genes,],1,sum)/856*100,main="k4mono",ylim=c(0,100))
boxplot(apply(bimo_bin[rownames(bimo_bin) %in% k27genes,],1,sum)/856*100,main="k27mono",ylim=c(0,100))

### optimizing s.d. value

optimize_res <- optimize_noise(log_sc_cell_tpm, 20, 30)

plot_opti(optimize_res)

opti_bimonum <- apply(optimize_res,2,sum)
mysd <- as.numeric(colnames(optimize_res)[opti_bimonum == max(opti_bimonum)])

bootstrap_opti(log_sc_cell_tpm,mysd,20,30)


test <- test(cbind,optimize_noise_sub(logexp,0.1))
