setwd("~/projects/bivalent/modalv")

source("GSE75748_data.R")
source("GSE75748_function.R")
source("pseudotime_function.R")

sc_cell_original <- sc_cell
sc_cell <- sc_cell[,sc_cell_coldata$cell != c("H9")]
sc_cell <- sc_cell[!apply(sc_cell,1,sum)==0,]
sc_cell_coldata_H1 <- sc_cell_coldata[sc_cell_coldata$cell != c("H9"),]
sc_cell_coldata_H1$cell <- as.factor(as.character(sc_cell_coldata_H1$cell))
sc_cell_coldata_H1$cell = factor(sc_cell_coldata_H1$cell,
    levels(sc_cell_coldata_H1$cell)[c(3,1,2,4,5,6)])

sc_cell_tpm <- get_tpm(sc_cell)
log_sc_cell_tpm <- log(sc_cell_tpm+1)

log_sc_cell_tpm_messup <- messup(log_sc_cell_tpm, 1e-05)
log_sc_cell_tpm_messup_bimodal <- apply(log_sc_cell_tpm_messup, 1, test_bimodal)

bimofit <- fit_bimodal_multi(log_sc_cell_tpm_messup,ncores=30)

bimocondition <- filter_condition(log_sc_cell_tpm_messup_bimodal,
                                  bimofit,
                                  cond_dip=0.01,
                                  cond_lambda=1,
                                  cond_mu=1)
bimocondition[is.na(bimocondition)] <- FALSE

#sum(bimocondition,na.rm=TRUE) #2604

mclass <- sc_cell_coldata_H1$cell
nclass <- table(mclass)
cols <- get_colors(nclass)

bimocondition_bi <- (rownames(log_sc_cell_tpm) %in% bigenes) & bimocondition
bimocondition_k4 <- (rownames(log_sc_cell_tpm) %in% k4genes) & bimocondition
bimocondition_k27 <- (rownames(log_sc_cell_tpm) %in% k27genes) & bimocondition

bimocondition_bi <- bimocondition_bi[!is.na(bimocondition_bi)]
bimocondition_k4 <- bimocondition_k4[!is.na(bimocondition_k4)]
bimocondition_k27 <- bimocondition_k27[!is.na(bimocondition_k27)]


# binarize

cutoff <- findcutoff_fit_all(bimofit[bimofit$gene %in% names(bimocondition)[bimocondition],])

log_sc_cell_tpm_messup_bin <- binarizeexp_fit_all(cutoff,
                         log_sc_cell_tpm_messup[names(bimocondition)[bimocondition],])


## TSNE
binmat <- log_sc_cell_tpm_messup_bin[!is.na(log_sc_cell_tpm_messup_bin[,1]),]

set.seed(100)
tsne_out_bi <- Rtsne(t(unique(binmat[rownames(binmat) %in%names(bimocondition_bi)[bimocondition_bi],])),
                  pca=FALSE,
                  perplexity=30,
                  theta=0.0)

tsne1_bi <- ggplot(data=data.frame("C1"=tsne_out_bi$Y[,1],
                    "C2"=tsne_out_bi$Y[,2],
           "type"=mclass),
       aes(C1,C2,col=type)) +
    geom_point() + theme_bw() + scale_color_manual(values=cols) + ggtitle("tSNE based on bivalent bimodal genes")

tsne_out_k4 <- Rtsne(t(unique(binmat[rownames(binmat) %in%names(bimocondition_k4)[bimocondition_k4],])),
                  pca=FALSE,
                  perplexity=30,
                  theta=0.0)

tsne1_k4 <- ggplot(data=data.frame("C1"=tsne_out_k4$Y[,1],
                    "C2"=tsne_out_k4$Y[,2],
           "type"=mclass),
       aes(C1,C2,col=type)) +
    geom_point() + theme_bw() + scale_color_manual(values=cols) + ggtitle("tSNE based on H3K4me3 mono bimodal genes")


rd1 <- tsne_out_bi$Y[,1:2]
rownames(rd1) <- mclass

library(mclust)

cl1 <- Mclust(rd1,length(names(nclass)))$classification
nclass_cluster <- table(cl1)

tsne2 <- ggplot(data=data.frame("C1"=tsne_out_bi$Y[,1],
                    "C2"=tsne_out_bi$Y[,2],
           "type"=as.factor(cl1)),
       aes(C1,C2,col=type)) +
    geom_point() + theme_bw() + scale_color_manual(values=cols) + ggtitle("6 clusters generated")


pdf("/mnt/gtklab01/ahjung/bivalent/figures/sc_cell_tsne.pdf",
    width=10, height=10)
grid.arrange(tsne1_bi,tsne1_k4,tsne2,ncol=2)
dev.off()


## ratio_per_cell <- function(nclass, exp) {
##     cstart <- c(1,rev(rev(cumsum(nclass)+1)[-1]))
##     cend <- cumsum(nclass)
##     ci <- sum(exp[cstart[1]:cend[1]])/nclass[1]
##     for (i in 2:length(nclass)) {
##         ci <- c(ci, sum(exp[cstart[i]:cend[i]])/nclass[i])
##     }
##     return(ci)
## }


## bn <- function(i, wsize) {
##     windowratio(ordered_cell[i,],table(as.numeric(cl1)),wsize)[,"ON"]    
## }

hmat_bin <- log_sc_cell_tpm_messup_bin[names(bimocondition_bi)[bimocondition_bi],]

## mim_gene_var <- sapply(1:nrow(hmat_bin),
##                    function(x) get_var(hmat_bin, x))

## mim_gene_var <-matrix(unlist(lapply(mim_gene_var,
##        function(i)
##            c(rep(NA, (nrow(hmat_bin)-length(i))),i))), nrow=nrow(hmat_bin))

## colnames(mim_gene_var) <- rownames(hmat_bin)
## rownames(mim_gene_var) <- rownames(hmat_bin)

## hr <-hclust(as.dist(mim_gene_var), method="complete")

hmat_bin_clratio <- t(apply(hmat_bin,1,function(x) clusterratio(x,cl1)))
hmat_bin_ratio <- t(apply(hmat_bin,1,function(x) split_cluster_fraction(mclass,x,32)))

tsne_out <- Rtsne(hmat_bin_ratio,
                  pca=FALSE,
                  perplexity=30,
                  theta=0.0)

rd1_wd <- tsne_out$Y[,1:2]
rownames(rd1_wd) <- rownames(hmat_bin)

cl_ratio <- Mclust(rd1_wd)$classification

pdf("/mnt/gtklab01/ahjung/bivalent/figures/sc_cell_heatmap.pdf",
    width=10, height=10)

hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(20)

heatmap.2(as.matrix(hmat_bin_ratio)[order(cl_ratio),],
          trace="none",
#ColSideColors=cols[cl1[order(cl1)]],
          col=hmcol,
          ColSideColors=bin_clusters(cols,mclass,32),
          Colv=FALSE,
          Rowv=FALSE,
          dendrogram = "none",
          RowSideColors=as.character(cl_ratio[order(cl_ratio)]),
          )

heatmap.2(log_sc_cell_tpm[rownames(hmat_bin_ratio),order(cl1)][order(cl_ratio),],
          trace="none",
ColSideColors=cols[cl1[order(cl1)]],
          col=hmcol,
          Colv=FALSE,
          Rowv=FALSE,
          dendrogram = "none",
          RowSideColors=as.character(cl_ratio[order(cl_ratio)]),
          )

dev.off()

## plot_cluster(table(cl1), log_sc_cell_tpm["SOX2",],"SOX2",cols,ylim=c(0,2),xlim=c(-2,8))

## plot_cluster(table(cl1), log_sc_cell_tpm["CER1",],"CER1",cols,ylim=c(0,1),xlim=c(-2,12))

## plot_cluster(table(cl1), log_sc_cell_tpm["DNMT3B",],"DNMT3B",cols,ylim=c(0,1),xlim=c(-2,8))

## plot_cluster(table(cl1), log_sc_cell_tpm["COCH",],"COCH",cols,ylim=c(0,0.5),xlim=c(-2,8))

## plot_cluster(table(cl1), log_sc_cell_tpm["BAG3",],"BAG3",cols,ylim=c(0,0.8),xlim=c(-2,8))

## plot_cluster(table(cl1), log_sc_cell_tpm["PRDM1",],"PRDM1",cols,ylim=c(0,0.8),xlim=c(-2,8))

## plot_cluster(table(cl1), log_sc_cell_tpm["CYP26A1",],"CYP26A1",cols,ylim=c(0,0.8),xlim=c(-2,11))



## par(mfrow=c(1,6))
## #color_check(cols,sc_cell_coldata$exp)
## color_check(cols,1:6)

## sapply(100:149, function(x)
##     plot_cluster(table(cl1), log_sc_cell_tpm[bimocondition_sub,][x,], rownames(log_sc_cell_tpm)[bimocondition_sub][x],
##                  cols))


## #rownames(ordered_cell_ratio_wd) <- rownames(ordered_cell)
## par(mfrow=c(1,6))

## sapply(c("BAG3","GATA6","NBPF15","BID","RARRES2","MRPS26"),
##        function(x) scatter.smooth(ordered_cell_ratio_wd[x,],
##                                   lpars = list(col="red", lwd=3),
##                                   main=x,
##                                   ylab="Fraction of ON",
##                                   xlab="windows"))

## #plot.new()

## cumsum(table(cl1))

## invisible(sapply(c("BAG3","BID","GATA6","NBPF15","RARRES2","MRPS26"),
##        function(x) plot(density(as.numeric(log_sc_cell_tpm[x,712:758])),main="")))

## table(cl1),log_sc_cell_tpm[x,],
##                                 x,
##                                 cols))
## color_check(cols,1:6)

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
plot_cv(log_sc_cell_tpm, TRUE, TRUE)
plot_cv(ordered_cell, TRUE, FALSE)

plot_cv(log_sc_cell_tpm[vargenes,],TRUE,FALSE)

m_cv <- get_cv(log_sc_cell_tpm,TRUE)
vargenes <- names(m_cv)[order(m_cv,decreasing=TRUE)][1:length(bimogenes)]



## boxplot(apply(log_sc_cell_tpm,1,mean))
## boxplot(apply(ordered_cell,1,mean))
## boxplot(apply(log_sc_cell_tpm[order(apply(log_sc_cell_tpm,1,sd), decreasing=TRUE),][1:1000,],1,mean))

### only bivalent genes

log_sc_cell_tpm_bi <- log_sc_cell_tpm[rownames(log_sc_cell_tpm) %in% bigenes,]

tsne_out <- Rtsne(t(unique(log_sc_cell_tpm[bimogenes[bimogenes %in% bigenes],])),
                  pca=FALSE,
                  perplexity=30,
                  theta=0.0)

rd1_cell <- tsne_out$Y[,1:2]

tsne_bimo <- ggplot(data=data.frame("C1"=tsne_out$Y[,1],"C2"=tsne_out$Y[,2],
           "type"=mclass),
       aes(C1,C2,col=type)) +
    geom_point()



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

bimo_bin <- binarizeexp(log_sc_cell_tpm[bimocondition,)]

par(mfrow=c(1,4))

plot((apply(bimo_bin,1,sum)/856*100),main="all bimodal",ylim=c(0,100))

points(apply(bimo_bin[rownames(bimo_bin) %in% k4genes,],1,sum)/856*100,col="red")




boxplot(apply(bimo_bin[rownames(bimo_bin) %in% k4genes,],1,sum)/856*100,main="k4mono",ylim=c(0,100))
boxplot(apply(bimo_bin[rownames(bimo_bin) %in% k27genes,],1,sum)/856*100,main="k27mono",ylim=c(0,100))

### optimizing s.d. value

optimize_res <- optimize_noise(log_sc_cell_tpm, 5, 30)

## pdf("/mnt/gtklab01/ahjung/Figures/modalv/optimize_noise.pdf",width=10,height=5)
## plot_opti(optimize_res)
## dev.off()

# best sd for sc_cell data is 1e-05


opti_bimonum <- apply(optimize_res,2,sum,na.rm=TRUE)
mysd <- as.numeric(colnames(optimize_res)[opti_bimonum == max(opti_bimonum)])

mysd <- 1e-05
bootstrap_num <- bootstrap_opti(log_sc_cell_tpm, mysd, 100, 20)

table(bootstrap_num)





###########3


plot_cdensity <- function(gene, nclass) {

d1 <- as.numeric(log_sc_cell_tpm[gene,])
df <- data.frame("logTPM"=d1,
"cwindow"=nclass)

q <- ggplot(df, aes(logTPM,
                    col=as.factor(cwindow),
                    fill=as.factor(cwindow))) +
                        geom_density(size=1) +
                            facet_grid(rows=vars(cwindow),
                                       scales="free") +
                                ylim(0,1) + xlim(-2,10) +
                                    theme_classic() +
                                        theme(strip.background = element_blank(),
#              strip.text.y = element_blank(),
legend.position="none") + ggtitle(gene) +
    scale_color_manual(values=cols) +
        scale_fill_manual(values=cols) + xlab("log(TPM+1)") ## + coord_flip() 

return(q)

}

plot_density <- function(gene, nclass) {

d1 <- as.numeric(log_sc_cell_tpm[gene,])
df <- data.frame("logTPM"=d1,
"cwindow"=nclass)

q <- ggplot(df, aes(logTPM)) +
                        geom_density(size=1) +
#                            facet_grid(rows=vars(cwindow)) +
xlim(-2,10) +
                                    theme_classic() +
                                        theme(strip.background = element_blank(),
#              strip.text.y = element_blank(),
legend.position="none",
                                              plot.margin = margin(0.5, 1, 0, 0.4, "cm")) + ggtitle(gene) +
    scale_color_manual(values=cols) +
        scale_fill_manual(values=cols) + xlab("log(TPM+1)")
## + coord_flip() 

return(q)

}

grob_gene <- function(gene) {
    return(list(ggplotGrob(plot_density(gene,cl1)),
                ggplotGrob(plot_cdensity(gene,cl1))))
}


pdf("/mnt/gtklab01/ahjung/bivalent/figures/sc_cell_examples.pdf",
    width=10, height=8)

grid.arrange(
    grobs = c(grob_gene("FZD7"),
grob_gene("DUSP1"),
grob_gene("TCF3"),
grob_gene("MSX1"),
grob_gene("KLF6")),
    heights = c(1,5),
  layout_matrix = rbind(c(1,3,5,7,9),
                        c(2,4,6,8,10))
)

dev.off()


pdf("/mnt/gtklab01/ahjung/bivalent/figures/sc_cell_examples2.pdf",
    width=10, height=8)


grid.arrange(
    grobs = c(grob_gene("TMEM106C"),
grob_gene("FZD7"),
grob_gene("LGALS3"),
grob_gene("PSAT1")
),
    heights = c(1,5),
  layout_matrix = rbind(c(1,3,5,7),
                        c(2,4,6,8))
)
dev.off()

