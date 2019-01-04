library(ggplot2)
library(gridExtra)
library(mixtools)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(diptest)
library(doParallel)
library(Rtsne)
library(mclust)
library(scatterplot3d)
library(ade4)
library(igraph)

############# Ensembl server is sometimes down. Prefer not to run this everytime.

## tx2gene <- function(){
##     mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
##                              dataset = "hsapiens_gene_ensembl",
##                              host="grch37.ensembl.org",
##                              path="/biomart/martservice")
##     t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
##                               "ensembl_gene_id",
##                               "hgnc_symbol",
##                               "refseq_mrna",
##                               "chromosome_name",
##                               "start_position",
##                               "end_position",
##                               "strand"), mart = mart)
##     t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
##                          ens_gene = ensembl_gene_id, ext_gene = hgnc_symbol,
##                          refseq = refseq_mrna,
##                          chr = chromosome_name)
##   return(t2g)			
## }
## t2g <- tx2gene()

## t2g <- cbind(t2g,t2g$end_position-t2g$start_position)
## t2g$chromosome <- as.integer(t2g$chr)
## t2g <- t2g[!is.na(t2g$chromosome),]

## t2g$chromosome2 <- paste("chr",t2g$chr,sep="")

## colnames(t2g) <- c("target_id","ens_gene","ext_gene","refseq","chrn","start","end","strand","length","chrn2","chr")

## t2g <- t2g[,c("target_id","ens_gene","ext_gene","refseq","chr","length")]

## save(t2g,file="/mnt/gtklab01/ahjung/bivalent/t2g_modalv.RData")
load("data/t2g_modalv.RData")

#================================

plot_cluster <- function(nclass,
                         expvalue,
                         gene,
                         cols,
                         xlim=c(-2,8),
                         ylim=c(0,0.6)) {

    l <- length(nclass)
    c <- cumsum(nclass)

    plot(density(expvalue[1:c[1]]),
         col=cols[1],
         xlim=xlim,
         ylim=ylim,
         main=gene,
         lwd=2)

    for (i in 2:(l-1)) {
        lines(density(expvalue[c[i]:(c[i+1]-1)]),
              col=cols[i],
              lwd=2)}
}

get_colors <- function(n) {
    return(brewer.pal(length(n),"Set2")) 
}

color_check <- function(cols,
                        coldata) {
    barplot(rep(1,length(nclass)),
            col=cols,
            names.arg=unique(coldata),
            las=2)
}

get_tpm <- function(df) {
## get gene lengths
    x <- t2g[t2g$ext_gene %in% rownames(df),]
    x <- unique(x[,c("ext_gene","length")])
    x <- x[!duplicated(x$ext_gene),]

    df2 <- df[unique(x$ext_gene),]
    xlength <- x$length/1000
    df2_rpk <- as.matrix(df2)/xlength
    rownames(df2_rpk) <- rownames(df2)

    df2_tpmfactor <- apply(df2_rpk,2, sum)

    df2_tpm <- t(t(df2_rpk) / df2_tpmfactor)
    df2_tpm <- df2_tpm[!apply(df2_tpm,1,sum)==0,]
    df2_tpm <- df2_tpm*1000000    

    return(df2_tpm)
}

test_bimodal <- function(exp) {
    a <- dip.test(as.numeric(exp))
    return(a$p.value)
}

fit_bimodal <- function(exp,name) {

    tryCatch(mixmdl <- normalmixEM(exp,
                                   k=2,
                                   mu=c(0,6),
                                   maxit = 100),
             error=function(e){NA})
    
    if (!exists("mixmdl")) {
        df <- data.frame("gene"=NA,
                         "mu1"=NA,
                         "mu2"=NA,
                         "lambda1"=NA,
                         "lambda2"=NA,
                         "converge"=NA)
    } else {
        df <- data.frame("gene"=name,
                         "mu1"=mixmdl$mu[1],
                         "mu2"=mixmdl$mu[2],
                         "lambda1"=mixmdl$lambda[1],
                         "lambda2"=mixmdl$lambda[2],
                         "converge"=rev(mixmdl$all.loglik)[1] == rev(mixmdl$all.loglik)[2] )
}

return(df)
}

messup <- function(sc_cell) {
    # sc_cell must be in log scale
    nnoise <- matrix(rnorm(nrow(sc_cell)*ncol(sc_cell),
                           mean=0, sd=0.1),ncol=ncol(sc_cell))
    sc_cell_messup <- sc_cell+nnoise
    return(sc_cell_messup)
}

fit_bimodal_multi <- function(logexp_messup,
                              ncores=detectCores()-2) {

cl <- makeCluster(ncores) # create a cluster with max-2 cores
registerDoParallel(cl) # register the cluster

res_tpm = foreach(i = 1:nrow(logexp_messup),
    .combine = "rbind",
    .packages = "mixtools",
    .export = c("fit_bimodal"))     %dopar% {
  # generate a bootstrap sample              
        fit_bimodal(logexp_messup[i,],rownames(logexp_messup)[i])
}

stopCluster(cl) # shut down the cluster
return(res_tpm)
}

filter_condition <- function(res_diptest,
                             res_tpm,
                             cond_dip=0.01,
                             cond_lambda=0.7,
                             cond_mu=1) {

    condition <- (res_diptest <= cond_dip) &
        sapply(1:nrow(res_tpm),
               function(i) min(res_tpm[i,]$mu1,
                               res_tpm[i,]$mu2) <= cond_mu) &
                   sapply(1:nrow(res_tpm),
                          function(i) max(res_tpm[i,]$mu1,
                                          res_tpm[i,]$mu2) > cond_mu)  &
                              sapply(1:nrow(res_tpm),
                                     function(i) res_tpm[i,]$lambda1 <= cond_lambda)
    return(condition)
}

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

findratio <- function(exp,idxrange) {
    binexp <- binarizeexp(exp)[idxrange]
    if (all(binexp)) { # all ON
        return(c(0,1))
    } else if (all(!binexp)) { # all OFF
        return(c(1,0))
    } else {
        return(prop.table(table(binexp*1)))
    }
}

clusterratio <- function(exp, nclass) {
    l <- length(nclass)
    c <- cumsum(nclass)

    cstart <- c(1,rev(rev(c)[-1])+1)
    cend <- c

    cratio <- t(sapply(1:l, function(x) findratio(exp, cstart[x]:cend[x])))
    rownames(cratio) <- names(nclass)
    colnames(cratio) <- c("OFF","ON")
    return(cratio)
}

windowratio <- function(exp, wsize) {
    cstart <- seq(1, length(exp), by=floor(length(exp)/wsize))
    cend <- c(cstart[-1]-1,length(exp))

    cratio <- t(sapply(1:wsize, function(x) findratio(exp, cstart[x]:cend[x])))
    rownames(cratio) <- 1:wsize
    colnames(cratio) <- c("OFF","ON")
    return(cratio)
}

plottrough <- function(exp){
    dens <- density(exp)
    trough <- findtrough(dens$y)
    plot(dens,main=rownames(exp))
    abline(v=dens$x[trough],col="red")
}

## plot_heatmap <- function(expmat,
##                          cocols=NULL, # column colors
##                          colv=FALSE) {
##     hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(20)
##     heatmap.2(expmat,
##           trace="none",
##           ColSideColors=cocols,
##           distfun = function(x) as.dist(1 - cor(t(x), method = 'sp') ^ 2),
##           col=hmcol,
##                                         Colv=colv
##           )


## }

get_bit <- function(nclass){
    bitcomb <- expand.grid(rep(list(0:1),length(nclass)))
    colnames(bitcomb) <- as.character(1:6)
    return(bitcomb)
}

expand_bit_s <- function(bitcomb_s,nclass) {
    return(unlist(sapply(1:length(nclass), function(i) rep(bitcomb_s[i],nclass[i]))))
}

expand_bit <- function(nclass) {
    bitcomb <- get_bit(nclass)
    idmat <- apply(bitcomb,1, function(i) expand_bit_s(i,nclass))
    return(idmat)
}


muin_one <- function(idmat, norm_bin_1) {
    muinx <- apply(idmat, 2,
                   function(i) aricode::NID(as.numeric(i),
                                              norm_bin_1))
    return(muinx)}

muin_all <- function(idmat, norm_bin, ncores=detectCores()-2) {

#    rnames <- rownames(norm_bin)
    cl <- makeCluster(ncores) # create a cluster with max-2 cores
    registerDoParallel(cl) # register the cluster

    res = foreach(i = 1:nrow(norm_bin),
        .combine = "rbind",
        .packages = "aricode",
        .export = c("muin_one")) %dopar% {
            muin_one(idmat, norm_bin[i,])
        }

    stopCluster(cl) # shut down the cluster
#    apply(norm_bin,1,function(i) muin_one(idmat, i))
    return(res)
}
