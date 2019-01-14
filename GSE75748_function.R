## for mixed model
library(diptest)
library(mixtools)
## for visualization
library(dplyr)
library(ggplot2)
library(gridExtra)
library(gplots)
library(RColorBrewer)
## for multi-thread processing
library(doParallel)
## for clustering and graph structures
library(Rtsne)
library(mclust)
library(igraph)
library(scatterplot3d) # optional
library(ade4) # optional


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

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

plot_mixmdl <- function(mixmdl) {

p <- data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 0.1, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.5) + xlim(0,12.5) +
                    ylab("Density") +
                        ggtitle(paste0("mu is ",
                                       round(mixmdl$mu[1],digits=2),
                                       ", ",
                                       round(mixmdl$mu[2],digits=2),
                                       "; lambda is ",
                                       round(mixmdl$lambda[1],digits=2),
                                       ", ",
                                       round(mixmdl$lambda[2],digits=2))) +
                                           theme_bw()
return(p)
}

get_mixmdl <- function(allexp) {

set.seed(1)
mixmdl <- normalmixEM(allexp, k = 2)

return(mixmdl)
}

plot_mg <- function(allexp) {
    mixmdl <- get_mixmdl(allexp)
    p <- plot_mixmdl(mixmdl)
    return(p)
}

raincloud_theme <- theme(
    text = element_text(size = 10),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 14),
#    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.title=element_text(size=16),
    legend.text=element_text(size=16),
    legend.position = "right",
    plot.title = element_text(lineheight=.8, face="bold", size = 16),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))


plot_cluster <- function(nclass,
                         expvalue,
                         gene,
                         cols,
                         xlim=c(-2,8),
                         ylim=c(0,0.6)) {
    l <- length(nclass)
    c <- cumsum(nclass)
    cstart <- c(1,rev(rev(c)[-1])+1)
    cend <- c

    plot(density(expvalue[cstart[1]:cend[1]]),
         col=cols[1],
         xlim=xlim,
         ylim=ylim,
         main=gene,
         lwd=2)

    for (i in 2:l) {
        lines(density(expvalue[cstart[i]:cend[i]]),
              col=cols[i],
              lwd=2)}
}

return_density <- function(expvalue_sub,
                           name){
    exptest <- density(expvalue_sub)
    return(data.frame("x"=exptest$x,"y"=exptest$y,"name"=name))
}

plot_time <- function(nclass,
                      expvalue,
                      gene) {
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    l <- length(nclass)
    c <- cumsum(nclass)
    cstart <- c(1,rev(rev(c)[-1])+1)
    cend <- c

    d <- return_density(expvalue[cstart[1]:cend[1]],
                        names(nclass)[1])
    
    for (i in 2:l) {
        d <- rbind(d,return_density(expvalue[cstart[i]:cend[i]],
                            names(nclass)[i]))}

    p <- ggplot(data=d,
                aes(x=x, y=y, col=name)) +
                    geom_line() +
                        facet_wrap(~name,ncol=1) +
                            scale_colour_manual(values=cbPalette)
return(p)
    ## plot(density(expvalue[1:ncate[2]]),
##          col="black",xlim=c(-2,10),ylim=c(0,1),main=gene)
## invisible(    sapply(2:(length(ncate)-1), function(x)
##         lines(density(expvalue[(ncate[x]+1):ncate[x+1]]),col=cols[x])))

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

fit_bimodal <- function(exp,
                        name) {

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

messup <- function(sc_cell, sd) {
    # sc_cell must be in log scale
    nnoise <- matrix(rnorm(nrow(sc_cell)*ncol(sc_cell),
                           mean=0, sd=sd),ncol=ncol(sc_cell))
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

fit_bimodal_one <- function(logexp_messup) {
    res_tpm <- fit_bimodal(logexp_messup[1,], rownames(logexp_messup)[1])

    sapply(2:nrow(logexp_messup), function(i)
        res_tpm <<- rbind(res_tpm,
                          fit_bimodal(logexp_messup[i,], rownames(logexp_messup)[i])))
return(res_tpm)
}

run_multi <- function(logexp, n, ncores=detectCores()-2) {

    cl <- makeCluster(ncores) # create a cluster with max-2 cores
    registerDoParallel(cl) # register the cluster

    optimize_res = foreach(i = 1:length(n),
        .combine = "cbind",
        .packages = "diptest",
        .export = c("messup",
            "dip.test",
            "optimize_noise_sub",
            "test_bimodal",
            "fit_bimodal",
            "fit_bimodal_one",
            "filter_condition")) %dopar% {
            optimize_noise_sub(logexp,n[i])
        }

    stopCluster(cl) # shut down the cluster

    ## optimize_res <- data.frame(optimize_noise_sub(logexp,n[1]),
    ##                            row.names=rownames(logexp))
    ## sapply(n[-1],
    ##        function(x) optimize_res <<- cbind(optimize_res,
    ##                                           optimize_noise_sub(logexp,x)))

    colnames(optimize_res) <- as.character(n)
    return(optimize_res)
}

optimize_noise_sub <- function(logexp, sd, ncores) {
    cat("noise model with sd= ",sd,"\n")
    logexp_m <- messup(logexp, sd)
    logexp_m_bimodal <- apply(logexp_m,
                              1,
                              test_bimodal)

    bimofit <- fit_bimodal_multi(logexp_m, ncores)
    bimocondition <- filter_condition(logexp_m_bimodal,
                                      bimofit)
    return(bimocondition)
}

optimize_noise <- function(logexp, niter=30, ncores=detectCores()-2) {
n <- c(0,1e-11,1e-10,1e-09,1.5e-08,1e-08,1e-07,1e-06,1.5e-5,1e-05,1.5e-04,1e-04,1e-03,1e-02,1e-01,2e-01,3e-01,4e-01,5e-01)
#    n <-seq(0.2,1,length.out=niter)
    optimize_res <- optimize_noise_sub(logexp, n[1], ncores)
    sapply(n[-1],
           function(x) optimize_res <<- cbind(optimize_res,
                                              optimize_noise_sub(logexp,x,ncores)))
    colnames(optimize_res) <- as.character(n)
    return(optimize_res)
}

plot_opti <- function(optimize_res) {
    optimize_num <- apply(optimize_res,2,sum,na.rm=TRUE)
    plot(as.numeric(colnames(optimize_res)),
         optimize_num,
         ylab="Number Of Identified Bimodal Genes",
         xlab="S.D. of Introduced Noise")
    lines(as.numeric(colnames(optimize_res)),
         optimize_num, col="red")
    ## scatter.smooth(1:length(optimize_num),optimize_num,
    ##                lpars = list(col="red", lwd=3))
    ##                ## main=,
    ##                ## ylab="Fraction of ON",
    ##                ## xlab="windows")
}

bootstrap_opti <- function(logexp, sd=0.1, nbootstrap=30, ncores=detectCores()-2) {
    n <- rep(sd,nbootstrap)
    bootstrap_res <- optimize_noise_sub(logexp, n[1], ncores)
    sapply(n[-1],
           function(x) bootstrap_res <<- cbind(bootstrap_res,
                                              optimize_noise_sub(logexp,x,ncores)))
    colnames(bootstrap_res) <- as.character(n)

#    bootstrap_num <- apply(bootstrap_res,1,sum)/nbootstrap
    return(bootstrap_num)
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
    if (is.na(cutoff)) { return(NA) # TODO: need a warning for NA
                     } else {
    return(exp > cutoff) }
}

findratio <- function(exp,idxrange) {
    binexp <- binarizeexp(exp)[idxrange]
    binexp <- binexp[!is.na(binexp)]
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

window_per_cluster <- function(csize, wsize, cinit) {
    # csize: size of cluster, wsize: size of window
    cend <- cumsum(rep(wsize, csize%/%wsize))
    if (rev(cend)[1] != csize) { cend <- c(cend, csize) }
    cstart <- c(1, rev(rev(cend)[-1])+1)
    clusterw <- data.frame("cstart"=cstart,"cend"=cend)
    clusterw <- clusterw + cinit
    return(clusterw)
}

window_all_clusters <- function(nclass, wsize) {
    clusterw <- window_per_cluster(nclass[1], wsize, 0)
    invisible(sapply(1:length(nclass),
           function(i)
               clusterw <<- rbind(clusterw, window_per_cluster(nclass[i], wsize, cumsum(nclass)[i-1]))))
           return(clusterw)
}

bin_clusters <- function(nclass,wsize) {
    clusterw <- window_all_clusters(nclass,wsize)
    b <- rep(1:length(nclass),
             diff(c(0,which(sapply(clusterw$cend, function(x) x %in% cumsum(nclass))))))
    return(b)
}

find_windowratio <- function(exp, clusterw) {
    ## cstart <- seq(1, length(exp), by=floor(length(exp)/wsize))
    ## cend <- c(cstart[-1]-1,length(exp))
    cstart <- clusterw$cstart
    cend <- clusterw$cend
    cratio <- t(sapply(1:nrow(clusterw), function(x) findratio(exp, cstart[x]:cend[x])))
    rownames(cratio) <- 1:nrow(clusterw)
    colnames(cratio) <- c("OFF","ON")
    return(cratio)
}

windowratio <- function(exp, nclass, wsize) {
    clusterw <- window_all_clusters(nclass, wsize)
    cratio <- find_windowratio(exp, clusterw)
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
