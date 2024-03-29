library(ggplot2)
library(gridExtra)
library(mixtools)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(diptest)
library(doParallel)
library(Rtsne)

########### import and process data

mydir <- "/mnt/gtklab01/ahjung/bivalent-tmp/GSE75748/"

## bulk_cell <- read.csv(paste0(mydir,"GSE75748_bulk_cell_type_ec.csv"),row.names=1)
## bulk_time <- read.csv(paste0(mydir,"GSE75748_bulk_time_course_ec.csv"),row.names=1)
sc_cell <- read.csv(paste0(mydir,"GSE75748_sc_cell_type_ec.csv"),row.names=1)
sc_time <- read.csv(paste0(mydir,"GSE75748_sc_time_course_ec.csv"),row.names=1)

extract_coldata_sub <- function(df,i){
df_split <- strsplit(colnames(df)[i],"_|[.]")[[1]]

coldata <- data.frame("cell"=df_split[1],
                              "exp"=df_split[2],
                      "idx"=df_split[3])
return(coldata)}

extract_coldata <- function(df) {

coldata <- extract_coldata_sub(df,1)
invisible(sapply(2:ncol(df),
       function(i) coldata <<- rbind(coldata,
                                     extract_coldata_sub(df,i))))
return(coldata)}

sc_cell_coldata <- extract_coldata(sc_cell)
sc_time_coldata <- extract_coldata(sc_time)
## bulk_time_coldata <- extract_coldata(bulk_time)
## bulk_cell_coldata <- data.frame("cell"=strsplit(colnames(bulk_cell)[1],"_")[[1]][1],
                              ## "rep"=strsplit(colnames(bulk_cell)[1],"_")[[1]][2])
## invisible(sapply(2:ncol(bulk_cell),
##        function(i) bulk_cell_coldata <<- rbind(bulk_cell_coldata,
##                                                data.frame("cell"=strsplit(colnames(bulk_cell)[i],"_")[[1]][1],
##                                                           "rep"=strsplit(colnames(bulk_cell)[i],"_")[[1]][2]))))

#hist(log(as.numeric(sc_cell[1,sc_cell_coldata$cell == "H1"])+1))

get_scexp <- function(sc_cell,sc_cell_coldata,celltype,mark,markid) {

    h1id <- unique(t2g$ext_gene[t2g$target_id %in% get(mark)])
#    h1id <- unique(t2g$ext_gene[t2g$target_id %in% a])
    h1id <- h1id[!h1id == ""]

    x <- log(as.numeric(unlist(
        sc_cell[rownames(sc_cell) %in% h1id,sc_cell_coldata$cell %in% celltype]))+1)
    
    cat("equal/below 2 is", length(x[x<=2]), "above 2 is", length(x[x>2]),"\n")

    x2 <- data.frame("exp"=x[x>2], ######### filter any counts below 2
                     "class"=rep(markid,length(x[x>2])))

return(x2) }

get_bulkexp <- function(celltype,mark,markid) {

    h1id <- unique(t2g$ext_gene[t2g$target_id %in% get(mark)])
#    h1id <- unique(t2g$ext_gene[t2g$target_id %in% a])
    h1id <- h1id[!h1id == ""]

    x <- apply(bulk_cell[rownames(bulk_cell) %in% h1id,
                         bulk_cell_coldata$cell %in% celltype],1,mean)
    x <- log(as.numeric(unlist(x))+1)

x2 <- data.frame("exp"=x[x>2], "class"=rep(markid,length(x[x>2])))

return(x2)}


#=====================================

## sc_test_tb <- ggplot(allexp, aes(x=exp, col=class)) + geom_density(size=2)+
##     ggtitle("sc: TB ; chip: BMP4T")

## grid.arrange(h1, bmp4m_dec, bmp4t, test,ncol=1)
## grid.arrange(bmp4t, test, h1,ncol=1)
## grid.arrange(bulk_test, bulk_test_tb, sc_test, sc_test_tb,ncol=2)
## grid.arrange(bulk_test, sc_test, ncol=1)




## plot_hist("EC","BMP4M_1_bivalent"),ylim=c(0,0.3),col="black")
## lines(plot_hist("EC","BMP4M_1_H3K4me3"),col="blue")
## lines(plot_hist("EC","BMP4M_1_H3K27me3"),col="red")
## #plot_hist(c("DEC","EC"),"H1hesc_bivalent")




## get_scexp2 <- function(celltype,mark, markid) {

## h1id <- unique(t2g$ext_gene[t2g$target_id %in% mark])
## h1id <- h1id[!h1id == ""]

## x <- log(as.numeric(unlist(
##     sc_cell[rownames(sc_cell) %in% h1id,sc_cell_coldata$cell %in% celltype]))+1)

## x2 <- data.frame("exp"=x[x>2], "class"=rep(markid,length(x[x>2])))

## return(x2)}

## a <- hESC_1_bivalent[hESC_1_H3K27me3 %in% BMP4M_1_bivalent]
## b <- hESC_1_bivalent[hESC_1_H3K27me3 %in% BMP4M_1_H3K4me3]
## c <- hESC_1_bivalent[hESC_1_H3K27me3 %in% BMP4M_1_H3K27me3]

## allexp <- rbind(get_scexp2("DEC", a,"bivalent"),
##                 get_scexp2("DEC", b,"H3K4me3"),
##                 get_scexp2("DEC", c,"H3K27me3"))

## test4a <- ggplot(allexp, aes(x=exp, col=class)) + geom_density(size=2)+
##     ggtitle("hESC H3K27me3 genes in BMP4M")

## grid.arrange(test2,test3,test4,
## test2a,test3a,test4a,ncol=3)

## df <- data.frame("x"=bulk_cell[bulk_cell[,5] > 0,5])
## ggplot(log(df),aes(x=x)) + geom_density()

## df <- data.frame("x"=sc_cell[sc_cell[,3] > 0,3])

## grid.arrange(ggplot(log(df),aes(x=x)) + geom_density(),
## h1,ncol=1)


############# mixed gaussian model

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
load("/mnt/gtklab01/ahjung/bivalent/t2g_modalv.RData")

#================================
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

plot_mixmdl <- function(mixmdl, mark) {

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
                        ggtitle(paste0(mark,"; mu is ",
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

get_mixmdl <- function(allexp, mark) {

set.seed(1)
wait <- allexp$exp[allexp$class %in% mark]
mixmdl <- normalmixEM(wait, k = 2)

return(mixmdl)
}

plot_mg <- function(allexp, mark) {
    mixmdl <- get_mixmdl(allexp, mark)
    p <- plot_mixmdl(mixmdl, mark)
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


plot_cluster <- function(nclass,expvalue,gene,cols,xlim=c(-2,8),ylim=c(0,0.6)){

length(nclass)
cumsum(nclass)

#cols <- brewer.pal(length(nclass),"Set1")
plot(density(expvalue[1:cumsum(nclass)[1]]),
     col=cols[1],
     xlim=xlim,
     ylim=ylim,
     main=gene)

for (i in 2:(length(nclass)-1)) {
    lines(density(expvalue[cumsum(nclass)[i]:(cumsum(nclass)[i+1]-1)]),
          col=cols[i])
}

#color <- ifelse(logivalue, "red", "blue")
## plot(density(expvalue[1:212]),col="black",xlim=c(-2,10),ylim=c(0,1),main=gene)
## lines(density(expvalue[213:350]),col="red")
## lines(density(expvalue[351:455]),col="blue")
## lines(density(expvalue[456:524]),col="green")
}

get_colors <- function(n) {
    return(brewer.pal(length(n),"Set2"))
}

color_check <- function(cols, coldata) {
    barplot(rep(1,length(nclass)),
            col=cols,
            names.arg=unique(coldata),
            las=2)
}

return_density <- function(expvalue_sub,name){
    exptest <- density(expvalue_sub)
    return(data.frame("x"=exptest$x,"y"=exptest$y,"name"=name))
}

plot_time <- function(expvalue,gene,coldatatbl) {
    ncate <- c(0,as.numeric(cumsum(table(coldatatbl))))
    tblnames <- names(table(coldatatbl))

    densitytbl <- return_density(expvalue[(ncate[1]+1):ncate[1+1]],tblnames[1])
    
    sapply(2:(length(ncate)-1), function(x)
        densitytbl <<- rbind(densitytbl,return_density(expvalue[(ncate[x]+1):ncate[x+1]],tblnames[x])))

p <- ggplot(data=densitytbl, aes(x=x,y=y, col=name)) + geom_line() + facet_wrap(~name,ncol=1)
return(p)
    ## plot(density(expvalue[1:ncate[2]]),
##          col="black",xlim=c(-2,10),ylim=c(0,1),main=gene)
## invisible(    sapply(2:(length(ncate)-1), function(x)
##         lines(density(expvalue[(ncate[x]+1):ncate[x+1]]),col=cols[x])))

}

exp_cluster <- function(expvalue,gene){
#color <- ifelse(logivalue, "red", "blue")
return(data.frame("gene"=gene,
                  "H1"=median(expvalue[1:212]),
                  "DEC"=median(expvalue[213:350]),
                  "EC"=median(expvalue[351:455]),
                  "TB"=median(expvalue[456:524])))
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

## exp <- sc_cell
## a <- exp[,1:212]

## c <- apply(a,1,test_bimodel)

## cellfactors <- as.numeric(sc_cell_coldata$cell[sc_cell_coldata$cell %in% c("H1","DEC","EC","TB")])
## cellfactors[cellfactors==3] <- 2
## cellfactors[cellfactors==4] <- 3
## cellfactors[cellfactors==7] <- 4


fit_bimodal <- function(exp,name) {

    tryCatch(mixmdl <- normalmixEM(exp, k=2, mu=c(0,6),maxit = 100),error=function(e){NA})
    
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

fit_bimodal_multi <- function(logexp_messup) {

cl <- makeCluster(detectCores()-2) # create a cluster with max-2 cores
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

filter_condition <- function(res_diptest, res_tpm, cond_dip=0.01, cond_lambda=0.7, cond_mu=1) {
## Default conditions are 0.01,0.7,1
    condition <- (res_diptest <= cond_dip) &
        sapply(1:nrow(res_tpm),
               function(i) min(res_tpm[i,]$mu1, res_tpm[i,]$mu2) <= cond_mu) &
                   sapply(1:nrow(res_tpm),
                          function(i) max(res_tpm[i,]$mu1, res_tpm[i,]$mu2) > cond_mu)  &
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

plottrough <- function(exp){
    dens <- density(exp)
    trough <- findtrough(dens$y)
    plot(dens,main=rownames(exp))
    abline(v=dens$x[trough],col="red")
}

