library(ggplot2)
library(gridExtra)
library(mixtools)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(diptest)

mydir <- "/mnt/gtklab01/ahjung/bivalent-tmp/GSE75748/"

bulk_cell <- read.csv(paste0(mydir,"GSE75748_bulk_cell_type_ec.csv"),row.names=1)
bulk_time <- read.csv(paste0(mydir,"GSE75748_bulk_time_course_ec.csv"),row.names=1)
sc_cell <- read.csv(paste0(mydir,"GSE75748_sc_cell_type_ec.csv"),row.names=1)
sc_time <- read.csv(paste0(mydir,"GSE75748_sc_time_course_ec.csv"),row.names=1)

sc_cell_coldata <- data.frame("cell"=strsplit(colnames(sc_cell)[1],"_|[.]")[[1]][1],
                              "exp"=strsplit(colnames(sc_cell)[1],"_|[.]")[[1]][2],
                              "idx"=strsplit(colnames(sc_cell)[1],"_|[.]")[[1]][3])

invisible(sapply(2:ncol(sc_cell),
       function(i) sc_cell_coldata <<- rbind(sc_cell_coldata,
                                             data.frame("cell"=strsplit(colnames(sc_cell)[i],"_|[.]")[[1]][1],
                                                        "exp"=strsplit(colnames(sc_cell)[i],"_|[.]")[[1]][2],
                                                        "idx"=strsplit(colnames(sc_cell)[i],"_|[.]")[[1]][3]))))

sc_time_coldata <- data.frame("cell"=strsplit(colnames(sc_time)[1],"_|[.]")[[1]][1],
                              "exp"=strsplit(colnames(sc_time)[1],"_|[.]")[[1]][2],
                              "idx"=strsplit(colnames(sc_time)[1],"_|[.]")[[1]][3])

invisible(sapply(2:ncol(sc_time),
       function(i) sc_time_coldata <<- rbind(sc_time_coldata,
                                             data.frame("cell"=strsplit(colnames(sc_time)[i],"_|[.]")[[1]][1],
                                                        "exp"=strsplit(colnames(sc_time)[i],"_|[.]")[[1]][2],
                                                        "idx"=strsplit(colnames(sc_time)[i],"_|[.]")[[1]][3]))))

bulk_time_coldata <- data.frame("cell"=strsplit(colnames(bulk_time)[1],"_")[[1]][1],
                              "exp"=strsplit(colnames(bulk_time)[1],"_")[[1]][2],
                              "rep"=strsplit(colnames(bulk_time)[1],"_")[[1]][3])

invisible(sapply(2:ncol(bulk_time),
       function(i) bulk_time_coldata <<- rbind(bulk_time_coldata,
                                             data.frame("cell"=strsplit(colnames(bulk_time)[i],"_")[[1]][1],
                                                        "exp"=strsplit(colnames(bulk_time)[i],"_")[[1]][2],
                                                        "rep"=strsplit(colnames(bulk_time)[i],"_")[[1]][3]))))

bulk_cell_coldata <- data.frame("cell"=strsplit(colnames(bulk_cell)[1],"_")[[1]][1],
                              "rep"=strsplit(colnames(bulk_cell)[1],"_")[[1]][2])


invisible(sapply(2:ncol(bulk_cell),
       function(i) bulk_cell_coldata <<- rbind(bulk_cell_coldata,
                                               data.frame("cell"=strsplit(colnames(bulk_cell)[i],"_")[[1]][1],
                                                          "rep"=strsplit(colnames(bulk_cell)[i],"_")[[1]][2]))))

hist(log(as.numeric(sc_cell[1,sc_cell_coldata$cell == "H1"])+1))


get_scexp <- function(sc_cell,sc_cell_coldata,celltype,mark,markid) {

    h1id <- unique(t2g$ext_gene[t2g$target_id %in% get(mark)])
#    h1id <- unique(t2g$ext_gene[t2g$target_id %in% a])
    h1id <- h1id[!h1id == ""]

    x <- log(as.numeric(unlist(
        sc_cell[rownames(sc_cell) %in% h1id,sc_cell_coldata$cell %in% celltype]))+1)
    
    cat("equal/below 2 is", length(x[x=<2]), "above 2 is", length(x[x>2]),"\n")

    x2 <- data.frame("exp"=x[x>2],
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

return(x2)


}

load("/mnt/gtklab01/ahjung/bivalent/bivalent-workspace.RData")

hESC_others <- t2g$target_id[!t2g$target_id %in% unique(c(hESC_H3K4me3, hESC_H3K27me3, hESC_bivalent))]

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


plot_cluster <- function(expvalue,gene){

#color <- ifelse(logivalue, "red", "blue")
plot(density(expvalue[1:212]),col="black",xlim=c(-2,10),ylim=c(0,1),main=gene)
lines(density(expvalue[213:350]),col="red")
lines(density(expvalue[351:455]),col="blue")
lines(density(expvalue[456:524]),col="green")
}

return_density <- function(expvalue_sub,name){
    exptest <- density(expvalue_sub)
    return(data.frame("x"=exptest$x,"y"=exptest$y,"name"=name))
}

plot_time <- function(expvalue,gene,coldatatbl) {
    ncate <- c(0,as.numeric(cumsum(table(coldatatbl))))
    tblnames <- names(table(coldatatbl))

    densitytbl <-         return_density(expvalue[(ncate[1]+1):ncate[1+1]],tblnames[1])
    
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



test_bimodal <- function(exp) {

a <- dip.test(log(as.numeric(exp)+1))
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

    tryCatch(mixmdl <- normalmixEM(log(as.numeric(exp)+1), k=2, mu=c(0,6),maxit = 10),error=function(e){NA})
    
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


messup <- function(exp,zeromin,zeromax) {
## The idea is to add a little variance to zeros
if (((sum(exp==0)/length(exp)*100) > 30) & ((sum(exp==0)/length(exp)*100) < 80)) {
    exp[sample(which(log(as.numeric(exp)+1)==0),floor(sum(exp==0)*0.5))] <- 0.3 } else { exp } # <- rep(NA,length(exp)) }
return(exp)
}


