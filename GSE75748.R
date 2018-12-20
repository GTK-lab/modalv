source("~/projects/bivalent/modalv/GSE75748_data_function.R")

## sc_test_tb <- ggplot(allexp, aes(x=exp, col=class)) + geom_density(size=2)+
##     ggtitle("sc: TB ; chip: BMP4T")

## grid.arrange(h1, bmp4m_dec, bmp4t, test,ncol=1)
## grid.arrange(bmp4t, test, h1,ncol=1)
## grid.arrange(bulk_test, bulk_test_tb, sc_test, sc_test_tb,ncol=2)
## grid.arrange(bulk_test, sc_test, ncol=1)




## plot_hist("EC","BMP4M_1_bivalent")
## ,ylim=c(0,0.3),col="black")


lines(plot_hist("EC","BMP4M_1_H3K4me3"),col="blue")
lines(plot_hist("EC","BMP4M_1_H3K27me3"),col="red")
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



## allexp <- rbind(get_scexp(sc_cell_messup,sc_cell_coldata[sc_cell_coldata$cell %in% c("H1","DEC","EC","TB"),],"H1", "hESC_1_bivalent", "bivalent"),
##                 get_scexp(sc_cell_messup,sc_cell_coldata[sc_cell_coldata$cell %in% c("H1","DEC","EC","TB"),],"H1", "hESC_1_H3K4me3", "H3K4me3"),
##                 get_scexp(sc_cell_messup,sc_cell_coldata[sc_cell_coldata$cell %in% c("H1","DEC","EC","TB"),],"H1", "hESC_1_H3K27me3", "H3K27me3"),
## get_scexp(sc_cell_messup,sc_cell_coldata[sc_cell_coldata$cell %in% c("H1","DEC","EC","TB"),],"H1", "hESC_others", "others"))

## plot(density(log(allexp$exp+1)))
## mixmdl <- normalmixEM(log(allexp$exp+1), k = 3)

## h1_k4 <- plot_mg(allexp, "H3K4me3")
## h1_k27 <- plot_mg(allexp, "H3K27me3")
## h1_bi <- plot_mg(allexp, "bivalent")
## h1_ot <- plot_mg(allexp, "others")

## grid.arrange(h1_k4, h1_k27, h1_bi, h1_ot, ncol=1)

## grid.arrange(h1_k4,h1_k27,h1_bi,h1_ot,
## bmp4m_k4,bmp4m_k27,bmp4m_bi,bmp4m_ot,
## bmp4t_k4,bmp4t_k27,bmp4t_bi, bmp4t_ot,ncol=4)

## BMP4M_H3K4me3 <- BMP4M_1_H3K4me3[BMP4M_1_H3K4me3 %in% BMP4M_2_H3K4me3]
## BMP4M_H3K27me3 <- BMP4M_1_H3K27me3[BMP4M_1_H3K27me3 %in% BMP4M_2_H3K27me3]
## BMP4M_bivalent <- BMP4M_1_bivalent[BMP4M_1_bivalent %in% BMP4M_2_bivalent]

## BMP4M_others <- t2g$target_id[!t2g$target_id %in%
## unique(c(BMP4M_1_H3K4me3, BMP4M_2_H3K4me3,
##   BMP4M_1_H3K27me3, BMP4M_2_H3K27me3,
##   BMP4M_1_bivalent, BMP4M_2_bivalent))]

k4genes <- unique(t2g$ext_gene[t2g$target_id %in% hESC_H3K4me3])
k4genes <- k4genes[!k4genes == ""]
k27genes <- unique(t2g$ext_gene[t2g$target_id %in% hESC_H3K27me3])
k27genes <- k27genes[!k27genes == ""]
bigenes <- unique(t2g$ext_gene[t2g$target_id %in% hESC_bivalent])
bigenes <- bigenes[!bigenes == ""]


######################### BULK

allexp <- rbind(get_bulkexp(c("EC","DEC"), "BMP4M_bivalent","bivalent"),
                get_bulkexp(c("EC","DEC"), "BMP4M_H3K4me3","H3K4me3"),
                get_bulkexp(c("EC","DEC"), "BMP4M_H3K27me3","H3K27me3"),
                get_bulkexp(c("EC","DEC"), "BMP4M_others","others"))

bulk_test_es <- ggplot(allexp, aes(y=exp,x=class, col=class)) +geom_jitter(alpha=0.3) + geom_violin(size=2,alpha=0.5)+    ggtitle("bulk: EC & DEC ; chip: BMP4M") + theme_bw() + raincloud_theme


g1 <- list(bmp4m_k4 + raincloud_theme,
bulk_test_es,
bmp4m_k27+ raincloud_theme,
bmp4m_bi+ raincloud_theme,
bmp4m_ot+ raincloud_theme)

grid.arrange(
  grobs = g1,
  widths = c(1, 1),
  layout_matrix = rbind(c(2, 1),
      c(2,3),
      c(2,4),
      c(2,5))
)

############################################## can we distinguish real zeros from noise
sc_cell_original <- sc_cell

sc_cell <- sc_cell[,sc_cell_coldata$cell != c("H9")]
sc_cell <- sc_cell[!apply(sc_cell,1,sum)==0,]
logsc_cell <- log(as.matrix(sc_cell)+1)
nnoise <- matrix(rnorm(nrow(sc_cell)*ncol(sc_cell),
                       mean=0, sd=0.1),ncol=ncol(sc_cell))

logsc_cell_messup <- logsc_cell+nnoise



#sc_cell <- sc_cell_original
## sc_cell_messup <- t(apply(sc_cell,1,messup))
#sc_cell_messup <- sc_cell_messup[!is.na(sc_cell_messup[,1]),]

## sc_cell_zero <- sc_cell[apply(sc_cell,1,function(x) any(x==0)),]
## sc_cell_nozero <- sc_cell[!apply(sc_cell,1,function(x) any(x==0)),]



#sc_cell_bimodal <- apply(sc_cell,1,test_bimodal)
sc_cell_bimodal <- apply(logsc_cell_messup,1,test_bimodal)

## sc_cell_fit <- fit_bimodal(sc_cell[1,],rownames(sc_cell[1,]))
## sapply(2:nrow(sc_cell),function(x) sc_cell_fit <<- rbind(sc_cell_fit,fit_bimodal(sc_cell[x,],rownames(sc_cell[x,]))))
## sc_cell_fit <- sc_cell_fit[!is.na(sc_cell_fit$gene),]

logsc_cell_messup_fit <- fit_bimodal(logsc_cell_messup[1,],rownames(logsc_cell_messup)[1])
invisible(sapply(2:nrow(logsc_cell_messup),function(x) logsc_cell_messup_fit <<- rbind(logsc_cell_messup_fit,fit_bimodal(logsc_cell_messup[x,],rownames(logsc_cell_messup)[x]))))
logsc_cell_messup_fit <- logsc_cell_messup_fit[!is.na(logsc_cell_messup_fit$gene),]

condition <- (sc_cell_bimodal < 0.01) &
    (min(logsc_cell_messup_fit$mu1, logsc_cell_messup_fit$mu2) < 1) &
        (logsc_cell_messup_fit$lambda1 < 0.6)

k4condition <- condition & rownames(logsc_cell_messup) %in% k4genes
bicondition <- condition & rownames(logsc_cell_messup) %in% bigenes
k27condition <- condition & rownames(logsc_cell_messup) %in% k27genes

cat(rownames(logsc_cell)[k4condition])

par(mfrow=c(5,10))

sapply(1:50, function(x)    plot_cluster(logsc_cell_messup[bicondition,][x,],rownames(logsc_cell_messup)[bicondition][x]))


## plot_cluster <- function(expvalue,gene){
## #color <- ifelse(logivalue, "red", "blue")
## lpot(density(expvalue[1:933]),col="black",xlim=c(-2,4),ylim=c(0,1),main=gene)
## lines(density(expvalue[934:1731]),col="red")
## }



## sapply(1:50, function(x)
##     plot_cluster(log(as.numeric(sc_cell_messup[condition,][x,])+1),
##                  rownames(sc_cell_messup)[condition][x]))


## sapply(1:50, function(x)
##     plot_cluster(log(as.numeric(sc_cell_messup[condition,][x,])+1),
                 ## rownames(sc_cell_messup)[x]))



## o6 <- ggplot(data=data.frame("exp"=log(as.numeric(sc_cell_messup[condition,][6,])+1),
##            "cell"=as.factor(as.character(sc_cell_coldata$cell[sc_cell_coldata$cell != "H9"]))), aes(x=exp,col=cell)) +
##                geom_density()

## grid.arrange(o1,o2,o3,o4,o5,o6,ncol=2)

plot(apply(log(sc_cell+1),1,mean),(apply(log(sc_cell+1),1,sd)^2/apply(log(sc_cell+1),1,mean)),col=  adjustcolor( "black", alpha.f = 0.3),pch=16)

points(apply(log(sc_cell[bicondition,]+1),1,mean),(apply(log(sc_cell[bicondition,]+1),1,sd)^2/apply(log(sc_cell[bicondition,]+1),1,mean)),col=adjustcolor("green",alpha.f=0.1),pch=16)

## plot_cluster_vio <- function(exp, name) {
##  p <- ggplot(data.frame("exp"=log(as.numeric(exp)+1),"cell"=sc_cell_coldata$cell[sc_cell_coldata$cell %in% c("H1","DEC","EC","TB")]),
##             aes(cell, exp,fill=cell)) + geom_violin()
## return(p)
## }
## heatmap

#cols_all <- palette(brewer.pal(8, "Dark2"))[as.numeric(sc_cell_coldata$cell)]
cols_all <- palette(brewer.pal(8, "Dark2"))[as.numeric(as.factor(as.character(sc_cell_coldata$cell[sc_cell_coldata$cell != "H9"])))]
## cols_all <- palette(brewer.pal(8, "Dark2"))[as.numeric(sc_cell_coldata$exp)]

## cols_all <- palette(brewer.pal(8, "Dark2"))[c(rep(1,933),rep(2,1731-933))]

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(50)

heatmap.2(logsc_cell_messup[k4condition,], #  & rownames(sc_cell_messup) %in% bigenes
trace="none",ColSideColors=cols_all,
#RowSideColors=cols2,
col=hmcol,
#                                        Colv=F
          )

legend(y=0.9, x=-4, xpd=TRUE,     
    legend = unique(sc_cell_coldata$exp),#[sc_cell_coldata$cell != "H9"]))),
    col = unique(cols_all), 
    lty= 1,             
    lwd = 5,           
    cex=1
    )


heatmap.2(log(as.matrix(sc_cell_messup[as.character(rownames(sc_cell_messup)[condition & rownames(sc_cell_messup) %in% k4genes]),])+1),
trace="none",ColSideColors=cols_all,
#RowSideColors=cols2,
col=hmcol,
#                                        Colv=F
          )

heatmap.2(log(as.matrix(sc_cell_original[rownames(sc_cell_original) %in% (rownames(sc_cell)[condition& rownames(sc_cell_messup) %in% bigenes]),])+1),
trace="none",ColSideColors=cols_all,
#RowSideColors=cols2,
col=hmcol,
#                                        Colv=F
          )



## x <- data.frame("mean"=c(apply(log(as.matrix(sc_cell_messup[as.character(rownames(sc_cell_messup)[condition & rownames(sc_cell_messup) %in% bigenes]),])+1),1,mean),apply(log(as.matrix(sc_cell_messup[as.character(rownames(sc_cell_messup)[condition & rownames(sc_cell_messup) %in% k4genes]),])+1),1,mean)),
## "class"=c(rep("bivalent",length(apply(log(as.matrix(sc_cell_messup[as.character(rownames(sc_cell_messup)[condition & rownames(sc_cell_messup) %in% bigenes]),])+1),1,mean))),rep("k4mono",length(apply(log(as.matrix(sc_cell_messup[as.character(rownames(sc_cell_messup)[condition & rownames(sc_cell_messup) %in% k4genes]),])+1),1,mean)))))

## ggplot(data=x, aes(x=class, y=mean, fill=class))+
##     geom_boxplot()



set.seed(1)
tsne_out <- Rtsne(t(unique(sc_cell[as.character(rownames(sc_cell)[bicondition]),])),pca=FALSE,perplexity=30,theta=0.0)

## tsne_out <- Rtsne(t(unique(sc_cell_original[rownames(sc_cell_original) %in% (rownames(sc_cell)[condition]),])),pca=FALSE,perplexity=30,theta=0.0)






## plot(tsne_out$Y,col=as.factor(as.character(sc_cell_coldata$cell[sc_cell_coldata$cell != "H9"])),asp=1)

## ggplot(data=data.frame("C1"=tsne_out$Y[,1],"C2"=tsne_out$Y[,2],"class"=class), aes(C1,C2,col=class)) +
##     geom_point()

ggplot(data=data.frame("C1"=tsne_out$Y[,1],"C2"=tsne_out$Y[,2],"type"=as.factor(as.character(sc_cell_coldata$cell[sc_cell_coldata$cell != "H9"]))), aes(C1,C2,col=type)) +
    geom_point()


ggplot(data=data.frame("C1"=tsne_out$Y[,1],"C2"=tsne_out$Y[,2],"type"=sc_cell_coldata$exp), aes(C1,C2,col=type)) +
    geom_point()

class <- rep("NA",nrow(unique(sc_cell[as.character(rownames(sc_cell)[condition]),])))
class[rownames(unique(sc_cell[as.character(rownames(sc_cell)[condition]),])) %in% k4genes] <- "K4"
class[rownames(unique(sc_cell[as.character(rownames(sc_cell)[condition]),])) %in% bigenes] <- "BI"
class[rownames(unique(sc_cell[as.character(rownames(sc_cell)[condition]),])) %in% k27genes] <- "K27"
class <- as.factor(class)

tsne_out <- Rtsne(t(unique(sc_cell[as.character(rownames(sc_cell)[condition]),])),pca=FALSE,perplexity=30,theta=0.0)




plot(density(log(as.matrix(sc_cell_original)[as.matrix(sc_cell_original)>2]+1)))

plot(
log(apply(as.matrix(sc_cell),1,mean)+1),log((apply(as.matrix(sc_cell),1,sd)/apply(as.matrix(sc_cell),1,mean))+1),pch=16,col=adjustcolor("black",alpha.f=0.1))


points(
log(apply(as.matrix(sc_cell),1,mean)+1)[condition],log((apply(as.matrix(sc_cell),1,sd)/apply(as.matrix(sc_cell),1,mean))+1)[condition],pch=16,col=adjustcolor("red",alpha.f=0.3))


####### BASICS
counts_original <- counts


devtools::install_github("catavallejos/BASiCS", build_vignettes = TRUE)

library(BASiCS)
## counts2 <- as.matrix(sc_cell)
## counts <- matrix(as.integer(counts2),nrow=nrow(counts2))
## rownames(counts) <- rownames(counts2)
counts <- esd0
counts <- counts[!apply(counts,1,sum)==0,]

Tech = rep(FALSE,nrow(counts))
set.seed(2)
## SpikeInput = rgamma(10,1,1)
## SpikeInfo <- data.frame("SpikeID" = paste0("Spike", 1:10), "SpikeInput" = SpikeInput)

library(SingleCellExperiment)
batchinfo <- as.numeric(data.frame("x"=sapply(as.character(colnames(counts)), function(x) strsplit(x,"[.]")[[1]][1]))$x)
## batchinfo <- as.numeric(table(sc_cell_coldata[sc_cell_coldata$cell != "H9",c(1,2)]))[!as.numeric(table(sc_cell_coldata[sc_cell_coldata$cell != "H9",c(1,2)]))==0][c(1,2,3,4,8,5,9,6,10,13,7,11,12)]

## batchinfoidx <- 1:length(batchinfo)

## batchinfo <- unlist(sapply(1:length(batchinfo), function(x) rep(x,batchinfo[x])))

DataNoSpikes <- SingleCellExperiment(assays = list(counts = counts),
                                    colData = data.frame(BatchInfo = batchinfo))

ChainNoSpikes <- BASiCS_MCMC(Data = DataNoSpikes, N = 1000, 
                             Thin = 10, Burn = 500, 
                             WithSpikes = FALSE,  Regression = FALSE,
                             PrintProgress = FALSE)

plot(ChainNoSpikes, Param = "mu", Gene = 1, log = "y")
plot(ChainNoSpikes, Param = "phi", Cell = 1)

par(mfrow = c(1,2))
HVG <- BASiCS_DetectHVG(ChainNoSpikes, VarThreshold = 0.9, Plot = TRUE)
LVG <- BASiCS_DetectLVG(ChainNoSpikes, VarThreshold = 0.2, Plot = TRUE)

HVG$Table[,"GeneName"] 


plot(apply(log(sc_cell+1),1,mean),(apply(log(sc_cell+1),1,sd)^2/apply(log(sc_cell+1),1,mean)),col=  adjustcolor( "black", alpha.f = 0.3),pch=16)

points(apply(log(sc_cell[condition,]+1),1,mean),(apply(log(sc_cell[condition,]+1),1,sd)^2/apply(log(sc_cell[condition,]+1),1,mean)),col=adjustcolor("red",alpha.f=0.3),pch=16)


plot(apply(log(sc_cell+1),1,mean),(apply(log(sc_cell+1),1,sd)^2/apply(log(sc_cell+1),1,mean)),col=  adjustcolor( "black", alpha.f = 0.3),pch=16)



points(apply(log(sc_cell[rownames(sc_cell) %in% HVG$Table[,"GeneName"],]+1),1,mean),(apply(log(sc_cell[rownames(sc_cell) %in% HVG$Table[,"GeneName"],]+1),1,sd)^2/apply(log(sc_cell[rownames(sc_cell) %in% HVG$Table[,"GeneName"],]+1),1,mean)),col=adjustcolor("red",alpha.f=0.3),pch=16)



#######################################
esd0 <- read.csv("/mnt/gtklab01/ahjung/bivalency_dropseq/mES/data/GSM1599494_ES_d0_main.csv", header=FALSE, stringsAsFactors=FALSE, row.names=1)
colnames(esd0) <- paste0("d0.", seq_len(ncol(esd0)))

head(counts)

poisregmixEM(counts[1,], lambda = NULL, beta = NULL, k = 2,
             addintercept = TRUE, epsilon = 1e-08, 
             maxit = 10000, verb = FALSE)


set.seed(100)
beta <- matrix(c(1, .5, .7, -.8), 2, 2)
x <- runif(50, 0, 10)
xbeta <- cbind(1, x)%*%beta
w <- rbinom(50, 1, .5)
y <- w*rpois(50, exp(xbeta[, 1]))+(1-w)*rpois(50, exp(xbeta[, 2]))
out <- poisregmixEM(y, x, verb = TRUE,  epsilon = 1e-03)
out

batchinfo <- c(rep(1,933),rep(2,))
