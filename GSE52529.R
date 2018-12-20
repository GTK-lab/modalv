source("~/projects/bivalent/modalv/GSE75748_data_function.R")

plot_cluster <- function(expvalue,gene,colname){
    colcum <- cumsum(table(colname))
                                        #color <- ifelse(logivalue, "red", "blue")
    plot(density(expvalue[1:colcum[1]]),col="black",main=gene)
    lines(density(expvalue[(colcum[1]+1):colcum[2]]),col="blue") # todo; wrap it in apply
    lines(density(expvalue[(colcum[2]+1):colcum[3]]),col="green")
    lines(density(expvalue[(colcum[3]+1):colcum[4]]),col="red")
}


################ load data
###### data in the very first Monocle paper
df <- read.table("/mnt/cbis/home/ahjung/projects/bivalent/modalv/data/GSE52529_fpkm_matrix.txt",
                 stringsAsFactors=FALSE)

logdf <- log(as.matrix(df)+1)
logdf <- logdf[!apply(logdf,1,sum)==0,]   ## getting rid of genes with all zeroes

logdf_messup <- messup(logdf) ## this is where the normally distributed noise around zero gets added

timepoints <- sapply(colnames(df), function(x) strsplit(x,"_")[[1]][1])
names(timepoints) <- NULL
timepoints <- as.factor(timepoints)

############### Hartigan's dip test

logdf_bimodal <- apply(logdf_messup,1,test_bimodal)

#_take a look at the result
par(mfrow=c(5,10))
# not multimodal 
sapply(1:50, function(x) plot_cluster(logdf[logdf_bimodal==1,][x,],"test",timepoints))
# multimodal
sapply(1:50, function(x) plot_cluster(logdf[logdf_bimodal==0,][x,],"test",timepoints))

############### fitting mixed gaussian (parallel run)

cl <- makeCluster(detectCores()-2) # create a cluster with max-2 cores
registerDoParallel(cl) # register the cluster

res = foreach(i = 1:nrow(logdf_messup),
    .combine = "rbind",
              .packages="mixtools")     %dopar% {
  # generate a bootstrap sample              
        fit_bimodal(logdf_messup[i,],rownames(logdf_messup)[i])
}

stopCluster(cl) # shut down the cluster

res <- res[!is.na(res$gene),]

############# 

condition <- (logdf_bimodal < 0.01) & # p-value filter for dip test
    (min(res$mu1, res$mu2) < 1) & # mean of the first distribution should be near zero
        (res$lambda1 <= 0.5) # ratio of the zeroes and non zeroes

sapply(1:50, function(x) plot_cluster(logdf[condition,][x,],
                                     rownames(logdf)[condition][x],timepoints))

############ heatmap
cols_all <- palette(brewer.pal(8, "Dark2"))[as.numeric(timepoints)]
hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(20)

heatmap.2(logdf[condition,], #  & rownames(sc_cell_messup) %in% bigenes
          trace="none",
          ColSideColors=cols_all,
          distfun = function(x) as.dist(1 - cor(t(x), method = 'sp') ^ 2),
          col=hmcol,
          #Colv=F
          )

coords <- locator(1)

legend(coords,#x="topright", y=100,
    legend = unique(timepoints),
    col = unique(cols_all), 
    ## lty= 1,             
    ## lwd = 5,           
    ## cex=1
    )


########### tsne
set.seed(1)
library(Rtsne)

tsne_out <- Rtsne(t(unique(logdf[condition,])),pca=FALSE,perplexity=30,theta=0.0)

ggplot(data=data.frame("C1"=tsne_out$Y[,1],"C2"=tsne_out$Y[,2],"time"=timepoints), aes(C1,C2,col=time)) +
    geom_point()
