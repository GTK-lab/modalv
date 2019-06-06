setwd("~/projects/bivalent/modalv")

source("GSE75748_data.R")
source("GSE75748_function.R")

load("data/20190601.RData")
sc_time_original <- sc_time
sc_time <- sc_time[!apply(sc_time,1,sum)==0,]

sc_time_tpm <- get_tpm(sc_time) # some genes get lost
log_sc_time_tpm <- log(sc_time_tpm+1)

log_sc_time_tpm_messup <- messup(log_sc_time_tpm, 1e-05)
log_sc_time_tpm_messup_bimodal <- apply(log_sc_time_tpm_messup, 1, test_bimodal)

bimofit <- fit_bimodal_multi(log_sc_time_tpm_messup,ncores=50)

bimocondition <- filter_condition(log_sc_time_tpm_messup_bimodal,
                                  bimofit,
                                  cond_dip=0.01,
                                  cond_lambda=1,
                                  cond_mu=1)

mclass <- sc_time_coldata$exp
nclass <- table(mclass)
cols <- get_colors(nclass)



bimocondition[is.na(bimocondition)] <- FALSE

bimocondition_bi <- (rownames(log_sc_time_tpm) %in% bigenes) & bimocondition
bimocondition_k4 <- (rownames(log_sc_time_tpm) %in% k4genes) & bimocondition
bimocondition_k27 <- (rownames(log_sc_time_tpm) %in% k27genes) & bimocondition

bimocondition_bi <- bimocondition_bi[!is.na(bimocondition_bi)]
bimocondition_k4 <- bimocondition_k4[!is.na(bimocondition_k4)]
bimocondition_k27 <- bimocondition_k27[!is.na(bimocondition_k27)]

library("Rtsne")
library("ggplot2")

binmat <- log_sc_time_tpm_messup_bin[!is.na(log_sc_time_tpm_messup_bin[,1]),]

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

rd1 <- tsne_out_bi$Y[,1:2]
rownames(rd1) <- mclass

library(mclust)

cl1 <- Mclust(rd1,6)$classification
## cl1[cl1==4] <- 6
## cl1[cl1==5] <- 4
## cl1[cl1==6] <- 5
#plot(rd1_pc, pch=16, asp = 1,col=cols[cl1])

tsne2 <- ggplot(data=data.frame("C1"=tsne_out_bi$Y[,1],
                    "C2"=tsne_out_bi$Y[,2],
           "type"=as.factor(cl1)),
       aes(C1,C2,col=type)) +
    geom_point() + theme_bw() + scale_color_manual(values=cols) + ggtitle("6 clusters generated")



library(gridExtra)

pdf("/mnt/gtklab01/ahjung/bivalent/figures/sc_time_tsne.pdf",    width=5, height=8)
grid.arrange(tsne1_bi, tsne2)
dev.off()

save.image("20190601.RData")

# binarize

log_sc_time_tpm_messup_bin <- binarizeexp(log_sc_time_tpm_messup[bimocondition,][1,])

invisible(sapply(2:sum(bimocondition),
                 function(i) log_sc_time_tpm_messup_bin <<- rbind(log_sc_time_tpm_messup_bin,
binarizeexp(log_sc_time_tpm_messup[bimocondition,][i,]))))

rownames(log_sc_time_tpm_messup_bin) <- names(bimocondition)[bimocondition]

#bi_bin <- log_sc_time_tpm_messup_bin[names(bimocondition_bi)[bimocondition_bi],]*1



## my.mst <- mstree(dist(rd1),1)
## mygraph <- graph_from_adjacency_matrix(neig2mat(my.mst))
## mygraph <- minimum.spanning.tree(mygraph)
## x <- get.shortest.paths(mygraph,1,758)
## a <- shortest.paths(mygraph, 1) # starting from any time from time 0
## a2 <- get.shortest.paths(mygraph,1) # need to keep the option to specify end point

## longest <- which(unlist(lapply(a2$vpath,length))==max(unlist(lapply(a2$vpath,length))))

## V(mygraph)$color <- cols[sc_time_coldata$exp]
## E(mygraph, path=a2$vpath[[longest]])$color <- "red"
## V(mygraph)[as.numeric(a2$vpath[[longest]])]$color <- "black"

## pdf("/mnt/gtklab01/ahjung/bivalent/figures/sc_time_pseudotime.pdf",
##     width=8, height=10)

## par(mfrow=c(1,1))
## plot(mygraph,
##      vertex.size=2,
## #     layout=rd1,
##      vertex.label=NA,
##      edge.arrow.size=0.5)

## dev.off()




library(infotheo)
library(minet)
library(Rgraphviz)

varinformation <- function(X,Y,method="emp") {
    condentropy(X,Y,method) + condentropy(Y,X,method)
}

get_var <- function(mat, idx) {
    x <- sapply(c(1:nrow(mat)),
                function(x) varinformation(mat[idx,], mat[x,]))
    return(x)
}

get_mi <- function(mat, idx) {
    x <- sapply(c(1:nrow(mat)),
                function(x) mutinformation(mat[idx,], mat[x,]))
    return(x)
}

hmmcondition <- (bimocondition_k27 | bimocondition_bi) & !bimocondition_k4

hmat <- matrix(as.numeric(log_sc_time_tpm[hmmcondition,]),
               nrow=sum(hmmcondition))
rownames(hmat) <- rownames(log_sc_time_tpm[hmmcondition,])
colnames(hmat) <- colnames(log_sc_time_tpm)


hmat_bin <- log_sc_time_tpm_messup_bin[names(hmmcondition)[hmmcondition],]*1
hmat_bin <- hmat_bin[!(is.na(hmat_bin[,1])),]


# mutual information matrix
mim_cell <- sapply(1:ncol(hmat_bin),
                   function(x) get_mi(t(hmat_bin), x))
colnames(mim_cell) <- colnames(hmat_bin)
rownames(mim_cell) <- colnames(hmat_bin)

mim_gene <- sapply(1:nrow(hmat_bin),
                   function(x) get_mi(hmat_bin, x))
colnames(mim_gene) <- rownames(hmat_bin)
rownames(mim_gene) <- rownames(hmat_bin)

# conditional entropy matrix as distance between cells
mim_cell_var <- sapply(1:ncol(hmat_bin),
                   function(x) get_var(t(hmat_bin), x))
colnames(mim_cell_var) <- colnames(hmat_bin)
rownames(mim_cell_var) <- colnames(hmat_bin)



## get_dpt <- function(cluster) {
##     covars <- data.frame("cells"=colnames(mim_cell)[cl1%in%cluster])
##     dm <- DiffusionMap(covars,                       distance=as.dist(mim_cell_var[cl1%in%cluster, cl1%in%cluster]))
##     dpt <- DPT(dm)
##     return(dpt)
## }

nclass_cluster <- table(cl1)



## ### within cluster order

## find_neig_tips <- function(mim_tips) {
##     mim_tips_m <- melt(mim_tips)
##     df <- mim_tips_m[order(mim_tips_m$value,decreasing=TRUE)[1:2],]
## #    which(colnames(mim_cell)==find_neig_tips(mim_tip_00_12)[1])
##     return(df)
## }

## mim_tip_2 <- mim_cell[tips(dpt_2)+cumsum(nclass_cluster)[1],
##                       tips(dpt_3)+cumsum(nclass_cluster)[3]]
## mim_tip_3 <- mim_cell[tips(dpt_3)+cumsum(nclass_cluster)[3],
##                       tips(dpt_4)+cumsum(nclass_cluster)[3]]
## #mim_tip_4_5 <- mim_cell[tips(dpt_4)+cumsum(nclass_cluster)[3], tips(dpt_5)+cumsum(nclass_cluster)[4]]


## tips_list <- list(find_neig_tips(mim_tip_2),
##                   find_neig_tips(mim_tip_3))
 

## ctips <- c("H9.24h_034", "H9.36h_128", "H9.36h_059","H9.36h_119")

## order_1 <- order(as.numeric(mim_cell[cl1==1,"H9.24h_034"]),decreasing=TRUE) # use how similar the cell is to the next tip and order in decreasing order
## order_2 <- order(dpt_2@dm@eigenvectors[,1],decreasing=FALSE)
## order_3 <- order(dpt_3@dm@eigenvectors[,1],decreasing=FALSE)



## cellorder <- c(order_1,
##                order_2+cumsum(nclass_cluster)[1],
##                order_3+cumsum(nclass_cluster)[3])

################### pseudotime ordering based on PST

get_connectingpair <- function(firstcluster,
                               secondcluster) {
    ## identify cells with the shortest varinfo between two clusters
    sub_var <- mim_cell_var[cl1 %in% c(firstcluster),
                            cl1 %in% c(secondcluster)]
    sub_var_order <- sub_var[order(sub_var)]
    sub_pairs <- sapply(1:(nrow(sub_var)*ncol(sub_var)),
                        function(x)
                            which(sub_var == sub_var_order[x],
                                  arr.ind = TRUE))
    sub_pairs <- sub_pairs[!duplicated(sub_pairs)]
    sub_pairs_df <- data.frame(
        "fromcluster"=firstcluster,
        "tocluster"=secondcluster,
        "fromcell" = as.numeric(unlist(lapply(sub_pairs,
            function(x)
                x[,1]))),
        "tocell" = as.numeric(unlist(lapply(sub_pairs,
            function(x)
                x[,2]))))
    return(sub_pairs_df)
}

get_cellorder <- function(firstcellidx,
                          cluster) {
## based on selected first cell, order the rest of the cells in the cluster
    mim_cell_var_sub <- mim_cell_var[cl1%in%c(cluster),
                                     cl1%in%c(cluster)]
    prev <- ifelse(cluster == 1,
                   0,
                   cumsum(nclass_cluster)[cluster-1])
    clusterorder <- order(mim_cell_var_sub[firstcellidx,],
                          decreasing=FALSE) #+ prev
    return(clusterorder)
}

get_ordered_cells <- function(cluster_orders) {
## order cells based on given ordering
    return(apply(cluster_orders,
                               2,
                 function(x) c(1:ncol(cluster_orders))[x]))
}

build_PST <- function(cluster_orders, cluster, n=100, cluster_orders2) {
    blackcol <- colorRampPalette(c("black","white"))

    if(missing(cluster_orders2)) {
        ordered_cells <- generate_pairorders(cluster_orders, n=100)
    } else {
        ordered_cells <- generate_pairorders(cluster_orders, n=100,
                                             cluster_orders2=cluster_orders2)
    }
    
    ## orderings used to build PST
    cellseq_sample <-TraMineR::seqdef(ordered_cells,
                                      labels = colnames(sc_time)[cl1%in%cluster],
                                      cpal=blackcol(ncol(cluster_orders)))

    C1 <- PST::pstree(cellseq_sample,
                      weighted=TRUE)
    return(C1)
}

predict_order <- function(C1, ordered_cell, cluster, n=100) {
        blackcol <- colorRampPalette(c("black","white"))
    ## orderings used to predict logloss
    cellseq_all <-TraMineR::seqdef(ordered_cell,
                                   labels = colnames(sc_time)[cl1%in%cluster],
                                   cpal=blackcol(ncol(ordered_cell)))

    C1_logloss <- PST::predict(C1,
                               cellseq_all,
                               output = "logloss")
    return(C1_logloss)
}

predicted_order <- function(tree, ordered_cell, cluster, n=100) {
    C1_logloss  <- predict_order(tree,
                                 ordered_cell,
                                 cluster,
                                 n=n)
    cellorder <- ordered_cell[which.min(C1_logloss),]
    cat("predicted order is from ",which.min(C1_logloss),"\n")
    return(cellorder)
}

generate_pairorders <- function(cluster_orders, n=100, cluster_orders2) {
    if(missing(cluster_orders2)) {
        ordered_cells_all <- get_ordered_cells(cluster_orders)[1:n,]
    } else {
        ordered_cells_all <- get_ordered_cells(cluster_orders)[1:n,]
        ordered_cells2_all <- get_ordered_cells(cluster_orders2)[1:n,]
        ordered_cells_all <- rbind(ordered_cells_all, ordered_cells2_all)
    }
    
    return(ordered_cells_all)

}

## fix_pos <- function(ordered_cells_all, n=1000) {
##     pos <-    apply(ordered_cells_all, 2,
##           function(x)
##               names(table(x)[which.max(table(x))]))
##     fixed_pos <- as.numeric(names(which(table(pos)==1)))
##     cellID <- seq(1,ncol(ordered_cells_all))
##     cellsequence <- data.frame(cellID,
##                                order=pos, stringsAsFactors=FALSE)
##     cellsequence$order[!pos%in%fixed_pos] <- NA

## cellsequence_shuffled <-  do.call( rbind,
##             replicate(n, shuffle_orders(cellsequence),
##                                          simplify=FALSE))
## return(cellsequence_shuffled)
        
## }

## shuffle_orders <- function(cellsequence) {
## df <- data.frame(cellID=cellsequence$cellID[is.na(cellsequence$order)],
##                  order=sample(c(1:nrow(cellsequence))[!1:length(cellID) %in% cellsequence$order]), stringsAsFactors=FALSE)

## df_shuffled <- rbind(cellsequence[!is.na(cellsequence$order),],
##                      df)

## return(as.numeric(df_shuffled$order))
## }

from1to2 <- get_connectingpair(1,2)
from2to3 <- get_connectingpair(2,3)
from3to4 <- get_connectingpair(3,4)
from4to5 <- get_connectingpair(4,5)
from5to6 <- get_connectingpair(5,6)

cluster1_order <- t(sapply(from1to2$fromcell, function(x) rev(get_cellorder(x,1))))
cluster2_order_from1 <- t(sapply(from1to2$tocell, function(x) get_cellorder(x,2)))
cluster2_order_from3 <- t(sapply(from2to3$fromcell, function(x) rev(get_cellorder(x,2))))
cluster3_order_from2 <- t(sapply(from2to3$tocell, function(x) get_cellorder(x,3)))
cluster3_order_from4 <- t(sapply(from3to4$fromcell, function(x) rev(get_cellorder(x,3))))
cluster4_order_from3 <- t(sapply(from3to4$tocell, function(x) get_cellorder(x,4)))
cluster4_order_from5 <- t(sapply(from4to5$fromcell, function(x) rev(get_cellorder(x,4))))
## cluster5_order_from4 <- t(sapply(from4to5$tocell, function(x) get_cellorder(x,5)))
## cluster5_order_from6 <- t(sapply(from5to6$fromcell, function(x) get_cellorder(x,5)))
## cluster6_order <- t(sapply(from5to6$tocell, function(x) get_cellorder(x,6)))
cluster5_order <- t(sapply(from4to5$tocell, function(x) get_cellorder(x,c(5,6))))


C1.pst <- build_PST(cluster1_order, 1, n=100)
C2.pst <- build_PST(cluster2_order_from1, 2, n=100,
                    cluster_orders2=cluster2_order_from3)
C3.pst <- build_PST(cluster3_order_from2, 3, n=100,
                    cluster_orders2=cluster3_order_from4)
C4.pst <- build_PST(cluster4_order_from3, 4, n=100,
                    cluster_orders2=cluster4_order_from5)
## C5.pst <- build_PST(cluster5_order_from4, 5, n=100,
##                     cluster_orders2=cluster5_order_from6)
## C6.pst <- build_PST(cluster6_order, 6, n=100)

C5_6.pst <- build_PST(cluster5_order, c(5,6), n=100)

generate_pst <- function(c.pst, vertices, cluster) {
    seq_length <- sum(cl1%in%cluster)
    ## gen.con <- TRUE
    ## while (gen.con) {
    c.gen <- PST::generate(c.pst,
                  seq_length,
                  1,
                  as.character(vertices[1]),
sum(vertices==vertices[1])/length(vertices),
                  method="prob")

    return(c.gen) }

library(foreach)
library(doParallel)

repeat_gen <- function(c.pst, vertices, cluster, nr, np=50) {
registerDoParallel(np)  # use multicore, set to the number of our cores
c.gen <- foreach (i=1:nr, .combine=rbind) %dopar% {
    cat(i,"\n")
    generate_pst(c.pst, vertices, cluster) }
return(c.gen)
}

C1.gen <- repeat_gen(C1.pst, from1to2[1:100,3],1,10,10)
C2.gen <- repeat_gen(C2.pst,
                     c(generate_pairorders(cluster2_order_from1,
                                           n=100)[,1],
                       generate_pairorders(cluster2_order_from3,n=100)[,as.numeric(nclass_cluster[2])]),
                     2,
                     10,
                     10)

C3.gen <- repeat_gen(C3.pst,
                     c(generate_pairorders(cluster3_order_from2, n=100)[,1],
                       generate_pairorders(cluster3_order_from4, n=100)[,as.numeric(nclass_cluster[3])]),
                     3,
                     10,
                     10)

C4.gen <- repeat_gen(C4.pst,
                     c(generate_pairorders(cluster4_order_from3, n=100)[,1],
                       generate_pairorders(cluster4_order_from5, n=100)[,as.numeric(nclass_cluster[4])]),
                     4,
                     10,
                     10)

C5_6.gen <- repeat_gen(C5_6.pst,
                       from4to5[1:100,4],
                       c(5,6),
                       10,
                       10)


## run_generate_pst <- function(c.pst,
##                              vertices,
##                              cluster,
##                              maxn=10) {
##     c.gen <- as.numeric(generate_pst(c.pst,
##                                      vertices,
##                                      cluster))

##     repeat {
##         c.gen <- rbind(c.gen,
##                        generate_pst(c.pst,
##                                     vertices,
##                                     cluster))
##         bcon <-  !any(duplicated(as.numeric(c.gen[nrow(c.gen),])))
##         if((nrow(c.gen)>=maxn)|(bcon)) { break }}
##     return(c.gen[!any(duplicated(as.numeric(c.gen[nrow(c.gen),]))),]) }

get_unique_order <- function(c.gen) {
    return(as.numeric(c.gen[which(apply(c.gen, 1, function(x) sum(duplicated(x))) == 0)[1],]))
}



which(apply(C5_6.gen, 1, function(x) sum(duplicated(x)))==0)

C1.order <- predicted_order(C1.pst, generate_pairorders(cluster1_order), 1)

C2.order <- as.numeric(C2.gen[which(apply(C2.gen, 1, function(x) sum(duplicated(x)))==0)[1],])
C3.order <- as.numeric(C3.gen[which(apply(C3.gen, 1, function(x) sum(duplicated(x)))==0)[1],])
C4.order <- as.numeric(C4.gen[which(apply(C4.gen, 1, function(x) sum(duplicated(x)))==0)[1],])
## C5.order <- as.numeric(C5.gen[which(apply(C5.gen, 1, function(x) sum(duplicated(x)))==0)[1],])
## C6.order <- as.numeric(C6.gen[which(apply(C6.gen, 1, function(x) sum(duplicated(x)))==0)[1],])
C5_6.order <- as.numeric(C5_6.gen[which(apply(C5_6.gen,1,function(x) sum(duplicated(x)))==0)[1],])

cellorder <- c(C1.order,
               C2.order+cumsum(nclass_cluster)[1],
               C3.order+cumsum(nclass_cluster)[2],
               C4.order+cumsum(nclass_cluster)[3],
               C5_6.order+cumsum(nclass_cluster)[4])

               ## C5.order+cumsum(nclass_cluster)[4],
               ## C6.order+cumsum(nclass_cluster)[5])



write.table(cellorder, "/mnt/gtklab01/ahjung/bivalent/results/sc_time_cellorder.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")
###################








hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(20)

hmat <- matrix(as.numeric(log_sc_time_tpm[bimocondition_bi,]),
               nrow=sum(bimocondition_bi))
rownames(hmat) <- rownames(log_sc_time_tpm[bimocondition_bi,])
colnames(hmat) <- colnames(log_sc_time_tpm)


## compare before and after cell ordering

pdf("/mnt/gtklab01/ahjung/bivalent/figures/compare_cellorder.pdf",
width=10, height=10)

heatmap.2(hmat_bin[,cellorder],
          trace="none",
ColSideColors=cols[sc_time_coldata$exp][cellorder],
                                        #as.numeric(a))],#cols[cl1][order(cl1)],
          col=hmcol,
          Colv=F,
#         Rowv=F,
          dendrogram = "none"
          )

dev.off()


#########################3 HMM
library(HMM)

init_HMM <- function(obs) {
    return(initHMM(States=c("OFF","ON"),
            Symbols=c("0","1"),
            startProbs=as.numeric(table(obs)/length(obs)),
            transProb=matrix(c(.99,.01,.01,.99),nrow=2),
            emission=matrix(c(.9,.1,.1,.9),nrow=2)))}

    
posteriors_bw <- function(obs) {
    bw <- baumWelch(init_HMM(obs), obs)
    print("done")
    return(as.numeric(HMM::posterior(bw$hmm,obs)[2,]))
}


posteriors_viterbi <- function(obs) {
    vit <- viterbiTraining(init_HMM(obs), obs)
    return(as.numeric(HMM::posterior(vit$hmm,obs)[2,]))
}


viterbi_path <- function(obs) {
    vit <- viterbiTraining(init_HMM(obs), obs)
    return(viterbi(vit$hmm,obs))
    }

run_viterbi <- function(exp, name) {
    obs <- as.character(exp)
    tryCatch(vit <- posteriors_viterbi(obs),
             error=function(e) { cat(name,"\n") })
    if (!exists("vit")) {
        vitc <- data.frame(rep(NA, length(obs)))
    } else {
        vitc <- data.frame(vit)
    }
    colnames(vitc) <- name
    return(vitc) }

run_viterbi_path <- function(exp, name) {
    obs <- as.character(exp)
    tryCatch(vit <- viterbi_path(obs),
             error=function(e) { cat(name,"\n") })
    if (!exists("vit")) {
        vitc <- data.frame(rep(NA, length(obs)))
    } else {
        vitc <- data.frame(vit)
    }
    colnames(vitc) <- name
    return(vitc) }


##### get matrix for HMM

hmat_bin_cellorder <- hmat_bin[,cellorder]

vit_df <- data.frame(sapply(1:nrow(hmat_bin_cellorder),
                            function(x) run_viterbi(hmat_bin_cellorder[x,],
                                                    rownames(hmat_bin_cellorder)[x])))
vit_df_path <- data.frame(sapply(1:nrow(hmat_bin_cellorder),
                            function(x) run_viterbi_path(hmat_bin_cellorder[x,],
                                                    rownames(hmat_bin_cellorder)[x])))

vit_idx <- (!is.na(vit_df[1,])) & (!is.na(vit_df_path[1,]))

vit_mat <- as.matrix(vit_df[,vit_idx])
rownames(vit_mat) <- colnames(hmat_bin)[cellorder]

vit_mat_path <- as.matrix(vit_df_path[,vit_idx])
rownames(vit_mat_path) <- colnames(hmat_bin)[cellorder]

## intupgenes <- colnames(vit_mat_path)[(vit_mat_path[1,]=="OFF")&(vit_mat_path[758,]=="OFF")]
##                                         #colnames(vit_mat)[(vit_mat[1,]<0.1)&(vit_mat[758,]<0.1)]
## upgenes <- colnames(vit_mat_path)[(vit_mat[1,]<0.001)&(vit_mat[758,]>0.999)]
## downgenes <- colnames(vit_mat_path)[(vit_mat[1,]>0.999)&(vit_mat[758,]<0.001)]


############# find connectors
mim_con <- mim_gene[colnames(mim_gene) %in% downstream,
                    !colnames(mim_gene) %in% downstream]

connectors <- unique(unlist(apply(mim_con, 1, function(x) names(which(x>=0.1)))))
connectors <- connectors[connectors %in% colnames(vit_mat)]
connectors <- connectors[!connectors %in% as.character(unique(gsa_mat_m_overlap$gene))]


find_switch <- function(genex, onoroff) {
    switch_on <- which(diff(as.numeric(vit_mat_path[,genex]=="ON")) == 1)+1
    switch_off <- which(diff(as.numeric(vit_mat_path[,genex]=="ON")) == (-1))+1
    return(ifelse(onoroff=="ON", data.frame("ON"=switch_on), data.frame("OFF"=switch_off)))
}

on_switches <- sapply(colnames(vit_mat_path),
                      function(x) find_switch(x, "ON"))

off_switches <- sapply(colnames(vit_mat_path),
                      function(x) find_switch(x, "OFF"))

filter_genes <- function(df) {
    start_bin <- as.logical(df[1])
    end_bin <- as.logical(df[2])
    on_switch_n <- as.numeric(df[3])
    off_switch_n <- as.numeric(df[4])
    start <- ifelse(start_bin, "ON","OFF")
    end <- ifelse(end_bin, "ON","OFF")
    return(colnames(vit_mat_path)[(vit_mat_path[1,]==start)&
                                      (vit_mat_path[758,]==end)&
                                          (as.numeric(lapply(on_switches, length)) == on_switch_n) &
                                              (as.numeric(lapply(off_switches, length) == off_switch_n))])}

gene_switches_df <- data.frame("genes"=names(on_switches),
                               "start_bin"=vit_mat_path[1,]=="ON",
                               "end_bin"=vit_mat_path[758,]=="ON",
                               "on_switches"=as.numeric(unlist(lapply(on_switches, length))),
                               "off_switches"=as.numeric(unlist(lapply(off_switches, length))),
                               stringsAsFactors=FALSE)

gene_switches_df_table <- as.data.frame(table(gene_switches_df[,-1]))
gene_switches_df_table <-  gene_switches_df_table[gene_switches_df_table$Freq != 0,]
gene_switches_df_table <-  cbind(gene_switches_df_table,
                                 "switch_group"= 1:nrow(gene_switches_df_table))

switches_genes <- apply(gene_switches_df_table[,1:4],1,filter_genes)

which_switch_group <- function(gene) {
    df <- gene_switches_df[gene_switches_df$genes == gene,]
    return(gene_switches_df_table[((gene_switches_df_table$start_bin == df$start_bin) &
         (gene_switches_df_table$end_bin == df$end_bin) &
             (gene_switches_df_table$on_switches == df$on_switches) &
                 (gene_switches_df_table$off_switches == df$off_switches)),"switch_group"])
}


which_change_window <- function(gene) {
    x <- change_genes[as.character(change_genes$gene) == gene,"change_window"]
    return(x)
}


order_by_switch <- function(switches_df_idx, whichswitch) {
    cat(switches_df_idx,"\n")

    gswhen <- as.numeric(unlist(lapply(whichswitch[switches_genes[[switches_df_idx]]],
                                       function(x) x[1])))

    gsorder <- order(gswhen)
    gsdf <- data.frame("gene"=switches_genes[[switches_df_idx]][gsorder],
                       "when"=gswhen[gsorder],
                       "when_cell"=colnames(hmat_bin_cellorder)[gswhen[gsorder]],
                       "switch_group"=switches_df_idx,
                       "switch_group_order"=1:length(gswhen))
    
    return(gsdf)
}


order_by_cluster <- function(switches_df_idx) {
    cat(switches_df_idx,"\n")

    sgenes <- switches_genes[[switches_df_idx]]

    gswhen <- as.numeric(unlist(lapply(on_switches[switches_genes[[switches_df_idx]]],
                                       function(x) ifelse(x[1]> table(cl1)[1], x[1], x[2]))))

    if (!length(sgenes) < 2) {
        colgenes <- gsub("-","[.]",colnames(vit_mat))
        
        intupclust <- hclust(dist(t(vit_mat)[colgenes %in% sgenes,]))
        gsorder <- intupclust$label[intupclust$order]
        gsdf <- data.frame("gene"=gsorder,
                           "when"=gswhen[intupclust$order],
                           "when_cell"=colnames(hmat_bin_cellorder)[gswhen[intupclust$order]],
                           "switch_group"=switches_df_idx,
                           "switch_group_order"=1:length(gsorder))
} else {
    gsdf <- data.frame("gene"=sgenes,
                       "when"=gswhen,
                       "when_cell"=colnames(hmat_bin_cellorder)[gswhen],
                       "switch_group"=switches_df_idx,
                       "switch_group_order"=1:length(sgenes))
}

    return(gsdf)
}

order_genes <- function(switches_df_idx) {
    if ((as.numeric(as.character(gene_switches_df_table[switches_df_idx,"on_switches"])) +
             as.numeric(as.character(gene_switches_df_table[switches_df_idx,"off_switches"]))) < 2) {

        if (as.logical(gene_switches_df_table[switches_df_idx,"start_bin"])) {
            whichswitch <- off_switches
        } else {
            whichswitch <- on_switches
        }
        
        return(order_by_switch(switches_df_idx, whichswitch))

    } else { return(order_by_cluster(switches_df_idx)) }}

find_change_genes <- function(idx) {
cat(idx)
    x <- gwhen[(gwhen$when >= gene_change_df_table[idx,"left"]) &
                   (gwhen$when <= gene_change_df_table[idx,"right"]),c("gene","when","when_cell")]

    x <- x[!is.na(x$gene),]
    x <- x[order(x$when),]
    x <- cbind(x, "change_window"=idx, "change_window_order"=1:nrow(x))
    
    return(x)
}

gene_change_df_table <- data.frame("change_group"=1,
                             "left"=c(1),
                             "right"=c(758),
                             stringsAsFactors=FALSE)


gwhen <- order_genes(1)
for (i in 2:nrow(gene_switches_df_table)) {
    gwhen <- rbind(gwhen, order_genes(i))
}

gwhen$gene <- as.character(gwhen$gene)
gwhen <- unique(gwhen)

change_genes <- find_change_genes(1)
## sapply(2:nrow(gene_change_df_table),
##        function(x)
##            change_genes <<- rbind(change_genes, find_change_genes(x)))


gwhen <- merge(gwhen, change_genes)
#gwhen <-gwhen[gwhen$switch_group == 2,]

get_hmap <- function(genes, binmat, whichorder, bulk=FALSE) {
#### vit_mat for _hmm, hmat_bin for _bin, and hmat for logtpm
#    genes <- unique(c(geneorder, connectors))
    if (any(!genes %in% rownames(binmat))) {
        genes <- gsub("[.]","-",genes)
    }

    hmap <- melt(binmat[genes,])
    colnames(hmap) <- c("gene","when_cell_all","value")

    gwhen_sub <- gwhen[gwhen$gene %in% genes,]

    if (whichorder == "switch_group") {
        gwhen_sub <- gwhen_sub[,c("gene","when","when_cell","switch_group","switch_group_order")]
    } else if (whichorder == "change_window") {
        gwhen_sub <- gwhen_sub[,c("gene","when","when_cell","change_window","change_window_order")]
    }

    colnames(gwhen_sub) <- c("gene","when","when_cell","gorder","gorder_order")

hmap_merge <- merge(hmap, gwhen_sub, all=T, by=c("gene"))
    if (!bulk) { 
hmap_merge$when_cell_all <- factor(hmap_merge$when_cell_all,
                                       levels=colnames(log_sc_time_tpm)[cellorder])
}


    hmap_res <- hmap_merge[!is.na(hmap_merge$gorder),]
    hmap_res$gene <- factor(hmap_res$gene,
                         levels=gwhen_sub$gene[order(gwhen_sub$gorder, gwhen_sub$gorder_order)])

    return(hmap_res)
}

plot_hmap <- function(myhmap, dofacet=TRUE) {
    p <- ggplot(myhmap, aes(y=gene, x=when_cell_all)) +
        geom_tile(aes(fill = value),
                  colour = ifelse(as.character(myhmap$when_cell_all) == as.character(myhmap$when_cell), "red", NA),
                  size=1) +
                      scale_fill_gradient(high=rev(brewer.pal(7,"YlGnBu"))[1],
                                          low="white",na.value="white")+
                                                             theme(axis.text.y = element_text(hjust = 1,
                                                                       size=12,
                                                                       face="bold"),
                                                                   plot.background=element_blank(),
                                                      axis.text.x=element_blank(),
                                                                   axis.ticks.x=element_blank(),
                                                                   legend.position="none") +
                                                                       xlab("Pseudotime Ordered Cells")+                                        scale_x_discrete(position = "top") + ylab("Genes")
    if (dofacet) {
        p +  facet_grid(
            scales="free",
            space="free",
            rows=vars(gorder)) }

return(p)
}


## x <- get_hmap(as.character(gwhen[,"gene"]),t(vit_mat), "switch_group")
## x2 <-  get_hmap(as.character(gwhen[,"gene"]),t(vit_mat), "change_window")


hmap <- get_hmap(unique(gwhen$gene), hmat[,cellorder], "switch_group")
hmap_bin <- get_hmap(unique(gwhen$gene), hmat_bin_cellorder, "switch_group")
hmap_hmm <- get_hmap(unique(gwhen$gene), t(vit_mat), "switch_group")
hmap_bulk <- get_hmap(unique(gwhen$gene), log(as.matrix(bulk_time_ave)+1,10), "switch_group", bulk=TRUE)


p <- plot_hmap(hmap)
p_bin <- plot_hmap(hmap_bin)
p_hmm <- plot_hmap(hmap_hmm)
p_bulk <- plot_hmap(hmap_bulk)

pdf("/mnt/gtklab01/ahjung/bivalent/figures/sc_time_switch.pdf",
width=20, height=50)

grid.arrange(p,p_bin,p_hmm,p_bulk,ncol=4)
          
dev.off()




################# DESEQ
library(DESeq2)
dds <- DESeqDataSetFromMatrix(round(bulk_time), bulk_time_coldata, design=~ idx + exp)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

res <- results(dds, name="exp_96h_vs_12h")
#res$padj < 0.05

dds_normalized <- counts(dds, normalized=TRUE)

bulk_time_ave <- cbind("12h"=apply(dds_normalized[,c(1,2,3)],1,mean),
"24h"=apply(dds_normalized[,c(4,5,6)],1,mean),
"36h"=apply(dds_normalized[,c(7,8,9)],1,mean),
"72h"=apply(dds_normalized[,c(10,11,12)],1,mean),
"96h"=apply(dds_normalized[,c(13,14,15)],1,mean))

################ network analysis

change_edge_dir <- function(edge) {
    fgene <- Rgraphviz::from(edge)
    tgene <- Rgraphviz::to(edge)
    if ((nrow(gwhen[gwhen$gene == fgene,])==0)|(nrow(gwhen[gwhen$gene == tgene,])==0)) {
        "none" } else if ((is.na(gwhen$when[gwhen$gene==fgene])) | (is.na(gwhen$when[gwhen$gene==tgene]))) {
            return("none")
        } else if ((gwhen$when[gwhen$gene==fgene] - gwhen$when[gwhen$gene==tgene])<0) {
            return("forward")
        } else if ((gwhen$when[gwhen$gene==fgene] - gwhen$when[gwhen$gene==tgene])>0) {
            return("back")
        } else {
            return("none")
    }}


check_edge_color <- function(fgene, tgene) {
    fcol <- ifelse(fgene %in% gwhen$gene[gwhen$switch_group==1],
                   "blue",
                   ifelse(fgene %in% gwhen$gene[gwhen$switch_group==2],
                                 "red",
                                 "green"))
    tcol <- ifelse(tgene %in% gwhen$gene[gwhen$switch_group==1],
                   "blue",
                   ifelse(tgene %in% gwhen$gene[gwhen$switch_group==2],
                          "red",
                          "green"))
    return(ifelse(fcol != tcol, fcol, "black"))
} 


change_edge_color <- function(edge) {
    fgene <- Rgraphviz::from(edge)
    tgene <- Rgraphviz::to(edge)
    if ((nrow(gwhen[gwhen$gene == fgene,])==0)|(nrow(gwhen[gwhen$gene == tgene,])==0)) {
        "grey" } else if ((is.na(gwhen$when[gwhen$gene==fgene])) | (is.na(gwhen$when[gwhen$gene==tgene]))) {
            return("grey")
        } else if (gwhen$when[gwhen$gene==fgene] - gwhen$when[gwhen$gene==tgene]<(-10)) {
            return(check_edge_color(fgene, tgene))
        } else if ((gwhen$when[gwhen$gene==fgene] - gwhen$when[gwhen$gene==tgene])>10) {
            return(check_edge_color(tgene, fgene))
        } else {
            return("grey")
        }}


change_edge_weight_switch <- function(edge) {
    fgene <- Rgraphviz::from(edge)
    tgene <- Rgraphviz::to(edge)
    if ((nrow(gwhen[gwhen$gene == fgene,])==0)|(nrow(gwhen[gwhen$gene == tgene,])==0)) {
        return(1) } else if ((is.na(gwhen$when[gwhen$gene==fgene])) | (is.na(gwhen$when[gwhen$gene==tgene]))) {
            return(10)
        } else if (gwhen$when[gwhen$gene==fgene] < gwhen$when[gwhen$gene==tgene]) {
            w <- (gwhen$when[gwhen$gene==tgene] - gwhen$when[gwhen$gene==fgene])
            return(w)
        } else if (gwhen$when[gwhen$gene==fgene] > gwhen$when[gwhen$gene==tgene]) {
            w <- (gwhen$when[gwhen$gene==fgene] - gwhen$when[gwhen$gene==tgene])
            return(w)
        } else {
            return(1)
        }}

change_edge_weight_mi <- function(edge) {
    fgene <- Rgraphviz::from(edge)
    tgene <- Rgraphviz::to(edge)
    return(mim[fgene,tgene])
}

mim_genes <- unique(c(downstream[downstream %in% rownames(mim_gene)], connectors))

mim <- mim_gene[mim_genes, mim_genes]
#mim <- mim_gene[unique(c(geneslist,"SP8","WNT5A","FOS","GATA3","COL1A2","GSC","WLS","GATA6")),
 #               unique(c(geneslist,"SP8","WNT5A","FOS","GATA3","COL1A2","GSC","WLS","GATA6"))]
net2_o <- aracne(mim)
mygraph_mim <- as(net2_o ,"graphNEL")

gdf <- data.frame(gene=gsub("[.]","-",gwhen$gene),
                  updown=ifelse(gwhen$switch_group==1,
                      "blue",
                      ifelse(gwhen$switch_group==2,
                             "red",
                             "green")),
                  stringsAsFactors=FALSE)



gvec <- as.vector(gdf$updown)
names(gvec) <- gdf$gene




gcolors <- brewer.pal(5,"Set3")
gcolor_df <- data.frame(name=c(mim_genes), col="white", stringsAsFactors=FALSE)
gcolor_df$col[gcolor_df$name %in% pathway_downstream_l[[1]]]<- gcolors[1]
gcolor_df$col[gcolor_df$name %in% pathway_downstream_l[[2]]]<- gcolors[2]
gcolor_df$col[gcolor_df$name %in% pathway_downstream_l[[3]]]<- gcolors[3]
gcolor_df$col[gcolor_df$name %in% pathway_downstream_l[[4]]]<- gcolors[4]



## gcolor_df$col[gcolor_df$name %in% nfat_geneslist] <- gcolors[2]
## gcolor_df$col[gcolor_df$name %in% matrisome] <- gcolors[3]
#gcolor_df$col[gcolor_df$name %in% colnames(mim)[as.numeric(sp5_all$res[[1]])]] <- gcolors[3]

nAttrs <- list()

nAttrs$fillcolor <- as.character(gcolor_df$col)
names(nAttrs$fillcolor) <- gcolor_df$name
nAttrs$fontsize <- rep(20, length(c(mim_genes)))
names(nAttrs$fontsize) <- c(mim_genes)
nAttrs$size <- rep(1, length(c(mim_genes)))
names(nAttrs$size) <- c(mim_genes)
nAttrs$color <- gvec

eAttrs <- list()

edges <- buildEdgeList(mygraph_mim)

testvec <- c(1:length(edges))
names(testvec) <- names(edges)
          
edge_dir <- data.frame("edgename" = names(edges),
                       "edgedir" = sapply(1:length(edges),
                           function(x)
                               change_edge_dir(edges[[x]])),
                       "edgesize" = order(sapply(1:length(edges),
                           function(x)
                               change_edge_weight_switch(edges[[x]]))),
                       "edgecolor" =  sapply(1:length(edges),
                           function(x)
                               change_edge_color(edges[[x]])),
                       "arrowhead" = "box",
                       stringsAsFactors=FALSE
                       )
#          edge_dir$edgecolor[edge_dir$edgedir == "none"] <- "grey"
edge_dir$arrowhead[edge_dir$edgedir == "none"] <- "open"

eAttrs$dir <- as.character(edge_dir$edgedir)
names(eAttrs$dir) <- edge_dir$edgename

eAttrs$color <- as.character(edge_dir$edgecolor)
names(eAttrs$color) <- edge_dir$edgename

eAttrs$arrowhead <- as.character(edge_dir$arrowhead)
names(eAttrs$arrowhead) <- edge_dir$edgename

eAttrs$weight <- edge_dir$edgesize
names(eAttrs$weight) <- edge_dir$edgename

## eAttrs$label <- 684-edge_dir$edgesize
## names(eAttrs$label) <- edge_dir$edgename

          
pdf("/mnt/gtklab01/ahjung/bivalent/figures/final_graph.pdf",
width=20, height=10)

Rgraphviz::plot(mygraph_mim, nodeAttrs = nAttrs, edgeAttrs = eAttrs)
          
dev.off()



## ############# find genes with similar switching time
## gwhen_sub <- gwhen[!gwhen$gene %in% downgenes,]

## genesA <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 90) & (gwhen_sub$when <= 105)])
## genesB <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 106) & (gwhen_sub$when <= 179)])
## genesC <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 180) & (gwhen_sub$when <= 200)])
## genesD <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 201) & (gwhen_sub$when <= 239)])
## genesE <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 240) & (gwhen_sub$when <= 280)])
## genesF <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 271) & (gwhen_sub$when <= 379)])
## genesG <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 380) & (gwhen_sub$when <= 440)])


## ## plotgenes <- c("T","EOMES","BAMBI","NOG","WNT3","WNT5B","FOS","DKK1","PITX2","SP5","LEF1","SP8","GSC","COL1A2","SLC24A4","WNT5A","NFATC1","PPP3CC","SMAD2","GATA3","POU5F1","NODAL","FZD8","FZD5","BMP2","BMP4","GPR83")

## plotgenes_bin <- plotgenes[plotgenes %in% colnames(vit_mat)]
## plotgenes_nobin <- plotgenes[!plotgenes %in% colnames(vit_mat)]

## plotgenes_bin_hmap <- get_hmap(plotgenes_bin, t(vit_mat))
## plotgenes_nobin_hmap <- get_hmap(plotgenes_nobin, log_sc_time_tpm[,cellorder])

## q_bin <- plot_hmap(plotgenes_bin_hmap)
## q_nobin <- plot_hmap(plotgenes_nobin_hmap)



## df <- data.frame("cells"=colnames(log_sc_time_tpm),
##            "value"=rep(1, ncol(log_sc_time_tpm)),
##            "cluster"=cl1,
##            "timepoint"=sc_time_coldata$exp)

## df_ordered <- df[cellorder,]

## q_cells <- ggplot(df_ordered, aes(x=1:758, y=value)) +
##         geom_tile(aes(fill = timepoint)) +
##             scale_fill_manual(values=cols) + theme_bw() +
##                 theme(legend.position="bottom",
##                       plot.margin=margin(l=0.3, r=-0.3,unit="cm"))

                      

## grid.arrange(q_bin, q_nobin, q_cells, ncol=1)


##################################### gsa analysis
gsa_res <- read.table("/mnt/gtklab01/ahjung/bivalent/results/GSA/sc_time_bik27_notk4_KEGG_q0.001_all.csv",
                      fill=TRUE, sep="\t" ,header=FALSE,skip=9,nrows=51)
gsa_res <- gsa_res[-1,1:7]
colnames(gsa_res) <- c("geneset_name","genes_in_geneset","description","genes_in_overlap","k_K","pvalue","FDRqvalue")

gsa_mat <- read.table("/mnt/gtklab01/ahjung/bivalent/results/GSA/sc_time_bik27_notk4_KEGG_q0.001_all.csv",
                      fill=TRUE, sep="\t" ,header=TRUE,skip=64,nrows=1200, stringsAsFactors=FALSE, quote="")

#pathways <- colnames(gsa_mat)[1:32][grep("SIGNALING", colnames(gsa_mat)[1:32])]
pathways <- colnames(gsa_mat)[unique(c(grep("PATHWAY", colnames(gsa_mat)),
                                       grep("SIGNALING", colnames(gsa_mat)),
                                       grep("PATHWAY",colnames(gsa_mat)),
                                       grep("ADHESION", colnames(gsa_mat)),
                                       grep("MATRISOME", colnames(gsa_mat))))]


gsa_mat_tf <- gsa_mat[,pathways] != ""
rownames(gsa_mat_tf) <- gsa_mat$Gene.Symbol

cols_pathways <- brewer.pal(length(pathways), "Set3")

gsa_mat_m <- melt(gsa_mat_tf[apply(gsa_mat_tf,1,sum)>0,])

gsa_mat_m_overlap <- merge(gsa_mat_m,
                           data.frame("Var1"=names(apply(gsa_mat_tf,1,sum)),
                                      "overlap"=as.numeric(apply(gsa_mat_tf,1,sum))))
colnames(gsa_mat_m_overlap) <- c("gene","pathway","value","overlap")


x2 <-  get_hmap(as.character(gwhen[,"gene"]),t(vit_mat), "switch_group")

gsa_hmap <- x2[x2$gene %in% colnames(vit_mat)[colnames(vit_mat) %in% as.character(unique(gsa_mat_m$Var1))],]


gsa_mat_m_overlap <- merge(gsa_mat_m_overlap, gwhen[,c("gene","switch_group")],all=T)

gsa_mat_m_overlap$gene <- factor(gsa_mat_m_overlap$gene,
                                 levels = levels(gsa_hmap$gene)[levels(gsa_hmap$gene) %in% levels(gsa_mat_m_overlap$gene)])

gsa_mat_m_overlap <- gsa_mat_m_overlap[!is.na(gsa_mat_m_overlap$gene),]
gsa_mat_m_overlap <- gsa_mat_m_overlap[!is.na(gsa_mat_m_overlap$pathway),]
gsa_mat_m_overlap <- gsa_mat_m_overlap[!is.na(gsa_mat_m_overlap$switch_group),]

gsa_overlap <- ggplot(gsa_mat_m_overlap[gsa_mat_m_overlap$pathway %in% pathways[c(3,5,8,2,9,13,21,23,26,27)],],#[gsa_mat_m_overlap$gene %in% gsa_mat_motif_m_overlap$gene,],
                      aes(y=gene, x=pathway)) +
    geom_tile(aes(fill = ifelse(value, overlap, 0))) + 
        scale_fill_gradient(low = "white", high = "black") +
            scale_x_discrete(
               labels=c(
                   "TGFB",
                   "SMAD2_3N",
                   "P53",
                   "P53_DOWNSTREAM",
                   "MAPK",
                   "WNT",
                   "WNT",
                   "CYTOKINE",
                   "FOCAL_ADHESION",
                   "MATRISOME"
                                      ),
                             position="top") +
    theme_bw() +
        theme(axis.text.y = element_text(hjust = 0.5,
                  size=10,
                  face="bold"),
              axis.text.x = element_text(hjust = 0, angle = 45, size=12),
              legend.position="none",
               plot.margin = margin(0, 1, 0.2, 0, "cm")) +
                      ## facet_grid(
                      ##     scales="free",
                      ##     space="free",
                      ##     rows=vars(change_window)
                      ##     )+  
                    ylab("") + xlab("")

######### MKNK1 and SOCS7 expression could not be binarized



pdf("/mnt/gtklab01/ahjung/bivalent/figures/gsa_heatmap_overlap.pdf",
width=10, height=20)

l <- plot_hmap(gsa_hmap,dofacet=FALSE)#[gsa_hmap$gene %in% gsa_mat_motif_m_overlap$gene,])
l <- l + theme( plot.margin = margin(3, 0, 0.2, 0, "cm"))
grid.arrange(l, gsa_overlap, ncol=2)

dev.off()

############## GSA motif
gsa_res_motif <- read.table("/mnt/gtklab01/ahjung/bivalent/results/GSA/sc_time_bik27_notk4_TFmotif_q0.001.csv",
                      fill=TRUE, sep="\t" ,header=FALSE,skip=9,nrows=101)
gsa_res_motif <- gsa_res_motif[-1,1:7]
colnames(gsa_res_motif) <- c("geneset_name","genes_in_geneset","description","genes_in_overlap","k_K","pvalue","FDRqvalue")

gsa_mat_motif <- read.table("/mnt/gtklab01/ahjung/bivalent/results/GSA/sc_time_bik27_notk4_TFmotif_q0.001.csv",
                      fill=TRUE, sep="\t" ,header=TRUE,skip=114,nrows=1200, stringsAsFactors=FALSE, quote="")

motifs <- colnames(gsa_mat_motif)[unique(c(grep("MAZ", colnames(gsa_mat_motif))[1],
                                           grep("SP1", colnames(gsa_mat_motif))[1],
                                           grep("LEF1", colnames(gsa_mat_motif))[1],
                                           grep("E12", colnames(gsa_mat_motif))[1],
                                           grep("AP1", colnames(gsa_mat_motif))[1],
                                           grep("NFAT", colnames(gsa_mat_motif))[1],
                                           grep("NFKB", colnames(gsa_mat_motif))[1]
    ##                                        grep("SMAD4", colnames(gsa_mat_motif)),
    ##                                        grep("FOXO1", colnames(gsa_mat_motif)),
    ##                                        grep("IRF1", colnames(gsa_mat_motif)),
    ##                           
#                                           grep("FOXO3", colnames(gsa_mat_motif)),
    ##                                        grep("ETS2", colnames(gsa_mat_motif)),
    ##                                        grep("SRF", colnames(gsa_mat_motif))[1],
    ## grep("PITX2", colnames(gsa_mat_motif)))
    ))]


gsa_hmap_motif <- x2[x2$gene %in% c(colnames(vit_mat)[colnames(vit_mat) %in%
                                                      as.character(unique(gsa_mat_m_motif$Var1))]),]

gsa_mat_motif_tf <- gsa_mat_motif[,motifs] != ""
rownames(gsa_mat_motif_tf) <- gsa_mat_motif$Gene.Symbol

#cols_motifs <- brewer.pal(length(motifs), "Set3")

gsa_mat_motif_m <- melt(gsa_mat_motif_tf[apply(gsa_mat_motif_tf,1,sum)>0,])

gsa_mat_motif_m_overlap <- merge(gsa_mat_motif_m,
                           data.frame("Var1"=names(apply(gsa_mat_motif_tf,1,sum)),
                                      "overlap"=as.numeric(apply(gsa_mat_motif_tf,1,sum))))
colnames(gsa_mat_motif_m_overlap) <- c("gene","motif","value","overlap")

gsa_mat_motif_m_overlap <- merge(gsa_mat_motif_m_overlap, gwhen[,c("gene","switch_group")],all=T)

gsa_mat_motif_m_overlap$gene <- factor(gsa_mat_motif_m_overlap$gene,
                                 levels = levels(x2$gene)[levels(x2$gene) %in% levels(gsa_mat_motif_m_overlap$gene)])

gsa_mat_motif_m_overlap <- gsa_mat_motif_m_overlap[!is.na(gsa_mat_motif_m_overlap$switch_group),]
gsa_mat_motif_m_overlap <- gsa_mat_motif_m_overlap[!is.na(gsa_mat_motif_m_overlap$motif),]


## gsa_hmap_motif$gene <- factor(gsa_hmap_motif$gene,
##                               levels = levels(gsa_mat_motif_m_overlap$gene)[levels(gsa_mat_motif_m_overlap$gene) %in% levels(x2$gene)])




## gsa_mat_tf_motif <- gsa_mat_motif[,motifs] != ""
## rownames(gsa_mat_tf_motif) <- gsa_mat_motif$Gene.Symbol

## #cols_pathways <- brewer.pal(length(pathways), "Set3")

## gsa_mat_m_motif <- melt(gsa_mat_tf_motif[apply(gsa_mat_tf_motif,1,sum)>0,])

## gsa_mat_m_overlap_motif <- merge(gsa_mat_m_motif,
##                            data.frame("Var1"=names(apply(gsa_mat_tf_motif,1,sum)),
##                                       "overlap"=as.numeric(apply(gsa_mat_tf_motif,1,sum))))
## colnames(gsa_mat_m_overlap_motif) <- c("gene","motif","value","overlap")

## gsa_mat_m_overlap_motif <- merge(gsa_mat_m_overlap_motif, gwhen[,c("gene","gorder")],all=T)

## gsa_mat_m_overlap_motif$gene <- factor(gsa_mat_m_overlap_motif$gene,
##                                  levels = levels(gsa_hmap_motif$gene)[levels(gsa_hmap_motif$gene) %in% levels(gsa_mat_m_overlap_motif$gene)])

## gsa_mat_m_overlap_motif <- gsa_mat_m_overlap_motif[!is.na(gsa_mat_m_overlap_motif$gene),]
## gsa_mat_m_overlap_motif <- gsa_mat_m_overlap_motif[!is.na(gsa_mat_m_overlap_motif$motif),]
## gsa_mat_m_overlap_motif <- gsa_mat_m_overlap_motif[!is.na(gsa_mat_m_overlap_motif$switch_group),]

gsa_overlap_motif <- ggplot(gsa_mat_motif_m_overlap,
    #gsa_mat_motif_m_overlap[gsa_mat_motif_m_overlap$gene %in% gsa_mat_m_overlap$gene,],
                      aes(y=gene, x=motif)) +
    geom_tile(aes(fill = ifelse(value, overlap, 0))) + 
        scale_fill_gradient(low = "white", high = "black") +
            scale_x_discrete(                             position="top") +
    theme_bw() +
        theme(axis.text.y = element_text(hjust = 0.5,
                  size=6,
                  face="bold"),
              axis.text.x = element_text(hjust = 0, angle = 45, size=10),
              legend.position="none")+
#               plot.margin = margin(0, 1, 1, 0, "cm")) +
                      ## facet_grid(
                      ##     scales="free",
                      ##     space="free",
                      ##     rows=vars(change_window)) +
                    ylab("") + xlab("")

lm <- plot_hmap(gsa_hmap_motif)
lm <- lm + theme( plot.margin = margin(3, 0, 0.2, 0, "cm"))

pdf("/mnt/gtklab01/ahjung/bivalent/figures/gsa_overlap_motif.pdf",
width=10, height=20)

grid.arrange(lm, gsa_overlap_motif, ncol=2)

dev.off()



################################################



edges <- buildEdgeList(mygraph_mim)
edge_df <- data.frame("from"=unlist(lapply(edges, Rgraphviz::from)),
                      "to"=unlist(lapply(edges, Rgraphviz::to)))

get_neigh_nodes <- function(gene) {
return(as.character(unique(unlist(c(edge_df[(edge_df$from==gene)|(edge_df$to==gene),])))))
}

get_neigh_nodes("LHX1") %in% 

gwhen[gwhen$gene %in% get_neigh_nodes("LHX1"),"when"]



####################################

sgenes <- as.character(unique(gsa_mat_m_overlap$gene))

on_path <- as.data.frame(matrix(rep(0, length(on_switches)*758), ncol=758))
rownames(on_path) <- names(on_switches)

for (x in sgenes) {
    on_path[x,on_switches[[x]]] <- 1 }

for (x in sgenes) {
    on_path[x,off_switches[[x]]] <- -1 }




get_switch_step <- function(htable, pathway) {
    on_path_sub <- on_path[colnames(vit_mat)[colnames(vit_mat) %in% names(which(htable[,pathway]))],]
    on_path_sub <- on_path_sub[!is.na(on_path_sub[,1]),]
    gene <- names(which(sapply(strsplit(pathway, "_")[[1]], function(x) x %in% rownames(on_path_sub))))
    if (length(rownames(on_path_sub) == gene)!=0) {
        on_path_sub <- on_path_sub[!rownames(on_path_sub)==gene,] }
    step_path <- cumsum(as.numeric(apply(on_path_sub,2,sum)))+
        sum((vit_mat_path[1,rownames(on_path_sub)]=="ON")*1)
    return(step_path)
}

df_step_motif <- data.frame("step_path"=matrix(sapply(motifs,
                                function(x) get_switch_step(gsa_mat_motif_tf, x)),
                                ncol=1),
                            "cells"=rep(1:758,length(motifs)),
                            "class"=matrix(sapply(motifs, function(x) rep(x, 758)), ncol=1),
                            "type"="motif")

df_step <- data.frame("step_path"=matrix(sapply(pathways, function(x) get_switch_step(gsa_mat_tf,x)), ncol=1),
                      "cells"=rep(1:758,length(pathways)),
                      "class"=matrix(sapply(pathways, function(x) rep(x, 758)), ncol=1),
                      "type"="pathway")

step_pathways <- ggplot(df_step[df_step$class %in% pathways[c(23,26,27)],],## c(3,5,8,2,21,9,13) rbind(df_step_motif[df_step_motif$class%in%motifs[c(4,5,6)],],
                     ## df_step[df_step$class %in% pathways[c(1,2,3,6)],]),
                     aes(x=cells, y=step_path, col=class)) +
                         geom_line(size=1.5) +
                      facet_grid(
                          scales="free",
                          space="free",
                          cols=vars(class))+
        theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.position = "right")

step_motifs <- ggplot(df_step_motif[df_step_motif$class %in% motifs,],## c(3,5,8,2,21,9,13) rbind(df_step_motif[df_step_motif$class%in%motifs[c(4,5,6)],],
                     ## df_step[df_step$class %in% pathways[c(1,2,3,6)],]),
                     aes(x=cells, y=step_path, col=class)) +
                         geom_line(size=1.5) +
                      facet_grid(
                          scales="free",
                          space="free",
                          cols=vars(class))+
        theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.position = "right")


grid.arrange(step_1,step_2, ncol=2)



##################### fractional ON
library(scales)





plot_fraction <- function(gene) {

d <- as.numeric(log_sc_time_tpm_messup_bin[gene,cellorder])
d <- unlist(lapply(split(d, ceiling(seq_along(d)/20)),
       function(x)
           sum(x)/length(x)))

p <- ggplot(data.frame("fractionON"=d,
                  "cwindow"=1:length(d)), aes(x=cwindow, y=fractionON, col=as.factor(cwindow))) +
                      geom_point() + geom_smooth(col="grey", se=FALSE, size=1) +
                          theme_classic() + ggtitle(gene) + scale_color_grey(start=0.8, end=0.2) + theme(              legend.position="none") + ylab("fraction of active cells") + xlab("pseudotime windows")  
return(p)
}

plot_logtpm <- function(gene) {

d <- as.numeric(log_sc_time_tpm[gene,cellorder])
p <- ggplot(data.frame("logTPM"=d,
                  "cells"=1:length(d)), aes(x=cells, y=logTPM)) +
                      geom_point(alpha=0.2) + 
                          theme_classic() + ggtitle(gene) + scale_color_grey(start=0.8, end=0.2) + theme(              legend.position="none") + ylab("log(TPM+1)") + xlab("cells in pseudotime order") 
return(p)
}


plot_twindow <- function(gene) {

d1 <- as.numeric(log_sc_time_tpm[gene,cellorder])
df <- data.frame("logTPM"=unlist(split(d1, ceiling(seq_along(d1)/130))),
                 "cwindow"=rep(1:6, lapply(split(d1, ceiling(seq_along(d1)/130)), length)),
                 "cells"=unlist(lapply(split(d1, ceiling(seq_along(d1)/130)), function(x) 1:length(x))))

q <- ggplot(df, aes(logTPM, col=as.factor(cwindow), fill=as.factor(cwindow))) +
    geom_density(size=1) + facet_grid(cols=vars(cwindow)) + xlim(-2,10) + theme_classic() +
        theme(strip.background = element_blank(),
              strip.text.y = element_blank(),
              legend.position="none") + scale_fill_grey(start=0.8, end=0.2)+ xlab("log(TPM+1)") + scale_color_grey(start=0.8, end=0.2) + coord_flip() 

return(q)

}

## plot_thmm <- function(gene) {
## d1 <- vit_mat_path[,gene]

## df <- data.frame("hmm_path"=unlist(split(d1, ceiling(seq_along(d1)/130))),
##                  "cwindow"=rep(1:6, lapply(split(d1, ceiling(seq_along(d1)/130)), length)),
##                  "cells"=unlist(lapply(split(d1, ceiling(seq_along(d1)/130)), function(x) 1:length(x))))

## q <- ggplot(df, aes(y=cells, x=hmm_path, col=hmm_path, fill=hmm_path)) +
##     geom_bar(position="dodge",stat="identity") + theme_classic() +
##         theme(strip.background = element_blank(),
##               strip.text.y = element_blank(),
##               legend.position="none") + facet_grid(                          scales="free",
##                           space="free",cols=vars(cwindow)) + scale_fill_grey(start=0.8, end=0.2)+  scale_color_grey(start=0.8, end=0.2)  

## return(q)

## }



pdf("/mnt/gtklab01/ahjung/bivalent/figures/sc_time_example.pdf",width=10, height=6)

grid.arrange(plot_fraction("PITX2"),
#             plot_fraction("LPGAT1"),
             plot_fraction("EOMES"),
             plot_twindow("PITX2"),
#             plot_twindow("LPGAT1"),
             plot_twindow("EOMES"),
             plot_fraction("PMAIP1"),
#             plot_fraction("KLF10"),
             plot_fraction("FOXH1"),
             plot_twindow("PMAIP1"),
#             plot_twindow("KLF10"),
             plot_twindow("FOXH1"),
             ncol=2)

dev.off()

#lapply(gene_switches, function(x) c("PITX2","LPGAT1","EOMES","PMAIP1","KLF10","FOXH1") %in% x)


############################3
all_motif <- as.character(unlist(sapply(colnames(gsa_mat_motif), function(x) strsplit(x, "_"))))
all_motif <- all_motif[-grep("^Q",all_motif)]
all_motif <- all_motif[-grep("^[0-9]",all_motif)]
all_motif <- all_motif[-grep("CC",all_motif)]
all_motif <- all_motif[-grep("AA",all_motif)]
all_motif <- all_motif[-grep("GG",all_motif)]
all_motif <- all_motif[-grep("TT",all_motif)]
all_motif <- all_motif[-grep("CT",all_motif)]
all_motif <- all_motif[-grep("GT",all_motif)]
all_motif <- all_motif[-grep("TG",all_motif)]
all_motif <- all_motif[-grep("UNKNOWN",all_motif)]
all_motif <- all_motif[-grep("Gene",all_motif)]

motif_genes <- rbind(data.frame("original"=unique(all_motif[all_motif %in% rownames(log_sc_time_tpm)]),
           "alias"=NA),
      data.frame("original"=unique(all_motif[!all_motif %in% rownames(log_sc_time_tpm)])[c(-8,-9,-22)],
                 "alias"=c("ELSPBP1", NA, "NFATC1", "TFAP4",NA, "FOXF2","JUN","VSX2","ESRRA",NA,"NFKB1","MYOD1", "GTF2A1", "ZEB1", "PDX1", NA, "NFYA","IRF1", "TFCP2","EBF1","MAMSTER", "CREBZF","GCM1","FOXA1","FOXA2","NFKB1","THRB", "GABPB1", "ZIC3", "NR1I2","PAXIP1","LMO2", NA, "NHLH1", "GTF3A", "FOXC1", "FOXJ1", "NFKB1", "GATA", "ELK1", "NFIL3", NA, "EREG", NA, "PATZ1")) )

bimo_motif <- unique(c(as.character(motif_genes$original[motif_genes$original %in% names(bimocondition)[bimocondition]]),
as.character(motif_genes$alias[motif_genes$alias %in% names(bimocondition)[bimocondition]])))



plot_motif_targets <- function(gene) {
    if (gene %in% motif_genes$alias) { gene <- as.character(motif_genes$original[as.character(motif_genes$alias) %in% gene]) }
    x <- data.frame("step_path"=df_step_motif[grep(gene, df_step_motif$class),"step_path"][1:758],
                    "cells"=1:758)
    noverlap <- sum(gsa_mat_motif[,as.character(unique(df_step_motif[grep(gene, df_step_motif$class)[1],"class"]))] != "")
    p <- ggplot(x,  aes(x=cells, y=step_path)) +
        geom_line(size=1) +
            theme_classic() + ggtitle(paste0("Genes with motif for ",gene, " (",noverlap,")"))+ ylab("number of active genes") + xlab("cells in psudotime order") +  scale_y_continuous(breaks= pretty_breaks())
    return(p)
}

pdf("/mnt/gtklab01/ahjung/bivalent/figures/motif_targets.pdf",   width=4, height=15)

grid.arrange(plot_logtpm("NFATC1"),
             plot_logtpm("NFATC2"),
             plot_logtpm("NFATC3"),
             plot_logtpm("NFATC4"),
             plot_logtpm("NFAT5"),
             plot_motif_targets("NFAT"),ncol=1)

dev.off()

pdf("/mnt/gtklab01/ahjung/bivalent/figures/motif_targets_2.pdf",
   width=14, height=6)

grid.arrange(plot_logtpm("MAZ"),
             plot_logtpm("SP1"),
             plot_logtpm("LEF1"),
             plot_logtpm("MYC"),
             plot_motif_targets("MAZ"),
             plot_motif_targets("SP1"),
             plot_motif_targets("LEF1"), 
            plot_motif_targets("MYC"),
             ncol=4)

dev.off()

grid.arrange(plot_logtpm("SOX9"),
             plot_motif_targets("SOX9"), ncol=1)



bimo_motif_select <- c("ETS2","IRF1","PITX2","FOXO3","FOXO1","TFAP4","FOXF2","JUN","ESRRA","GABPB1","GCM1","GTF3A","EREG","PATZ1")

bimo_motif_select <- c("IRF1","PITX2","FOXO1","FOXO3","GCM1")
bimo_motif_select_names <- motifs[sapply(as.character(motif_genes$original)[(as.character(motif_genes[,1]) %in% bimo_motif) | (as.character(motif_genes[,2]) %in% bimo_motif)], function(x) grep(x, motifs)[1])]

bimo_motif_select_names <- motifs[c(101)]

pdf("/mnt/gtklab01/ahjung/bivalent/figures/bimo_motif_targets.pdf",
    width=17, height=9)

grid.arrange(plot_logtpm(bimo_motif_select[1]),
             plot_logtpm(bimo_motif_select[2]),
             plot_logtpm(bimo_motif_select[3]),
             plot_logtpm(bimo_motif_select[4]),
             plot_logtpm(bimo_motif_select[5]),
             plot_fraction(bimo_motif_select[1]),
             plot_fraction(bimo_motif_select[2]),
             plot_fraction(bimo_motif_select[3]),
             plot_fraction(bimo_motif_select[4]),
             plot_fraction(bimo_motif_select[5]),
             plot_motif_targets(bimo_motif_select[1]),
             plot_motif_targets(bimo_motif_select[2]),
             plot_motif_targets(bimo_motif_select[3]),
             plot_motif_targets(bimo_motif_select[4]),
             plot_motif_targets(bimo_motif_select[5]),
             ncol=5)

dev.off()

downstream <- unique(unlist(apply(gsa_mat_motif[,bimo_motif_select_names],2, function(x) gsa_mat_motif[,2][x!=""])))

motif_downstream <- gsa_mat_motif[,2][gsa_mat_motif[,bimo_motif_select_names]!=""]

downstream_bin <- log_sc_time_tpm_messup_bin[downstream,cellorder][!is.na(log_sc_time_tpm_messup_bin[downstream,cellorder][,1]),]*1

## downstream_hmap <- get_hmap(unique(rownames(downstream_bin)),
##                             log_sc_time_tpm_messup_bin, "change_window")


pathway_downstream_l <-apply(gsa_mat_tf[,pathways[c(3,5,9,13)]],2,function(x) rownames(gsa_mat_tf)[x])
pathway_downstream <- unique(unlist(apply(gsa_mat_tf[,pathways[c(3,5,9,13)]],2,function(x) rownames(gsa_mat_tf)[x])))

downstream <- unique(c(pathway_downstream))
