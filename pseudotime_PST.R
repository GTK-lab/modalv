setwd("~/projects/bivalent/modalv")

print("pseudotime script\n\n")

source("GSE75748_data.R")
source("GSE75748_function.R")
source("pseudotime_function.R")

## sc_time_original <- sc_time
## sc_time <- sc_time[!apply(sc_time,1,sum)==0,]

## sc_time_tpm <- get_tpm(sc_time) # some genes get lost
## log_sc_time_tpm <- log(sc_time_tpm+1)

## log_sc_time_tpm_messup <- messup(log_sc_time_tpm, 1e-05)
## log_sc_time_tpm_messup_bimodal <- apply(log_sc_time_tpm_messup, 1, test_bimodal)

## set.seed(100)
## bimofit <- fit_bimodal_multi(log_sc_time_tpm_messup,ncores=30)

## bimocondition <- filter_condition(log_sc_time_tpm_messup_bimodal,
##                                   bimofit,
##                                   cond_dip=0.01,
##                                   cond_lambda=1,
##                                   cond_mu=1)
## bimocondition[is.na(bimocondition)] <- FALSE

## print("bimodality fitting done\n\n\n")


mclass <- sc_time_coldata$exp
nclass <- table(mclass)
cols <- get_colors(nclass)


## bimocondition_bi <- (rownames(log_sc_time_tpm) %in% bigenes) & bimocondition
## bimocondition_k4 <- (rownames(log_sc_time_tpm) %in% k4genes) & bimocondition
## bimocondition_k27 <- (rownames(log_sc_time_tpm) %in% k27genes) & bimocondition

## bimocondition_bi <- bimocondition_bi[!is.na(bimocondition_bi)]
## bimocondition_k4 <- bimocondition_k4[!is.na(bimocondition_k4)]
## bimocondition_k27 <- bimocondition_k27[!is.na(bimocondition_k27)]

## save(sc_time,
##      sc_time_tpm,
##      log_sc_time_tpm,
##      log_sc_time_tpm_messup,
##      log_sc_time_tpm_messup_bimodal,
##      bimofit,
##      bimocondition,
##      bimocondition_bi,
##      bimocondition_k4,
##      bimocondition_k27,
##      file = "files/sc_time_basic.RData")

load("files/sc_time_basic.RData")

library("Rtsne")
library("ggplot2")

# binarize

## log_sc_time_tpm_messup_bin <- binarizeexp(log_sc_time_tpm_messup[bimocondition,][1,])

## invisible(sapply(2:sum(bimocondition),
##                  function(i) log_sc_time_tpm_messup_bin <<- rbind(log_sc_time_tpm_messup_bin,
## binarizeexp(log_sc_time_tpm_messup[bimocondition,][i,]))))

## rownames(log_sc_time_tpm_messup_bin) <- names(bimocondition)[bimocondition]


## save(log_sc_time_tpm_messup_bin,
##      file="files/log_sc_time_tpm_messup_bin.RData")

load("files/log_sc_time_tpm_messup_bin.RData")

print("binarization done \n\n")

## binmat <- log_sc_time_tpm_messup_bin[!is.na(log_sc_time_tpm_messup_bin[,1]),]

## set.seed(100)
## tsne_out_bi <- Rtsne(t(unique(binmat[rownames(binmat) %in%names(bimocondition_bi)[bimocondition_bi],])),
##                   pca=FALSE,
##                   perplexity=30,
##                   theta=0.0)

## tsne1_bi <- ggplot(data=data.frame("C1"=tsne_out_bi$Y[,1],
##                     "C2"=tsne_out_bi$Y[,2],
##            "type"=mclass),
##        aes(C1,C2,col=type)) +
##     geom_point() + theme_bw() + scale_color_manual(values=cols) + ggtitle("tSNE based on bivalent bimodal genes")

## rd1 <- tsne_out_bi$Y[,1:2]
## rownames(rd1) <- mclass

## library(mclust)

## cl1 <- Mclust(rd1,6)$classification
## nclass_cluster <- table(cl1)

## tsne2 <- ggplot(data=data.frame("C1"=tsne_out_bi$Y[,1],
##                     "C2"=tsne_out_bi$Y[,2],
##            "type"=as.factor(cl1)),
##        aes(C1,C2,col=type)) +
##     geom_point() + theme_bw() + scale_color_manual(values=cols) + ggtitle("6 clusters generated")

## save(tsne_out_bi,
##      tsne1_bi,
##      rd1,
##      cl1,
##      nclass_cluster,
##      tsne2,
##      file="files/sc_time_tsne.RData")

load("files/sc_time_tsne.RData")

library(gridExtra)

## pdf("/mnt/gtklab01/ahjung/bivalent/figures/sc_time_tsne.pdf",
##     width=5, height=8)
## grid.arrange(tsne1_bi, tsne2)
## dev.off()

print("tsne done\n\n\n")


## hmmcondition <- (bimocondition_k27 | bimocondition_bi) & !bimocondition_k4

## hmat <- matrix(as.numeric(log_sc_time_tpm[hmmcondition,]),
##                nrow=sum(hmmcondition))
## rownames(hmat) <- rownames(log_sc_time_tpm[hmmcondition,])
## colnames(hmat) <- colnames(log_sc_time_tpm)


## hmat_bin <- log_sc_time_tpm_messup_bin[names(hmmcondition)[hmmcondition],]*1
## hmat_bin <- hmat_bin[!(is.na(hmat_bin[,1])),]


## save(hmmcondition,
##      hmat,
##      hmat_bin,
##      file="files/hmmcondition_matrix.RData")

load("files/hmmcondition_matrix.RData")

# mutual information matrix
## mim_cell <- sapply(1:ncol(hmat_bin),
##                    function(x) get_mi(t(hmat_bin), x))
## colnames(mim_cell) <- colnames(hmat_bin)
## rownames(mim_cell) <- colnames(hmat_bin)

## # conditional entropy matrix as distance between cells
## mim_cell_var <- sapply(1:ncol(hmat_bin),
##                    function(x) get_var(t(hmat_bin), x))
## colnames(mim_cell_var) <- colnames(hmat_bin)
## rownames(mim_cell_var) <- colnames(hmat_bin)

## mim_gene <- sapply(1:nrow(hmat_bin),
##                    function(x) get_mi(hmat_bin, x))
## colnames(mim_gene) <- rownames(hmat_bin)
## rownames(mim_gene) <- rownames(hmat_bin)

## # conditional entropy matrix as distance between genes
## mim_gene_var <- sapply(1:nrow(hmat_bin),
##                    function(x) get_var(hmat_bin, x))
## colnames(mim_gene_var) <- rownames(hmat_bin)
## rownames(mim_gene_var) <- rownames(hmat_bin)


## save(mim_cell,
##      mim_cell_var,
##      mim_gene,
##      mim_gene_var,
##      file="files/mim.RData")

load("files/mim.RData")

print("mutual information matrix done \n\n")
################### pseudotime ordering based on PST


## from1to2 <- get_connectingpair(1,2)
## from2to3 <- get_connectingpair(2,3)
## from3to4 <- get_connectingpair(3,4)
## from4to5 <- get_connectingpair(4,5)
## from5to6 <- get_connectingpair(5,6)

## cluster1_order <- t(sapply(from1to2$fromcell, function(x) rev(get_cellorder(x,1))))
## cluster2_order_from1 <- t(sapply(from1to2$tocell, function(x) get_cellorder(x,2)))
## cluster2_order_from3 <- t(sapply(from2to3$fromcell, function(x) rev(get_cellorder(x,2))))
## cluster3_order_from2 <- t(sapply(from2to3$tocell, function(x) get_cellorder(x,3)))
## cluster3_order_from4 <- t(sapply(from3to4$fromcell, function(x) rev(get_cellorder(x,3))))
## cluster4_order_from3 <- t(sapply(from3to4$tocell, function(x) get_cellorder(x,4)))
## cluster4_order_from5 <- t(sapply(from4to5$fromcell, function(x) rev(get_cellorder(x,4))))
## cluster5_order <- t(sapply(from4to5$tocell, function(x) get_cellorder(x,c(5,6))))

## save(from1to2,
##      from2to3,
##      from3to4,
##      from4to5,
##      from5to6,
##      cluster1_order,
##      cluster2_order_from1,
##      cluster2_order_from3,
##      cluster3_order_from2,
##      cluster3_order_from4,
##      cluster4_order_from3,
##      cluster4_order_from5,
##      cluster5_order,
##      file="files/cellorder-1.RData")
load("files/cellorder-1.RData")

print("cluster orderings done \n\n")



## C1.pst <- build_PST(cluster1_order, 1, n=100)
## C2.pst <- build_PST(cluster2_order_from1, 2, n=100,
##                     cluster_orders2=cluster2_order_from3)
## C3.pst <- build_PST(cluster3_order_from2, 3, n=100,
##                     cluster_orders2=cluster3_order_from4)
## C4.pst <- build_PST(cluster4_order_from3, 4, n=100,
##                     cluster_orders2=cluster4_order_from5)
## C5_6.pst <- build_PST(cluster5_order, c(5,6), n=100)



## save(C1.pst,
##      C2.pst,
##      C3.pst,
##      C4.pst,
##      C5_6.pst,
##      file="files/cellorder-2.RData")

load("files/cellorder-2.RData")

print("PST building done \n\n")

## C1.gen <- repeat_gen(C1.pst, from1to2[1:100,3],1,30)
## C2.gen <- repeat_gen(C2.pst,
##                      c(generate_pairorders(cluster2_order_from1, n=100)[,1],
##                        generate_pairorders(cluster2_order_from3, n=100)[,as.numeric(nclass_cluster[2])]),
##                      2,
##                      30)


## C3.gen <- repeat_gen(C3.pst,
##                      c(generate_pairorders(cluster3_order_from2, n=100)[,1],
##                        generate_pairorders(cluster3_order_from4, n=100)[,as.numeric(nclass_cluster[3])]),
##                      3,
##                      30)

## C4.gen <- repeat_gen(C4.pst,
##                      c(generate_pairorders(cluster4_order_from3, n=100)[,1],
##                        generate_pairorders(cluster4_order_from5, n=100)[,as.numeric(nclass_cluster[4])]),
##                      4,
##                      30)

## C5_6.gen <- repeat_gen(C5_6.pst,
##                        from4to5[1:100,4],
##                        c(5,6),
##                        30)

## save(C1.gen,
##      C2.gen,
##      C3.gen,
##      C4.gen,
##      C5_6.gen,
##      file="files/cellorder-3.RData")

load("files/cellorder-3.RData")


print("C.gen done \n\n")

## C1.order <- get_unique_order(C1.gen)
## C2.order <- get_unique_order(C2.gen)
## C3.order <- get_unique_order(C3.gen)
## C4.order <- get_unique_order(C4.gen)
## C5_6.order <- get_unique_order(C5_6.gen)

## C1.order <- predicted_order(C1.pst, generate_pairorders(cluster1_order), 1)

## cellorder <- c(C1.order,
##                C2.order+cumsum(nclass_cluster)[1],
##                C3.order+cumsum(nclass_cluster)[2],
##                C4.order+cumsum(nclass_cluster)[3],
##                C5_6.order+cumsum(nclass_cluster)[4])

## write.table(cellorder, "/mnt/gtklab01/ahjung/bivalent/results/sc_time_cellorder.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\n")

## save(C1.order,
##      C2.order,
##      C3.order,
##      C4.order,
##      C5_6.order,
##      cellorder,
##      file="files/cellorder-4.RData")

load("files/cellorder-4.RData")

print("Cellorder all done \n\n")

###################

## hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(20)

## ## compare before and after cell ordering

## pdf("/mnt/gtklab01/ahjung/bivalent/figures/compare_cellorder.pdf",
## width=10, height=10)

## gplots::heatmap.2(hmat_bin,
##           trace="none",
## ColSideColors=cols[sc_time_coldata$exp],
##                                         #as.numeric(a))],#cols[cl1][order(cl1)],
##           col=hmcol,
##           Colv=F,
## #         Rowv=F,
##           dendrogram = "none"
##           )

## gplots::heatmap.2(hmat_bin[,cellorder],
##           trace="none",
## ColSideColors=cols[sc_time_coldata$exp][cellorder],
##                                         #as.numeric(a))],#cols[cl1][order(cl1)],
##           col=hmcol,
##           Colv=F,
## #         Rowv=F,
##           dendrogram = "none"
##           )

## dev.off()



######################### HMM

##### get matrix for HMM

hmat_bin_cellorder <- hmat_bin[,cellorder]

## vit_df <- data.frame(sapply(1:nrow(hmat_bin_cellorder),
##                             function(x) run_viterbi(hmat_bin_cellorder[x,],
##                                                     rownames(hmat_bin_cellorder)[x])))

## vit_df_path <- data.frame(sapply(1:nrow(hmat_bin_cellorder),
##                             function(x) run_viterbi_path(hmat_bin_cellorder[x,],
##                                                     rownames(hmat_bin_cellorder)[x])))

## vit_idx <- (!is.na(vit_df[1,])) & (!is.na(vit_df_path[1,]))

## vit_mat <- as.matrix(vit_df[,vit_idx])
## rownames(vit_mat) <- colnames(hmat_bin)[cellorder]

## vit_mat_path <- as.matrix(vit_df_path[,vit_idx])
## rownames(vit_mat_path) <- colnames(hmat_bin)[cellorder]


## save(vit_df,
##      vit_df_path,
##      vit_mat,
##      vit_mat_path,
##      file="files/viterbi.RData")

load("files/viterbi.RData")

print("Viterbi done \n\n")


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



## x <- get_hmap(as.character(gwhen[,"gene"]),t(vit_mat), "switch_group")
## x2 <-  get_hmap(as.character(gwhen[,"gene"]),t(vit_mat), "change_window")




################# DESEQ
## library(DESeq2)
## dds <- DESeqDataSetFromMatrix(round(bulk_time), bulk_time_coldata, design=~ idx + exp)
## dds <- DESeq(dds)
## resultsNames(dds) # lists the coefficients

## res <- results(dds, name="exp_96h_vs_12h")
## #res$padj < 0.05

## dds_normalized <- counts(dds, normalized=TRUE)

## bulk_time_ave <- cbind("12h"=apply(dds_normalized[,c(1,2,3)],1,mean),
## "24h"=apply(dds_normalized[,c(4,5,6)],1,mean),
## "36h"=apply(dds_normalized[,c(7,8,9)],1,mean),
## "72h"=apply(dds_normalized[,c(10,11,12)],1,mean),
## "96h"=apply(dds_normalized[,c(13,14,15)],1,mean))

## save(bulk_time_ave,
##      file="files/bulk_time_ave.RData")

load("files/bulk_time_ave.RData")

hmap <- get_hmap(unique(gwhen$gene), hmat[,cellorder], "switch_group")
hmap_bin <- get_hmap(unique(gwhen$gene), hmat_bin_cellorder, "switch_group")
hmap_hmm <- get_hmap(unique(gwhen$gene), t(vit_mat), "switch_group")
hmap_bulk <- get_hmap(unique(gwhen$gene), log(as.matrix(bulk_time_ave)+1,10), "switch_group", bulk=TRUE)


p <- plot_hmap(hmap, TRUE)
p_bin <- plot_hmap(hmap_bin, TRUE)
p_hmm <- plot_hmap(hmap_hmm, TRUE)
p_bulk <- plot_hmap(hmap_bulk, TRUE)

pdf("/mnt/gtklab01/ahjung/bivalent/figures/sc_time_switch.pdf",
width=20, height=50)

grid.arrange(p,p_bin,p_hmm,p_bulk,ncol=4)
          
dev.off()


save(gene_switches_df,
     switches_genes,
     gene_change_df_table,
     gwhen,
     hmap,
     hmap_bin,
     hmap_hmm,
     hmap_bulk,
     p,
     p_bin,
     p_hmm,
     p_bulk,
     file="files/hmap.RData")


#save.image(file="files/20190615.RData")
#load("files/20190615.RData")
