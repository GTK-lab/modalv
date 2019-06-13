library(RColorBrewer)
library(reshape2)
library(infotheo)
library(minet)


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

repeat_gen <- function(c.pst, vertices, cluster, nr, ncores=30) {
registerDoParallel(ncores)  # use multicore, set to the number of our cores
c.gen <- foreach (i=1:nr, .combine=rbind) %dopar% {
    cat(i,"\n")
    generate_pst(c.pst, vertices, cluster) }
return(c.gen)
}


get_unique_order <- function(c.gen) {
    return(as.numeric(c.gen[which(apply(c.gen, 1, function(x) sum(duplicated(x))) == 0)[1],]))
}

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
        p <- p +  facet_grid(
            scales="free",
            space="free",
            rows=vars(gorder)) }

return(p)
}
