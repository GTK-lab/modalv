setwd("~/projects/bivalent/modalv")

source("GSE75748_data.R")
source("GSE75748_function.R")

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

set.seed(100)

bimocondition[is.na(bimocondition)] <- FALSE

bimocondition_bi <- (rownames(log_sc_time_tpm) %in% bigenes) & bimocondition
bimocondition_k4 <- (rownames(log_sc_time_tpm) %in% k4genes) & bimocondition
bimocondition_k27 <- (rownames(log_sc_time_tpm) %in% k27genes) & bimocondition

bimocondition_bi <- bimocondition_bi[!is.na(bimocondition_bi)]
bimocondition_k4 <- bimocondition_k4[!is.na(bimocondition_k4)]
bimocondition_k27 <- bimocondition_k27[!is.na(bimocondition_k27)]


tsne_out_bi <- Rtsne(mim_cell,
                  pca=FALSE,
                  perplexity=30,
                  theta=0.0)

tsne1_bi <- ggplot(data=data.frame("C1"=tsne_out_bi$Y[,1],
                    "C2"=tsne_out_bi$Y[,2],
           "type"=mclass),
       aes(C1,C2,col=type)) +
    geom_point() + theme_bw() + scale_color_manual(values=cols)

rd1 <- tsne_out_bi$Y[,1:2]
rownames(rd1) <- mclass

cl1 <- Mclust(rd1,5)$classification
## cl1[cl1==4] <- 6
## cl1[cl1==5] <- 4
## cl1[cl1==6] <- 5
#plot(rd1_pc, pch=16, asp = 1,col=cols[cl1])

tsne2 <- ggplot(data=data.frame("C1"=tsne_out_bi$Y[,1],"C2"=tsne_out_bi$Y[,2],
           "type"=as.factor(cl1)),
       aes(C1,C2,col=type)) +
    geom_point() + theme_bw() + scale_color_manual(values=cols)

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

get_mi <- function(mat, idx) {
    x <- sapply(c(1:nrow(mat))[-idx],
                function(x) mutinformation(mat[idx,], mat[x,]))
    return(append(x, 0, after=idx-1))
}


hmmcondition <- (bimocondition_k27 | bimocondition_bi) & !bimocondition_k4

hmat <- matrix(as.numeric(log_sc_time_tpm[hmmcondition,]),
               nrow=sum(hmmcondition))
rownames(hmat) <- rownames(log_sc_time_tpm[hmmcondition,])
colnames(hmat) <- colnames(log_sc_time_tpm)


hmat_bin <- log_sc_time_tpm_messup_bin[names(hmmcondition)[hmmcondition],]*1
hmat_bin <- hmat_bin[!(is.na(hmat_bin[,1])),]



mim_cell <- sapply(1:ncol(hmat_bin), function(x) get_mi(t(hmat_bin), x))
#mim_ts <- sapply(1:nrow(hmat_bin_both), get_ts_mi)
colnames(mim_cell) <- colnames(hmat_bin)
rownames(mim_cell) <- colnames(hmat_bin)

get_dpt <- function(cluster) {
    dm <- DiffusionMap(mim_cell[cl1%in%cluster, cl1%in%cluster])
    dpt <- DPT(dm)
    return(dpt)
}

nclass_cluster <- table(cl1)
### within cluster order
#dpt_1 <- get_dpt(1)
dpt_2 <- get_dpt(2)
dpt_3 <- get_dpt(3)
dpt_4_5 <- get_dpt(c(4,5))

plot(dpt_2, col=sc_time_coldata$exp[cl1%in%c(2)])
plot(dpt_3)
plot(dpt_4_5, col=as.character(sc_time_coldata$exp[cl1%in%c(4,5)]))


find_neig_tips <- function(mim_tips) {
    mim_tips_m <- melt(mim_tips)
    df <- mim_tips_m[order(mim_tips_m$value,decreasing=TRUE)[1:2],]
#    which(colnames(mim_cell)==find_neig_tips(mim_tip_00_12)[1])
    return(df)
}

mim_tip_2_3 <- mim_cell[tips(dpt_2)+cumsum(nclass_cluster)[1], tips(dpt_3)+cumsum(nclass_cluster)[2]]
mim_tip_3_4_5 <- mim_cell[tips(dpt_3)+cumsum(nclass_cluster)[2], tips(dpt_4)+cumsum(nclass_cluster)[3]]
#mim_tip_4_5 <- mim_cell[tips(dpt_4)+cumsum(nclass_cluster)[3], tips(dpt_5)+cumsum(nclass_cluster)[4]]


tips_list <- list(find_neig_tips(mim_tip_2_3),
                  find_neig_tips(mim_tip_3_4_5))
 

## connect_clusters <- function(tips1, tips2) {
##     tips1[1,2] == tips2[1,2]

## }

ctips <- c("H9.12h_098", "H9.24h_063", "H9.24h_043", "H9.36h_119")

order_1 <- order(as.numeric(mim_cell[cl1==1,"H9.12h_098"]),decreasing=FALSE) # use how similar the cell is to the next tip and order in increasing order
order_2 <- order(dpt_2@dm@eigenvectors[,1],decreasing=TRUE)
order_3 <- order(dpt_3@dm@eigenvectors[,1],decreasing=TRUE)
order_4_5 <- order(dpt_4_5@dm@eigenvectors[,1],decreasing=FALSE)


cellorder <- c(order_1,
               order_2+cumsum(nclass_cluster)[1],
               order_3+cumsum(nclass_cluster)[2],
               order_4_5+cumsum(nclass_cluster)[3])


hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(20)

hmat <- matrix(as.numeric(log_sc_time_tpm[bimocondition_bi,]),
               nrow=sum(bimocondition_bi))
rownames(hmat) <- rownames(log_sc_time_tpm[bimocondition_bi,])
colnames(hmat) <- colnames(log_sc_time_tpm)



heatmap.2(t(vit_mat),
          trace="none",
          ColSideColors=cols[sc_time_coldata$exp][cellorder],
                                        #as.numeric(a))],#cols[cl1][order(cl1)],
          col=hmcol,
          Colv=F,
#          Rowv=F,
          dendrogram = "none"
          )


## geneslist <- LEF1[LEF1 %in% as.character(gwhen$gene[(gwhen$when>150)&(gwhen$when<250)]) ]
## geneslist <- geneslist[geneslist %in% c(upgenes, intupgenes)]

## geneslist <- NFAT[NFAT %in% as.character(gwhen$gene[(gwhen$when>250)])]
## geneslist <- geneslist[geneslist %in% c(upgenes, intupgenes)]

## c(uporder, intuporder, connectororder)[c(uporder, intuporder, connectororder)%in%matrisome[!matrisome %in% c(tgfb,LEF1,NFAT)]]


## +


## geneslist <- c(c(uporder, intuporder, connectororder)[c(uporder, intuporder, connectororder)%in%lef1_geneslist],
##                c(uporder, intuporder, connectororder)[c(uporder, intuporder, connectororder)%in%nfat_geneslist],
##                c(uporder, intuporder, connectororder)[c(uporder, intuporder, connectororder)%in%tgfb])
## #c(uporder, intuporder, connectororder)[c(uporder, intuporder, connectororder)%in%matrisome[!matrisome %in% c(tgfb,LEF1,NFAT)]])
                                       
## pdf("/mnt/gtklab01/ahjung/bivalent/figures/test2.pdf",
## width=10, height=10)

## plotgenes <- c(rownames(log_sc_time_tpm)[str_detect(rownames(log_sc_time_tpm), "FZD")],rownames(log_sc_time_tpm)[str_detect(rownames(log_sc_time_tpm), "WNT")])

## plotgenes <- rownames(log_sc_time_tpm)[grep("^SP[0-9]",rownames(log_sc_time_tpm))]#rownames(log_sc_time_tpm)[str_detect(rownames(log_sc_time_tpm), "SP")]

## pdf("/mnt/gtklab01/ahjung/bivalent/figures/final_sc_time_plotgenes.pdf",
## width=10, height=10)

## heatmap.2(rbind(t(vit_mat)[c(unique(plotgenes)),],
##           (log_sc_time_tpm[c(plotgenes2),cellorder]>0)*1)
##          ,
##           trace="none",
##           ColSideColors=cols[sc_time_coldata$exp][cellorder],
##                                         #as.numeric(a))],#cols[cl1][order(cl1)],
##           col=hmcol,
##           Colv=F,
##         Rowv=F,
##           dendrogram = "none"
##           )


## heatmap.2(hmat_bin[as.character(unique(gsa_mat_m$Var1))[as.character(unique(gsa_mat_m$Var1)) %in% rownames(hmat_bin)],],
##           trace="none",
##           ColSideColors=cols[sc_time_coldata$exp][cellorder],
##                                         #as.numeric(a))],#cols[cl1][order(cl1)],
##           col=hmcol,
##           Colv=F,
##         Rowv=F,
##           dendrogram = "none"
##           )


## df <-   log_sc_time_tpm_messup_bin[names(bimocondition)[(bimocondition_bi | bimocondition_k27) & bimocondition_k4],]

## df <- log_sc_time_tpm_messup_bin[colnames(vit_df)[which(is.na(vit_df[1,])) [which(is.na(vit_df[1,])) %in% which(is.na(vit_df_path[1,]))]][1:100],]

## heatmap.2(df[,cellorder]*1,
##           trace="none",
##           Colv=F,
##           col=hmcol, dendrogram="none")


## dev.off()

## plotgenes <- c("MSX2","NODAL","GATA3","NOG","MSX1","BAMBI","EOMES","MIXL1","WLS","FST","SP5","CPM","BMP8A","BMP8B","ACVRL1","IL12RB1","KCNJ5","KCNA7","SLC24A4","CACNG8","RHBG","EREG","SLC27A1","SEC14L4","EPHA2","GATA6","GSC","SHOX2","BMP2","AMIGO2","SHOX2","BMP4","FN1","PLAU","OTX2","LZTS1","CDK6","COL1A1","COL4A1","COL1A2","ERBB4","SALL1","SALL3","PMAIP1","CCDC81","ID2","PRDM1","SOX17","WNT5A","CAMK2D","GLIPR2","SERHL2","PDCD4","LZTS1","FZD5","FZD8", "DUSP4","DUSP5","SFRP1","CDH11","C11orf84","CTDSP1","BMPER","GPR83")

## plotgenes2 <- c("LEF1","NFATC1","SMAD2","SMAD5","POU5F1","WNT3","WNT5B","T","DVL3","PPP3CC")

## plotgenes_pathways <-c() 

## TGFb <- c("CREBBP","BMP2","BMP4","THBS1","INHBE","ACVRL1","BMP8A","BMP8B","FST","ID2","NODAL","PITX2","NOG"),
## ERBB <- c("PIK3CD","PIK3R5","PIK3R2","AKT1","HRAS","BAD","CDKN1B","CAMK2D","ERBB4","EREG"),
## WNT <- c("CCND1","CREBBP","WNT5A","FZD2","FZD5","FZD8","CAMK2D","CHP2","PPP2R5C","DKK1","VANGL2","SFRP1","SFRP2","SOX17","VANGL1"),
## MAPK <- c("RRAS2","CACNG8","DUSP10","DUSP3","DUSP4","DUSP5","SRF","TAOK2")



## heatmap.2(log_sc_time_tpm[c(rownames(log_sc_time_tpm)[grep("MAPK",rownames(log_sc_time_tpm))],
##                             rownames(log_sc_time_tpm)[grep("DUSP",rownames(log_sc_time_tpm))]),
##                           cellorder]
##          ,
##           trace="none",
##           ColSideColors=cols[sc_time_coldata$exp][cellorder],
##                                         #as.numeric(a))],#cols[cl1][order(cl1)],
##           col=hmcol,
##           Colv=F,
## #        Rowv=F,
##           dendrogram = "none"
##           )

## dev.off()


############


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

intupgenes <- colnames(vit_mat_path)[(vit_mat_path[1,]=="OFF")&(vit_mat_path[758,]=="OFF")]
                                        #colnames(vit_mat)[(vit_mat[1,]<0.1)&(vit_mat[758,]<0.1)]
upgenes <- colnames(vit_mat_path)[(vit_mat[1,]<0.001)&(vit_mat[758,]>0.999)]
downgenes <- colnames(vit_mat_path)[(vit_mat[1,]>0.999)&(vit_mat[758,]<0.001)]


############# find connectors
mim_gene <- sapply(1:nrow(hmat_bin_cellorder), function(x) get_mi(hmat_bin_cellorder, x))
colnames(mim_gene) <- rownames(hmat_bin_cellorder)
rownames(mim_gene) <- rownames(hmat_bin_cellorder)

mim_con <- mim_gene[c(upgenes, intupgenes),!colnames(mim_gene) %in% c(upgenes, intupgenes)]

connectors <- unique(unlist(apply(mim_con, 1, function(x) names(which(x>0.1)))))
connectors <- connectors[connectors %in% colnames(vit_mat)]
connectors <- connectors[!connectors %in% c(upgenes, intupgenes, downgenes)]


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


order_by_switch <- function(switches_df_idx, whichswitch) {
    cat(switches_df_idx,"\n")

    gswhen <- as.numeric(unlist(lapply(whichswitch[switches_genes[[switches_df_idx]]],
                                       function(x) x[1])))

    gsorder <- order(gswhen)
    gsdf <- data.frame("gene"=switches_genes[[switches_df_idx]][gsorder],
                       "when"=gswhen[gsorder],
                       "when_cell"=colnames(hmat_bin_cellorder)[gswhen[gsorder]],
                       "switch_group"=switches_df_idx)
    
    return(gsdf)
}


order_by_cluster <- function(switches_df_idx) {
    cat(switches_df_idx,"\n")

    sgenes <- switches_genes[[switches_df_idx]]
    if (!length(sgenes) < 2) {

    colgenes <- gsub("-","[.]",colnames(vit_mat))

    intupclust <- hclust(dist(t(vit_mat)[colgenes %in% sgenes,]))
    gsorder <- intupclust$label[intupclust$order]
    gsdf <- data.frame("gene"=gsorder,
                       "when"=NA,
                       "when_cell"=colnames(hmat_bin_cellorder)[gsorder],
                       "switch_group"=switches_df_idx)
} else {

    gsdf <- data.frame("gene"=sgenes,
                       "when"=NA,
                       "when_cell"=colnames(hmat_bin_cellorder)[gsorder],
                       "switch_group"=switches_df_idx)
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

gwhen <- order_genes(1)
for (i in 2:nrow(gene_switches_df_table)) {
    gwhen <- rbind(gwhen, order_genes(i))
}

gwhen$gene <- as.character(gwhen$gene)
gwhen <- unique(gwhen)

## gwhen <- data.frame(gene=c(geneorder),
##                      when = c(geneorder_num),
##                      when_cell=c(colnames(hmat_bin_cellorder)[geneorder_num]))

get_hmap <- function(genes, binmat) {
#### vit_mat for _hmm, hmat_bin for _bin, and hmat for logtpm
#    genes <- unique(c(geneorder, connectors))
    genes <- gsub("[.]","-",genes)

    hmap <- melt(binmat[genes,])
    colnames(hmap) <- c("gene","when_cell_all","value")

    gwhen_sub <- gwhen[gwhen$gene %in% genes,]

    hmap_merge <- merge(hmap, gwhen_sub, all=T, by=c("gene"))
    hmap_merge <- merge(hmap_merge, gene_switches_df_table[,c("start_bin","end_bin","switch_group")])
    hmap_merge$when_cell <- factor(hmap_merge$when_cell, levels=levels(hmap_merge$when_cell)[cellorder])
    return(hmap_merge)
}


plot_hmap <- function(myhmap) {
    p <- ggplot(myhmap, aes(y=gene, x=when_cell_all)) +
        geom_tile(aes(fill = value),
                  colour = ifelse(as.character(myhmap$when_cell_all) == as.character(myhmap$when_cell), "red", NA),
                  size=1) +
                      scale_fill_gradient(high=rev(brewer.pal(7,"YlGnBu"))[1],
                                          low="white",na.value="white")+
                                              facet_grid(
                                                         scales="free",
                                                         space="free",
                                                  rows=vars(myhmap$switch_group)
                                                         )+
                                                             theme(axis.text.y = element_text(hjust = 1,
                                                                       size=12,
                                                                       face="bold"),
                                                                   plot.background=element_blank(),
                                                      axis.text.x=element_blank(),
                                                                   axis.ticks.x=element_blank(),
                                                                   legend.position="none") +
                                                                       xlab("Pseudotime Ordered Cells")+                                        scale_x_discrete(position = "top") + ylab("Genes")

return(p)

}

grid.arrange(plot_hmap(get_hmap(gwhen[gwhen$switch_group%in%gene_switches_df_table[as.logical(gene_switches_df_table$start_bin),"switch_group"],"gene"],t(vit_mat))),
             plot_hmap(get_hmap(gwhen[gwhen$switch_group%in%gene_switches_df_table[!as.logical(gene_switches_df_table$start_bin),"switch_group"],"gene"],t(vit_mat))),ncol=2)

plot_hmap(get_hmap(gwhen[gwhen$switch_group%in%c(1,2,3),"gene"],t(vit_mat)))


hmap_hmm <- get_hmap(gwhen$gene,
                     t(vit_mat))

hmap <- get_hmap(unique(c(geneorder, connectors)), hmat[,cellorder])
hmap_bin <- get_hmap(unique(c(geneorder, connectors)), hmat_bin_cellorder)
#rownames(vit_mat) <- colnames(hmat_bin)
hmap_hmm <- get_hmap(unique(c(geneorder, connectors)), t(vit_mat))
hmap_bulk <- get_hmap(unique(c(geneorder, connectors)), log(as.matrix(bulk_time_ave)+1, 10))
#hmap_bulk[(hmap_bulk$Var1 %in% dds_significant) & (hmap_bulk$Var2 == "96h"), "when"] <- 1


p <- plot_hmap(hmap)
p_bin <- plot_hmap(hmap_bin)
p_hmm <- plot_hmap(hmap_hmm)
p_bulk <- plot_hmap(hmap_bulk)

pdf("/mnt/gtklab01/ahjung/bivalent/figures/sc_time_switch.pdf",
width=20, height=15)

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

## par(mfrow=c(2,1))
## plotMA(res)
## plotMA(res[colnames(vit_mat)[colnames(vit_mat) %in% rownames(res)],])

## dds_significant <-  rownames(res)[res$padj < 0.05]

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
    fcol <- ifelse(fgene %in% upgenes,
                   "blue",
                   ifelse(fgene %in% intupgenes,
                          "green",
                          ifelse(fgene %in% connectors,
                                 "black",
                                 "red")))
    tcol <- ifelse(tgene %in% upgenes,
                   "blue",
                   ifelse(tgene %in% intupgenes,
                          "green",
                          ifelse(tgene %in% connectors,
                                 "black",
                                 "red")))
    return(ifelse(fcol != tcol, fcol, "black"))
} 


change_edge_color <- function(edge) {
    fgene <- Rgraphviz::from(edge)
    tgene <- Rgraphviz::to(edge)
    if ((nrow(gwhen[gwhen$gene == fgene,])==0)|(nrow(gwhen[gwhen$gene == tgene,])==0)) {
        "grey" } else if ((is.na(gwhen$when[gwhen$gene==fgene])) | (is.na(gwhen$when[gwhen$gene==tgene]))) {
            return("grey")
        } else if ((gwhen$when[gwhen$gene==fgene] - gwhen$when[gwhen$gene==tgene])<(-10)) {
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


mim_genes <- as.character(unique(gsa_mat_m_overlap$gene))

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

## NFAT <- c(
## "BMP4","PDLIM3","DLC1","CDH6","ERBB4","HAS2","GATA6","PRDM1","FN1","COLEC12","VAMP5","H1F0","EOMES","FOS","VANGL1","C1orf61","BCAR3","AMOTL1","B3GNT5","BHLHE40","BMP2","GATA3","LGALS1","CNTNAP2","PPT2","VEGFB","ST8SIA4","NOG","PMP22","PPFIBP1","PITX2","ZNF521","CDK6","NFIB","FST","CYR61","SHOX2","ST8SIA1","PDCD4","COL1A2","CNN1","WNT5A")

## LEF1 <- c("PITX2","TMEM57","OTX2","ZNF521","CDK6","MSX1","DUSP4","NFIB","NRP2","C11orf84","FST","CYR61","SHOX2","MAF","NEFM","COMMD3","ST8SIA1","PDCD4","CDH11","GNG11","CDKN1C","DDIT4","BTG2","H1F0","CDH2","DLL3","EOMES","FOS","VANGL1","DKK1","C1orf61","BCAR3","AMOTL1","GLIPR2","CTDSP1","TSKU","RBM24","TET2","CDH3","PLK2","B3GNT5","PLOD2","GABRB2","CYP1B1","PIK3R5","MET","TMOD2")

## lef1_geneslist <- LEF1[LEF1 %in% as.character(gwhen$gene[(gwhen$when>150)&(gwhen$when<300)]) ]
## #lef1_geneslist <- lef1_geneslist[lef1_geneslist %in% c(upgenes, intupgenes)]

## nfat_geneslist <- NFAT[NFAT %in% as.character(gwhen$gene[(gwhen$when>300)])]
## #nfat_geneslist <- nfat_geneslist[nfat_geneslist %in% c(upgenes, intupgenes)]

## tgfb <- c("BMP2","BMP4","FST","NODAL","NOG","PITX2","ACVRL1")

## matrisome <- c("FN1","LAMB1","COL4A1","BMP2","BMP4","WNT5A","VEGFB","CYR61","COL1A2","PLAU","BMPER","COCH","FBN2","TSKU","EFEMP2","MATN3","FST","NODAL","LGALS1","PLOD2","LEPREL1","MMP15","SERPINB6","SERPINE2","TNFSF9","COLEC12")



gcolors <- brewer.pal(5,"Set3")
gcolor_df <- data.frame(name=c(mim_genes), col="white", stringsAsFactors=FALSE)
gcolor_df$col[gcolor_df$name %in% as.character(levels(gsa_mat_m_overlap$Var1))] <- gcolors[1]
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

          
pdf("/mnt/gtklab01/ahjung/bivalent/figures/test.pdf",
width=20, height=10)

Rgraphviz::plot(mygraph_mim, nodeAttrs = nAttrs, edgeAttrs = eAttrs)
          
dev.off()



############# find genes with similar switching time
gwhen_sub <- gwhen[!gwhen$gene %in% downgenes,]

genesA <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 90) & (gwhen_sub$when <= 105)])
genesB <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 106) & (gwhen_sub$when <= 179)])
genesC <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 180) & (gwhen_sub$when <= 200)])
genesD <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 201) & (gwhen_sub$when <= 239)])
genesE <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 240) & (gwhen_sub$when <= 280)])
genesF <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 271) & (gwhen_sub$when <= 379)])
genesG <- as.character(gwhen_sub$gene[(gwhen_sub$when >= 380) & (gwhen_sub$when <= 440)])


## plotgenes <- c("T","EOMES","BAMBI","NOG","WNT3","WNT5B","FOS","DKK1","PITX2","SP5","LEF1","SP8","GSC","COL1A2","SLC24A4","WNT5A","NFATC1","PPP3CC","SMAD2","GATA3","POU5F1","NODAL","FZD8","FZD5","BMP2","BMP4","GPR83")

plotgenes_bin <- plotgenes[plotgenes %in% colnames(vit_mat)]
plotgenes_nobin <- plotgenes[!plotgenes %in% colnames(vit_mat)]

plotgenes_bin_hmap <- get_hmap(plotgenes_bin, t(vit_mat))
plotgenes_nobin_hmap <- get_hmap(plotgenes_nobin, log_sc_time_tpm[,cellorder])

q_bin <- plot_hmap(plotgenes_bin_hmap)
q_nobin <- plot_hmap(plotgenes_nobin_hmap)



df <- data.frame("cells"=colnames(log_sc_time_tpm),
           "value"=rep(1, ncol(log_sc_time_tpm)),
           "cluster"=cl1,
           "timepoint"=sc_time_coldata$exp)

df_ordered <- df[cellorder,]

q_cells <- ggplot(df_ordered, aes(x=1:758, y=value)) +
        geom_tile(aes(fill = timepoint)) +
            scale_fill_manual(values=cols) + theme_bw() +
                theme(legend.position="bottom",
                      plot.margin=margin(l=0.3, r=-0.3,unit="cm"))

                      

grid.arrange(q_bin, q_nobin, q_cells, ncol=1)


##################################### GSA analysis
gsa_res <- read.table("/mnt/gtklab01/ahjung/bivalent/results/GSA/sc_time_bik27_notk4_KEGG_q0.001_all.csv",
                      fill=TRUE, sep="\t" ,header=FALSE,skip=9,nrows=51)
gsa_res <- gsa_res[-1,1:7]
colnames(gsa_res) <- c("geneset_name","genes_in_geneset","description","genes_in_overlap","k_K","pvalue","FDRqvalue")

gsa_mat <- read.table("/mnt/gtklab01/ahjung/bivalent/results/GSA/sc_time_bik27_notk4_KEGG_q0.001_all.csv",
                      fill=TRUE, sep="\t" ,header=TRUE,skip=64,nrows=1200, stringsAsFactors=FALSE, quote="")

#pathways <- colnames(gsa_mat)[1:32][grep("SIGNALING", colnames(gsa_mat)[1:32])]
pathways <- colnames(gsa_mat)[unique(c(grep("PATHWAY", colnames(gsa_mat)),
                                grep("SIGNALING", colnames(gsa_mat)),
                                grep("ADHESION", colnames(gsa_mat)),
                                grep("MATRISOME", colnames(gsa_mat))))][c(3,5,8,2,9,13,23,26,27)]


gsa_mat_tf <- gsa_mat[,pathways] != ""
rownames(gsa_mat_tf) <- gsa_mat$Gene.Symbol

cols_pathways <- brewer.pal(length(pathways), "Set3")

gsa_mat_m <- melt(gsa_mat_tf[apply(gsa_mat_tf,1,sum)>0,])

gsa_mat_m_overlap <- merge(gsa_mat_m,
                           data.frame("Var1"=names(apply(gsa_mat_tf,1,sum)),
                                      "overlap"=as.numeric(apply(gsa_mat_tf,1,sum))))
colnames(gsa_mat_m_overlap) <- c("gene","pathway","value","overlap")

gsa_mat_m_overlap <- merge(gsa_mat_m_overlap, gwhen[,c("gene","switch_group")],all=T)

gsa_mat_m_overlap$gene <- factor(gsa_mat_m_overlap$gene,
                                 levels = levels(gsa_hmap$gene)[levels(gsa_hmap$gene) %in% levels(gsa_mat_m_overlap$gene)])

gsa_mat_m_overlap <- gsa_mat_m_overlap[!is.na(gsa_mat_m_overlap$gene),]
gsa_mat_m_overlap <- gsa_mat_m_overlap[!is.na(gsa_mat_m_overlap$pathway),]

gsa_overlap <- ggplot(gsa_mat_m_overlap[!is.na(gsa_mat_m_overlap$switch_group),],
                      aes(y=gene, x=pathway)) +
    geom_tile(aes(fill = ifelse(value, overlap, 0))) + 
        scale_fill_gradient(low = "white", high = "black") +
            scale_x_discrete(labels=c(
                                 "TGFB",
                                 "SMAD2_3N",
                                 "P53",
                                 "P53_DOWNSTREAM",
                                 "MAPK",
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
               plot.margin = margin(0, 1, 1, 0, "cm")) +
                      facet_grid(
                          scales="free",
                          space="free",
                          rows=vars(gsa_mat_m_overlap[!is.na(gsa_mat_m_overlap$switch_group),"switch_group"]),
                          )+  
                    ylab("") + xlab("")

######### MKNK1 and SOCS7 expression could not be binarized
gsa_hmap <- hmap_hmm[hmap_hmm$gene %in% colnames(vit_mat)[colnames(vit_mat) %in% as.character(unique(gsa_mat_m$Var1))],]


## gsa_hmap$Var1 <- factor(gsa_hmap$gene,
##                         levels = genelevels)

## gsa_heatmap <- ggplot(gsa_hmap,
##                       aes(x=Var2, y=Var1)) +
##         geom_tile(aes(fill = value),
##                   colour = c(NA,"red")[as.factor(gsa_hmap[gsa_hmap$Var1 %in% c(upgenes, downgenes, intupgenes, connectors),]$when)],size=1) +
##                       scale_fill_gradient(high=rev(brewer.pal(7,"YlGnBu"))[1],
##                                           low="white",na.value="white")+
##                                                              theme(axis.text.y = element_text(hjust = 1,
##                                                                        size=12,
##                                                                        face="bold"),
##                                                                    plot.background=element_blank(),
##                                                       axis.text.x=element_blank(),
##                                                                    axis.ticks.x=element_blank(),
##                                                                    legend.position="none") +
##                                                                        xlab("Pseudotime Ordered Cells")+                                        scale_x_discrete(position = "top") + ylab("Genes")



plot_hmap(gsa_hmap)



## gsa_heatmap <- ggplot(gsa_hmap,
##                       aes(x=when_cell_all, y=gene)) +
##     geom_tile(aes(fill = value)) +
##                   scale_fill_gradient(high=rev(brewer.pal(7,"YlGnBu"))[1],
##                                       low="white",na.value="white")+
##                                           theme(axis.text.y = element_text(hjust = 1,
##                                                     size=10,
##                                                     face="bold"),
##                                                 plot.background=element_blank(),
##                                                 ## axis.ticks.x=element_blank(),
##                                                 legend.position="none") +
##                                                     xlab("Pseudotime Ordered Cells")+
##                                                         ylab("") +



pdf("/mnt/gtklab01/ahjung/bivalent/figures/test.pdf",
width=10, height=20)

l <- plot_hmap(gsa_hmap)
l <- l + theme( plot.margin = margin(3, 0, 0, 0, "cm"))
grid.arrange(l, gsa_overlap, ncol=2)

dev.off()



################################################

edges_switch <- buildEdgeList(mygraph_mim_switch)
edge_df_switch <- data.frame("from"=unlist(lapply(edges_switch, Rgraphviz::from)),
                      "to"=unlist(lapply(edges_switch, Rgraphviz::to)))

get_neigh_nodes <- function(gene) {
return(as.character(unique(unlist(c(edge_df_switch[(edge_df_switch$from==gene)|(edge_df_switch$to==gene),])))))
}

cat(get_neigh_nodes("PMAIP1"))

gwhen[gwhen$gene %in% get_neigh_nodes("ID2"),"when"]




popular_nodes <- names(sapply(colnames(mim_switch),
                              function(x) length(get_neigh_nodes(x)))[order(sapply(colnames(mim_switch), function(x) length(get_neigh_nodes(x))),decreasing=TRUE)[1:30]])

gsa_genes <- unique(gsa_mat_m_overlap$Var1) #c("CREBBP","BMP2","BMP4","BMP7","INHBE","ACVRL1","BMP8A","BMP8B","FST","ID2","ID3","NODAL","PITX2","LEFTY2","NOG","CDK4","CDK6","FOXO1","TCF7","PLCG1","DVL2","WNT5A","FZD2","FZD5","FZD8","FOS","BID","CAMK2D","CHP2","PPP2R5C","DKK1","VANGL2","SFRP1","SFRP2","SOX17","VANGL1","HSPB1","GADD45B","MAP3K1","RRAS2","MRAS","CACNG8","DUSP10","DUSP1","DUSP3","DUSP4","DUSP5","DUSP7","RASGRF2","SRF","TAOK2","CDKN1B","CCNE1","CCNE2","THBS1","KDR","FLT1","FLT4","FLNC","ITGB5","COL1A1","COL1A2","COL6A2","CAPN2","ZYX","CD40","HLA-B","HLA-C","ICAM1","SDC1","SDC3","HLA-F","HLA-DPB1","CLDN7","CDH2","MADCAM1","CDH3","CNTNAP2","ICAM3","PVR","PTPRF","MKNK1","AKT1","PIK3CD","PIK3R5","PIK3R2","CCND1","ITGA2","FN1","COL4A1","LAMB1","HRAS","BAD","PDGFA","VEGFA","VEGFB","MET","DDIT4","PMAIP1","ZMAT3","SPRY2","SOCS7","ERBB4","EREG","RHEB","SOCS3","SOCS2","SLC2A4","PHKG2","FLOT1","EXOC7","RHOQ","RICTOR", "IFNGR2","IL15","IL4R","IL11","IL12RB1")
mim_genes <- unique(c(popular_nodes, as.character(gsa_genes)))
