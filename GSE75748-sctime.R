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
cl1[cl1==4] <- 6
cl1[cl1==5] <- 4
cl1[cl1==6] <- 5
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


bi_bin <- log_sc_time_tpm_messup_bin[names(bimocondition_bi)[bimocondition_bi],]*1



my.mst <- mstree(dist(rd1),1)
mygraph <- graph_from_adjacency_matrix(neig2mat(my.mst))
mygraph <- minimum.spanning.tree(mygraph)
x <- get.shortest.paths(mygraph,1,758)
a <- shortest.paths(mygraph, 1) # starting from any time from time 0
a2 <- get.shortest.paths(mygraph,1) # need to keep the option to specify end point

longest <- which(unlist(lapply(a2$vpath,length))==max(unlist(lapply(a2$vpath,length))))

V(mygraph)$color <- cols[sc_time_coldata$exp]
E(mygraph, path=a2$vpath[[longest]])$color <- "red"
V(mygraph)[as.numeric(a2$vpath[[longest]])]$color <- "black"

pdf("/mnt/gtklab01/ahjung/bivalent/figures/sc_time_pseudotime.pdf",
    width=8, height=10)

par(mfrow=c(1,1))
plot(mygraph,
     vertex.size=2,
#     layout=rd1,
     vertex.label=NA,
     edge.arrow.size=0.5)

dev.off()

















library(infotheo)
library(minet)
library(Rgraphviz)

get_mi <- function(mat, idx) {
    x <- sapply(c(1:nrow(mat))[-idx],
                function(x) multiinformation(t(mat[c(idx,x),])))
    return(append(x, 0, after=idx-1))
}




mim_cell <- sapply(1:ncol(bi_bin), function(x) get_mi(t(bi_bin), x))
#mim_ts <- sapply(1:nrow(hmat_bin_both), get_ts_mi)
colnames(mim_cell) <- colnames(bi_bin)
rownames(mim_cell) <- colnames(bi_bin)



nclass_cluster <- table(cl1)
### within cluster order
mim_00 <- mim_cell[cl1==1,cl1==1]#[sc_time_coldata$exp == "00hb4s",sc_time_coldata$exp == "00hb4s"]
dm_00 <- DiffusionMap(mim_00)
dpt_00 <- DPT(dm_00)
tips_00 <- tips(dpt_00)


mim_12_24 <- mim_cell[cl1%in%c(2,3),cl1%in%c(2,3)]#[sc_time_coldata$exp == "12h",sc_time_coldata$exp == "12h"]
dm_12_24 <- DiffusionMap(mim_12_24)
dpt_12_24 <- DPT(dm_12_24)
tips_12_24 <- tips(dpt_12_24)


mim_36 <- mim_cell[cl1%in%c(4),cl1%in%c(4)]#[sc_time_coldata$exp == "24h",sc_time_coldata$exp == "24h"]
dm_36 <- DiffusionMap(mim_36)
dpt_36 <- DPT(dm_36)
tips_36 <- tips(dpt_36) # no tips.


mim_72_96 <- mim_cell[cl1==5,cl1==5]#[sc_time_coldata$exp == "72h",sc_time_coldata$exp == "72h"]
dm_72_96 <- DiffusionMap(mim_72_96)
dpt_72_96 <- DPT(dm_72_96)
tips_72_96 <- tips(dpt_72_96)


find_neig_tips <- function(mim_tips) {
    mim_tips_m <- melt(mim_tips)
    df <- mim_tips_m[order(mim_tips_m$value,decreasing=TRUE)[1:2],]
#    which(colnames(mim_cell)==find_neig_tips(mim_tip_00_12)[1])
    return(df)
}

mim_tip_00_12_24 <- mim_cell[tips_00, tips_12_24+cumsum(nclass_cluster)[1]]
mim_tip_12_24_36 <- mim_cell[tips_12_24+cumsum(nclass_cluster)[1], tips_36+cumsum(nclass_cluster)[3]]
mim_tip_36_72_96 <- mim_cell[tips_36+cumsum(nclass_cluster)[3], tips_72_96+cumsum(nclass_cluster)[4]]
#mim_tip_72_96 <- mim_cell[tips_72+cumsum(nclass_cluster)[4], tips_96+cumsum(nclass_cluster)[5]]

tips_list <- list(find_neig_tips(mim_tip_00_12_24),find_neig_tips(mim_tip_12_24_36),find_neig_tips(mim_tip_36_72_96))

## connect_clusters <- function(tips1, tips2) {
##     tips1[1,2] == tips2[1,2]

## }

ctips <- c("H9.00hb4s_090", "H9.12h_053","H9.24h_020", "H9.36h_119","H9.36h_111", "H9.72h_101")
which(colnames(mim_cell) == ctips[1])
which(colnames(mim_cell) == ctips[2])-as.numeric(cumsum(nclass_cluster)[1])
which(colnames(mim_cell) == ctips[3])-as.numeric(cumsum(nclass_cluster)[1])
which(colnames(mim_cell) == ctips[4])-as.numeric(cumsum(nclass_cluster)[3])
which(colnames(mim_cell) == ctips[5])-as.numeric(cumsum(nclass_cluster)[3])
which(colnames(mim_cell) == ctips[6])-as.numeric(cumsum(nclass_cluster)[4])

order_00 <- order(dpt_00@dm@eigenvectors[,1],decreasing=TRUE)
order_12_24 <- order(dpt_12_24@dm@eigenvectors[,1],decreasing=FALSE)
order_36 <- order(dpt_36@dm@eigenvectors[,2],decreasing=FALSE)
order_72_96 <- order(dpt_72_96@dm@eigenvectors[,1],decreasing=FALSE)

cellorder <- c(order_00,
order_12_24+cumsum(nclass_cluster)[1],
order_36+cumsum(nclass_cluster)[3],
order_72_96+cumsum(nclass_cluster)[4])


hmcol <- colorRampPalette(brewer.pal(9, "BuPu"))(20)

hmat <- matrix(as.numeric(log_sc_time_tpm[bimocondition_bi,]),
               nrow=sum(bimocondition_bi))
rownames(hmat) <- rownames(log_sc_time_tpm[bimocondition_bi,])
colnames(hmat) <- colnames(log_sc_time_tpm)

hmat_bin <- bi_bin[!(is.na(bi_bin[,1])),]


heatmap.2(hmat_bin[,cellorder],
          trace="none",
          ColSideColors=cols[sc_time_coldata$exp][cellorder],
                                        #as.numeric(a))],#cols[cl1][order(cl1)],
          col=hmcol,
          Colv=F,
#          Rowv=F,
          dendrogram = "none"
          )


geneslist <- LEF1[LEF1 %in% as.character(gwhen$gene[(gwhen$when>150)&(gwhen$when<250)]) ]
geneslist <- geneslist[geneslist %in% c(upgenes, intupgenes)]

geneslist <- NFAT[NFAT %in% as.character(gwhen$gene[(gwhen$when>250)])]
geneslist <- geneslist[geneslist %in% c(upgenes, intupgenes)]

c(uporder, intuporder, connectororder)[c(uporder, intuporder, connectororder)%in%matrisome[!matrisome %in% c(tgfb,LEF1,NFAT)]]


+


geneslist <- c(c(uporder, intuporder, connectororder)[c(uporder, intuporder, connectororder)%in%lef1_geneslist],
               c(uporder, intuporder, connectororder)[c(uporder, intuporder, connectororder)%in%nfat_geneslist],
               c(uporder, intuporder, connectororder)[c(uporder, intuporder, connectororder)%in%tgfb])
#c(uporder, intuporder, connectororder)[c(uporder, intuporder, connectororder)%in%matrisome[!matrisome %in% c(tgfb,LEF1,NFAT)]])
                                       
pdf("/mnt/gtklab01/ahjung/bivalent/figures/test2.pdf",
width=10, height=10)

plotgenes <- c(rownames(log_sc_time_tpm)[str_detect(rownames(log_sc_time_tpm), "FZD")],rownames(log_sc_time_tpm)[str_detect(rownames(log_sc_time_tpm), "WNT")])

plotgenes <- rownames(log_sc_time_tpm)[grep("^SP[0-9]",rownames(log_sc_time_tpm))]#rownames(log_sc_time_tpm)[str_detect(rownames(log_sc_time_tpm), "SP")]

pdf("/mnt/gtklab01/ahjung/bivalent/figures/final_sc_time_plotgenes.pdf",
width=10, height=10)

heatmap.2(rbind(t(vit_mat)[c(unique(plotgenes)),],
          (log_sc_time_tpm[c(plotgenes2),cellorder]>0)*1)
         ,
          trace="none",
          ColSideColors=cols[sc_time_coldata$exp][cellorder],
                                        #as.numeric(a))],#cols[cl1][order(cl1)],
          col=hmcol,
          Colv=F,
        Rowv=F,
          dendrogram = "none"
          )


heatmap.2(hmat_bin[as.character(unique(gsa_mat_m$Var1))[as.character(unique(gsa_mat_m$Var1)) %in% rownames(hmat_bin)],],
          trace="none",
          ColSideColors=cols[sc_time_coldata$exp][cellorder],
                                        #as.numeric(a))],#cols[cl1][order(cl1)],
          col=hmcol,
          Colv=F,
        Rowv=F,
          dendrogram = "none"
          )


df <-   log_sc_time_tpm_messup_bin[names(bimocondition)[(bimocondition_bi | bimocondition_k27) & bimocondition_k4],]

df <- log_sc_time_tpm_messup_bin[colnames(vit_df)[which(is.na(vit_df[1,])) [which(is.na(vit_df[1,])) %in% which(is.na(vit_df_path[1,]))]][1:100],]

heatmap.2(df[,cellorder]*1,
          trace="none",
          Colv=F,
          col=hmcol, dendrogram="none")


dev.off()

plotgenes <- c("MSX2","NODAL","GATA3","NOG","MSX1","BAMBI","EOMES","MIXL1","WLS","FST","SP5","CPM","BMP8A","BMP8B","ACVRL1","IL12RB1","KCNJ5","KCNA7","SLC24A4","CACNG8","RHBG","EREG","SLC27A1","SEC14L4","EPHA2","GATA6","GSC","SHOX2","BMP2","AMIGO2","SHOX2","BMP4","FN1","PLAU","OTX2","LZTS1","CDK6","COL1A1","COL4A1","COL1A2","ERBB4","SALL1","SALL3","PMAIP1","CCDC81","ID2","PRDM1","SOX17","WNT5A","CAMK2D","GLIPR2","SERHL2","PDCD4","LZTS1","FZD5","FZD8", "DUSP4","DUSP5","SFRP1","CDH11","C11orf84","CTDSP1","BMPER","GPR83")

plotgenes2 <- c("LEF1","NFATC1","SMAD2","SMAD5","POU5F1","WNT3","WNT5B","T","DVL3","PPP3CC")

plotgenes_pathways <-c() 

TGFb <- c("CREBBP","BMP2","BMP4","THBS1","INHBE","ACVRL1","BMP8A","BMP8B","FST","ID2","NODAL","PITX2","NOG"),
ERBB <- c("PIK3CD","PIK3R5","PIK3R2","AKT1","HRAS","BAD","CDKN1B","CAMK2D","ERBB4","EREG"),
WNT <- c("CCND1","CREBBP","WNT5A","FZD2","FZD5","FZD8","CAMK2D","CHP2","PPP2R5C","DKK1","VANGL2","SFRP1","SFRP2","SOX17","VANGL1"),
MAPK <- c("RRAS2","CACNG8","DUSP10","DUSP3","DUSP4","DUSP5","SRF","TAOK2")



heatmap.2(log_sc_time_tpm[c(rownames(log_sc_time_tpm)[grep("MAPK",rownames(log_sc_time_tpm))],
                            rownames(log_sc_time_tpm)[grep("DUSP",rownames(log_sc_time_tpm))]),
                          cellorder]
         ,
          trace="none",
          ColSideColors=cols[sc_time_coldata$exp][cellorder],
                                        #as.numeric(a))],#cols[cl1][order(cl1)],
          col=hmcol,
          Colv=F,
#        Rowv=F,
          dendrogram = "none"
          )

dev.off()


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
hmmcondition <- (bimocondition_k27 | bimocondition_bi) & !bimocondition_k4

hmat <- matrix(as.numeric(log_sc_time_tpm[hmmcondition,]),
               nrow=sum(hmmcondition))
rownames(hmat) <- rownames(log_sc_time_tpm[hmmcondition,])
colnames(hmat) <- colnames(log_sc_time_tpm)


hmat_bin <- log_sc_time_tpm_messup_bin[names(hmmcondition)[hmmcondition],]*1
hmat_bin <- hmat_bin[!(is.na(hmat_bin[,1])),]

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
    return(ifelse(onoroff=="ON", switch_on[1], switch_off[1]))
}

## find_switch <- function(genex) {
## ## plot(as.numeric(bimodals[bimodals[,1]==genex,-1]))
## ## abline(v=which(diff(vit_mat[,genex])==max(diff(vit_mat[,genex]))))
## return(which(abs(diff(vit_mat[,genex]))==max(abs(diff(vit_mat[,genex])))))
## }

## find_switch_int <- function(genex) {
## ## plot(as.numeric(bimodals[bimodals[,1]==genex,-1]))
##     ## abline(v=which(diff(vit_mat[,genex])==max(diff(vit_mat[,genex]))))

##     if (any(((vit_mat[,genex]>0.4) &
##                          (vit_mat[,genex]<0.6) &
##                              (c(0, diff(vit_mat[,genex]))>0)))) {
##         return(min(which((vit_mat[,genex]>0.4) &
##                          (vit_mat[,genex]<0.6) &
##                              (c(0, diff(vit_mat[,genex]))>0))))
##     } else {
##         return(min(which(diff(vit_mat[,genex])>0.4)))
##     }
## }

connectorswitch <- sapply(connectors, function(x) find_switch(x, "ON"))
intupswitch <- sapply(intupgenes, function(x) find_switch(x, "ON"))
upswitch <- sapply(upgenes, function(x) find_switch(x, "ON"))
downswitch <- sapply(downgenes, function(x) find_switch(x, "OFF"))

connectorclust <- hclust(dist(t(vit_mat)[c(connectors),]))
connectororder <- connectorclust$label[connectorclust$order]
intupclust <- hclust(dist(t(vit_mat)[c(intupgenes),]))
intuporder <- intupclust$label[intupclust$order]
uporder <- upgenes[order(upswitch)]
downorder <- downgenes[order(downswitch)]
gorder <- c(uporder, downorder, intuporder)


geneorder <- gorder
geneorder_num <- c(as.numeric(upswitch)[order(upswitch)], as.numeric(downswitch)[order(downswitch)], as.numeric(intupswitch)[intupclust$order])
geneorder[geneorder == "HLA.B"] <- "HLA-B"   # dot and dash conversion
geneorder[geneorder == "HLA.DQB1"] <- "HLA-DQB1"
geneorder[geneorder == "HLA.DPB1"] <- "HLA-DPB1"
geneorder[geneorder == "ERVMER34.1"] <- "ERVMER34-1"

colnames(vit_mat)[colnames(vit_mat) == "HLA.B"] <- "HLA-B"   # dot and dash conversion
colnames(vit_mat)[colnames(vit_mat) == "HLA.DQB1"] <- "HLA-DQB1"
colnames(vit_mat)[colnames(vit_mat) == "HLA.DPB1"] <- "HLA-DPB1"
colnames(vit_mat)[colnames(vit_mat) == "ERVMER34.1"] <- "ERVMER34-1"

gwhen <- data.frame(gene=c(geneorder, connectors),
                     when = c(geneorder_num,connectorswitch),
                     when_cell=c(colnames(hmat_bin_cellorder)[c(geneorder_num,connectorswitch)]))

get_hmap <- function(genes, binmat) {
#### vit_mat for _hmm, hmat_bin for _bin, and hmat for logtpm
#    genes <- unique(c(geneorder, connectors))
    hmap <- melt(binmat[genes,])

    cluster_conditional <- ifelse(hmap$Var1 %in% upgenes, "up",
                              ifelse(hmap$Var1 %in% intupgenes,
                                     "intup",
                                     ifelse(hmap$Var1 %in% connectors,
                                            "connectors",
                                     "down")))

    hmap <- cbind(hmap, cluster=cluster_conditional)

    sapply(1:nrow(gwhen), function(x) 
        hmap$when[as.character(hmap$Var1)==gwhen$gene[x] & as.character(hmap$Var2)==gwhen$when_cell[x]] <<- 1)
    hmap$when[is.na(hmap$when)] <- 0
    hmap$cluster <- factor(hmap$cluster, levels = c("down", "intup", "up", "connectors"))
    return(hmap)
}


plot_hmap <- function(hmap) {
    p <- ggplot(hmap, aes(x=Var2, y=Var1)) +
        geom_tile(aes(fill = value),
                  colour = c(NA,"red")[as.factor(hmap$when)],size=1) +
                      scale_fill_gradient(high=rev(brewer.pal(7,"YlGnBu"))[1],
                                          low="white",na.value="white")+
                                              facet_grid(rows = vars(cluster),
                                                         scales="free",
                                                         space="free")+
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
    if ((is.na(gwhen$when[gwhen$gene==fgene])) | (is.na(gwhen$when[gwhen$gene==tgene]))) {
        return("none")
    } else if ((gwhen$when[gwhen$gene==fgene] - gwhen$when[gwhen$gene==tgene])<0) {
        edge@attrs$dir <- "forward"
    } else if ((gwhen$when[gwhen$gene==fgene] - gwhen$when[gwhen$gene==tgene])>0) {
        edge@attrs$dir <- "back"
    } else {
        edge@attrs$arrowhead <- "none"
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
    if ((is.na(gwhen$when[gwhen$gene==fgene])) | (is.na(gwhen$when[gwhen$gene==tgene]))) {
        return("grey")
    } else if ((gwhen$when[gwhen$gene==fgene] - gwhen$when[gwhen$gene==tgene])<(-20)) {
        return(check_edge_color(fgene, tgene))
    } else if ((gwhen$when[gwhen$gene==fgene] - gwhen$when[gwhen$gene==tgene])>20) {
        return(check_edge_color(tgene, fgene))
    } else {
        return("grey")
    }}


change_edge_weight_switch <- function(edge) {
    fgene <- Rgraphviz::from(edge)
    tgene <- Rgraphviz::to(edge)
    if ((is.na(gwhen$when[gwhen$gene==fgene])) | (is.na(gwhen$when[gwhen$gene==tgene]))) {
        return(10)
    } else if (gwhen$when[gwhen$gene==fgene] < gwhen$when[gwhen$gene==tgene]) {
        w <- (gwhen$when[gwhen$gene==tgene] - gwhen$when[gwhen$gene==fgene])
        edge@attrs$weight <- w
    } else if (gwhen$when[gwhen$gene==fgene] > gwhen$when[gwhen$gene==tgene]) {
        w <- (gwhen$when[gwhen$gene==fgene] - gwhen$when[gwhen$gene==tgene])
        edge@attrs$weight <- w
    } else {
        edge@attrs$weight <- 1
    }}

change_edge_weight_mi <- function(edge) {
    fgene <- Rgraphviz::from(edge)
    tgene <- Rgraphviz::to(edge)
    return(mim[fgene,tgene])
}


mim <- mim_gene[unique(plotgenes_pathways),
                unique(plotgenes_pathways)]
#mim <- mim_gene[unique(c(geneslist,"SP8","WNT5A","FOS","GATA3","COL1A2","GSC","WLS","GATA6")),
 #               unique(c(geneslist,"SP8","WNT5A","FOS","GATA3","COL1A2","GSC","WLS","GATA6"))]
net2_o <- aracne(mim)
mygraph_mim <- as(net2_o ,"graphNEL")
gdf <- data.frame(gene=geneorder,
                  updown=ifelse(geneorder %in% upgenes,
                      "blue",
                      ifelse(geneorder %in% intupgenes,
                             "green",
                             "red")))

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
gcolor_df <- data.frame(name=c(geneorder, connectors), col="white", stringsAsFactors=FALSE)
gcolor_df$col[gcolor_df$name %in% lef1_geneslist] <- gcolors[1]
gcolor_df$col[gcolor_df$name %in% nfat_geneslist] <- gcolors[2]
gcolor_df$col[gcolor_df$name %in% matrisome] <- gcolors[3]
#gcolor_df$col[gcolor_df$name %in% colnames(mim)[as.numeric(sp5_all$res[[1]])]] <- gcolors[3]

nAttrs <- list()

## nAttrs$fillcolor <- as.character(gcolor_df$col)
## names(nAttrs$fillcolor) <- gcolor_df$name
nAttrs$fontsize <- rep(20, length(c(geneorder, connectors)))
names(nAttrs$fontsize) <- c(geneorder, connectors)
nAttrs$size <- rep(1, length(c(geneorder, connectors)))
names(nAttrs$size) <- c(geneorder, connectors)
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

          
pdf("/mnt/gtklab01/ahjung/bivalent/figures/final_sc_time_graph_plotgenes.pdf",
width=20, height=10)

Rgraphviz::plot(mygraph_mim, nodeAttrs = nAttrs, edgeAttrs = eAttrs)
          
dev.off()

par(mfrow=c(4,1))
plot(log_sc_time_tpm["LEF1",],main="LEF1")
plot(log_sc_time_tpm["MAZ",],main="MAZ")
plot(log_sc_time_tpm["SP1",],main="SP1")
plot(log_sc_time_tpm["NFATC1",],main="NFATC1")

par(mfrow=c(1,1))

image(t(rbind(c(1:758) %in% gwhen$when[gwhen$gene %in% LEF1],
      c(1:758) %in% gwhen$when[gwhen$gene %in% NFAT],
      c(1:758) %in% gwhen$when[gwhen$gene %in% tgfb],
      c(1:758) %in% gwhen$when[gwhen$gene %in% matrisome])))


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


edge_df <- data.frame("from"=unlist(lapply(edges, Rgraphviz::from)), "to"=unlist(lapply(edges, Rgraphviz::to)))

get_neigh_nodes <- function(gene) {
return(as.character(unique(unlist(c(edge_df[(edge_df$from==gene)|(edge_df$to==gene),])))))
}

cat(get_neigh_nodes("SFRP1"))

gwhen[gwhen$gene %in% get_neigh_nodes("ID2"),"when"]



##################################### GSA analysis
gsa_res <- read.table("/mnt/gtklab01/ahjung/bivalent/results/GSA/sc_time_bik27_notk4_KEGG_q0.001.csv",
                      fill=TRUE, sep="\t" ,header=FALSE,skip=9,nrows=51)
gsa_res <- gsa_res[-1,1:7]
colnames(gsa_res) <- c("geneset_name","genes_in_geneset","description","genes_in_overlap","k_K","pvalue","FDRqvalue")

gsa_mat <- read.table("/mnt/gtklab01/ahjung/bivalent/results/GSA/sc_time_bik27_notk4_KEGG_q0.001.csv",
                      fill=TRUE, sep="\t" ,header=TRUE,skip=64,nrows=1200, stringsAsFactors=FALSE, quote="")

#pathways <- colnames(gsa_mat)[1:32][grep("SIGNALING", colnames(gsa_mat)[1:32])]
pathways <- colnames(gsa_mat)[1:32][
c(grep("SIGNALING", colnames(gsa_mat)[1:32]),grep("ADHESION", colnames(gsa_mat)[1:32]))][c(1,3,4,10,11,2,5,6,7,8,9)]

gsa_mat_tf <- gsa_mat[,pathways] != ""
rownames(gsa_mat_tf) <- gsa_mat$Gene.Symbol

rownames(gsa_mat_tf)[gsa_mat_tf[,1]]


cols_pathways <- brewer.pal(length(pathways), "Set3")

gsa_mat_m <- melt(gsa_mat_tf[apply(gsa_mat_tf,1,sum)>0,])
gsa_mat_m_overlap <- merge(gsa_mat_m,
                           data.frame("Var1"=names(apply(gsa_mat_tf,1,sum)),
                                      "overlap"=as.numeric(apply(gsa_mat_tf,1,sum))))


#gsa_mat_m_overlap$Var1 <- factor(gsa_mat_m_overlap$Var1,
#                                 levels = genelevels)

gsa_overlap <- ggplot(gsa_mat_m_overlap## [gsa_mat_m_overlap$Var1 %in% genelevels[genelevels %in% c(upgenes, downgenes, intupgenes, connectors)],]
                     ,
                      aes(y=Var1, x=Var2)) +
    geom_tile(aes(fill = ifelse(value, overlap, 0))) + 
        scale_fill_gradient(low = "white", high = "black") +
            scale_x_discrete(labels=c("TGFb",
                                 "MAPK",
                                 "WNT",
                                 "FOCAL_ADHESION",
                                 "CELL_ADHESION",
                                 "P53",
                                 "Insulin",
                                 "VEGF",
                                 "JAK/STAT",
                                 "ERBB",
                                 "NEUROTROPHIN"
                                      )) +
    theme_bw() +
        theme(axis.text.y = element_text(hjust = 0.5,
                  size=10,
                  face="bold")) +
                    ylab("") + xlab("")

######### MKNK1 and SOCS7 expression could not be binarized
gsa_hmap <- get_hmap(as.character(unique(gsa_mat_m$Var1)),
                     rbind(hmat_bin_cellorder,
                           "SOCS7"=(log_sc_time_tpm["SOCS7",cellorder] > 1)*1,
                           "MKNK1"=(log_sc_time_tpm["MKNK1",cellorder] > 1)*1,
                           "PVR"=(log_sc_time_tpm["PVR",cellorder] > 1)*1)[as.character(unique(gsa_mat_m$Var1)),]                     )

gsa_hmap$Var1 <- factor(gsa_hmap$Var1,
                        levels = genelevels)


gsa_heatmap <- ggplot(gsa_hmap[gsa_hmap$Var1 %in% genelevels[genelevels %in% c(upgenes, downgenes, intupgenes, connectors)],],
                      aes(x=Var2, y=Var1)) +
    geom_tile(aes(fill = value)) +
                  scale_fill_gradient(high=rev(brewer.pal(7,"YlGnBu"))[1],
                                      low="white",na.value="white")+
                                          theme(axis.text.y = element_text(hjust = 1,
                                                    size=10,
                                                    face="bold"),
                                                plot.background=element_blank(),
                                                ## axis.ticks.x=element_blank(),
                                                legend.position="none") +
                                                    xlab("Pseudotime Ordered Cells")+
                                                        ylab("") +
                                                            scale_x_discrete(
breaks=levels(gsa_hmap$Var2)[c(92,197,259,421)],
labels=as.character(c(92,197,259,421)))


grid.arrange(gsa_heatmap, gsa_overlap, ncol=2)







genelevels <- unique(gsa_mat_m_overlap$Var1) #c("CREBBP","BMP2","BMP4","BMP7","INHBE","ACVRL1","BMP8A","BMP8B","FST","ID2","ID3","NODAL","PITX2","LEFTY2","NOG","CDK4","CDK6","FOXO1","TCF7","PLCG1","DVL2","WNT5A","FZD2","FZD5","FZD8","FOS","BID","CAMK2D","CHP2","PPP2R5C","DKK1","VANGL2","SFRP1","SFRP2","SOX17","VANGL1","HSPB1","GADD45B","MAP3K1","RRAS2","MRAS","CACNG8","DUSP10","DUSP1","DUSP3","DUSP4","DUSP5","DUSP7","RASGRF2","SRF","TAOK2","CDKN1B","CCNE1","CCNE2","THBS1","KDR","FLT1","FLT4","FLNC","ITGB5","COL1A1","COL1A2","COL6A2","CAPN2","ZYX","CD40","HLA-B","HLA-C","ICAM1","SDC1","SDC3","HLA-F","HLA-DPB1","CLDN7","CDH2","MADCAM1","CDH3","CNTNAP2","ICAM3","PVR","PTPRF","MKNK1","AKT1","PIK3CD","PIK3R5","PIK3R2","CCND1","ITGA2","FN1","COL4A1","LAMB1","HRAS","BAD","PDGFA","VEGFA","VEGFB","MET","DDIT4","PMAIP1","ZMAT3","SPRY2","SOCS7","ERBB4","EREG","RHEB","SOCS3","SOCS2","SLC2A4","PHKG2","FLOT1","EXOC7","RHOQ","RICTOR", "IFNGR2","IL15","IL4R","IL11","IL12RB1")
