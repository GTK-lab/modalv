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

############# find connectors
mim_con <- mim_gene[colnames(mim_gene) %in% downstream,
                    !colnames(mim_gene) %in% downstream]

connectors <- unique(unlist(apply(mim_con, 1, function(x) names(which(x>=0.1)))))
connectors <- connectors[connectors %in% colnames(vit_mat)]
connectors <- connectors[!connectors %in% as.character(unique(gsa_mat_m_overlap$gene))]




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
d <- unlist(lapply(split(d, ceiling(seq_along(d)/10)),
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

plot_thmm <- function(gene) {
d1 <- vit_mat_path[,gene]

df <- data.frame("hmm_path"=unlist(split(d1, ceiling(seq_along(d1)/130))),
                 "cwindow"=rep(1:6, lapply(split(d1, ceiling(seq_along(d1)/130)), length)),
                 "cells"=unlist(lapply(split(d1, ceiling(seq_along(d1)/130)), function(x) 1:length(x))))

q <- ggplot(df, aes(y=cells, x=hmm_path, col=hmm_path, fill=hmm_path)) +
    geom_bar(position="dodge",stat="identity") + theme_classic() +
        theme(strip.background = element_blank(),
              strip.text.y = element_blank(),
              legend.position="none") + facet_grid(                          scales="free",
                          space="free",cols=vars(cwindow)) + scale_fill_grey(start=0.8, end=0.2)+  scale_color_grey(start=0.8, end=0.2)  

return(q)

}



pdf("/mnt/gtklab01/ahjung/bivalent/figures/sc_time_example.pdf",width=10, height=6)

grid.arrange(plot_fraction("PITX2"),
#             plot_fraction("LPGAT1"),
#             plot_fraction("EOMES"),
             plot_twindow("PITX2"),
#             plot_twindow("LPGAT1"),
#             plot_twindow("EOMES"),
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
