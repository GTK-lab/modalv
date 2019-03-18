#load("~/projects/bivalent/modalv/010319.RData")

bimonames <- names(bimocondition_bi)[bimocondition_bi][names(bimocondition_bi)[bimocondition_bi] %in% rownames(log_sc_cell_tpm_messup)]

hmat_bi <- log_sc_cell_tpm_messup[bimonames,]
clusters <- cl1[order(cl1)]
ctable <- table(clusters)

log_sc_cell_tpm_messup_bin <- binarizeexp(log_sc_cell_tpm_messup[bimonames,][1,])

invisible(sapply(2:length(bimonames),
                 function(i) log_sc_cell_tpm_messup_bin <<- rbind(log_sc_cell_tpm_messup_bin,
binarizeexp(log_sc_cell_tpm_messup[bimonames,][i,]))))

hmat_bin <- log_sc_cell_tpm_messup_bin[,order(cl1)]
rownames(hmat_bin) <- bimonames

hmat_bin <- hmat_bin[!is.na(hmat_bin[,1]),]

getratio <- function(binexp){
    if (all(binexp)) { # all ON
        return(c(0,length(binexp)))
    } else if (all(!binexp)) { # all OFF
        return(c(length(binexp),0))
    } else {
        return(table(binexp*1))
    }
}

clusterratio <- function(binexp, nclass) {
    l <- length(nclass)
    c <- cumsum(nclass)

    cstart <- c(1,rev(rev(c)[-1])+1)
    cend <- c

    cratio <- t(sapply(1:l, function(x) getratio(binexp[cstart[x]:cend[x]])))
    rownames(cratio) <- names(nclass)
    colnames(cratio) <- c("OFF","ON")
    return(cratio)
}

cblist <- list(c(1,2),c(1,3),c(1,4),c(1,5),c(1,6))

cres <- clusterratio(hmat_bin[1,],ctable)
hmat_pval <- unlist(lapply(cblist, function(x) fisher.test(cres[x,],alternative="two.sided")$p.value))
hmat_odds <- unlist(lapply(cblist, function(x) fisher.test(cres[x,],alternative="two.sided")$estimate))


for (i in 2:nrow(hmat_bin)) {
    cres <- clusterratio(hmat_bin[i,],ctable)

    hmat_pval <- rbind(hmat_pval,unlist(lapply(cblist, function(x) fisher.test(cres[x,][c(2,1),c(2,1)],alternative="two.sided")$p.value
                                             )))
    hmat_odds <- rbind(hmat_odds,unlist(lapply(cblist, function(x) fisher.test(cres[x,][c(2,1),c(2,1)],alternative="two.sided")$estimate
                                               )))
}

hmat_qval <- matrix(p.adjust(hmat_pval, method = "BH", n = length(hmat_pval)),nrow=nrow(hmat_pval))

rownames(hmat_pval) <- rownames(hmat_bin)
rownames(hmat_qval) <- rownames(hmat_bin)
rownames(hmat_odds) <- rownames(hmat_bin)

hmat_cratio <- list(hmat_pval, hmat_qval, hmat_odds)


cres <- clusterratio(hmat_bin[1,],ctable)
npc_fraction <- data.frame("h1"=cres[1,2]/(cres[1,1]+cres[1,2]), "npc"=cres[5,2]/(cres[5,1]+cres[5,2]))

for (i in 2:nrow(hmat_bin)) {
    cres <- clusterratio(hmat_bin[i,],ctable)
    npc_fraction <- rbind(npc_fraction,
                          data.frame("h1"=cres[1,2]/(cres[1,1]+cres[1,2]), "npc"=cres[5,2]/(cres[5,1]+cres[5,2])))
}

rownames(npc_fraction) <- rownames(hmat_bin)

pdf("/mnt/gtklab01/ahjung/bivalent/figures/ft_between_cellclusters_0.01.pdf",
    width=5,height=5)

library(gplots)

for (i in 1:5) {
## subg <- names(which(hmat_cratio[[2]][,i] < 0.01))
## subg <- names(apply(hmat_cratio < 1e-10,1,any))[apply(hmat_cratio < 1e-10,1,any)]

subg <- names(which(hmat_cratio[[2]][,i]<0.01))
print(length(subg))
subg_o <- order(hmat_cratio[[3]][,i][hmat_cratio[[2]][,i] < 0.01],decreasing=TRUE)

upr <- as.numeric(hmat_cratio[[3]][,i][hmat_cratio[[2]][,i] < 0.01][order(hmat_cratio[[3]][,i][hmat_cratio[[2]][,i] < 0.01],decreasing=TRUE)])>1
print(sum(upr))

row_col <- c("blue","red")[c(rep(1,sum(upr)),rep(2,sum(!upr)))]

subcl <- which(cl1 %in% c(1,i+1))
subcl_o <- order(cl1[cl1 %in% c(1,i+1)])

heatmap.2(log_sc_cell_tpm[subg,subcl][subg_o, subcl_o],
          trace="none",
          ColSideColors=cols[cl1][subcl][subcl_o],
          RowSideColors=row_col,
          col=hmcol,
          Colv=F,
          Rowv=F,
          dendrogram = "none",
          density.info="none"
          )



}

dev.off()

pdf("/mnt/gtklab01/ahjung/bivalent/figures/ft_npc_fraction.pdf",
    width=5,height=5)


heatmap.2(as.matrix(npc_fraction[subg,][subg_o,]),
          trace="none",
          ColSideColors=cols[c(1,5)],
          RowSideColors=row_col,
          col=rev(gray.colors(10, start = 0, end = 1, gamma = 2.2, alpha = NULL)),
          Colv=F,
          Rowv=F,
          dendrogram = "none",
          density.info="none"
          )



dev.off()

write.csv(hmat_cratio[[1]], file="/mnt/gtklab01/ahjung/bivalent/figures/hmat_pval_biv.csv", quote=FALSE)

write.csv(hmat_cratio[[2]], file="/mnt/gtklab01/ahjung/bivalent/figures/hmat_qval_biv.csv", quote=FALSE)

write.csv(hmat_cratio[[3]], file="/mnt/gtklab01/ahjung/bivalent/figures/hmat_oddsratio_biv.csv", quote=FALSE)



sc_cell_gcluster <- read.table("/mnt/gtklab01/ahjung/bivalent/figures/sc_cell_bimocondition_bi_geneclusters.txt")




upgenes <- list(
    names(which((hmat_cratio[[2]][,1]<0.01)&(hmat_cratio[[3]][,1]>1))),
    names(which((hmat_cratio[[2]][,2]<0.01)&(hmat_cratio[[3]][,2]>1))),
    names(which((hmat_cratio[[2]][,3]<0.01)&(hmat_cratio[[3]][,3]>1))),
#    names(which((hmat_cratio[[2]][,4]<0.01)&(hmat_cratio[[3]][,4]>1))),
    names(which((hmat_cratio[[2]][,5]<0.01)&(hmat_cratio[[3]][,5]>1))))

commongenes <- upgenes[[4]][upgenes[[4]]%in% upgenes[[3]][upgenes[[3]] %in% upgenes[[1]][upgenes[[1]] %in% upgenes[[2]]]]]


cat(upgenes[upgenes%in%k4genes])

library("venn")



venn(upgenes, ilab=TRUE, zcolor="style")


onetype <- names(which((hmat_cratio[[2]][,4]<0.01)&(hmat_cratio[[3]][,4]>1)))
cat(onetype[!onetype %in% commongenes])
