library(ggplot2)
library(gridExtra)
library(mixtools)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(diptest)
library(doParallel)

test_bimodal <- function(exp) {
a <- dip.test(as.numeric(exp))
return(a$p.value)
}

fit_bimodal <- function(exp,name) {
    tryCatch(mixmdl <- normalmixEM(exp, k=2, mu=c(0,6),maxit = 10),error=function(e){NA})
    
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


messup <- function(sc_cell) {
    # sc_cell must be in log scale
    nnoise <- matrix(rnorm(nrow(sc_cell)*ncol(sc_cell),
                           mean=0, sd=0.1),ncol=ncol(sc_cell))
    sc_cell_messup <- sc_cell+nnoise
    return(sc_cell_messup)
}

fit_bimodal_multi <- function(ncores, logexp_messup) {
#ncores <- detectCores()-2
cl <- makeCluster(ncores) # create a cluster with max-2 cores
registerDoParallel(cl) # register the cluster

res_tpm = foreach(i = 1:nrow(logexp_messup),
    .combine = "rbind",
              .packages="mixtools")     %dopar% {
  # generate a bootstrap sample              
        fit_bimodal(logexp_messup[i,],rownames(logexp_messup)[i])
}

stopCluster(cl) # shut down the cluster

return(res_tpm)

}

filter_condition <- function(res_diptest, res_tpm, cond_dip, cond_labmda) {
condition <- (logtpm_messup_bimodal <= cond_dip) &
    sapply(1:nrow(res_tpm), function(i) min(res_tpm[i,]$mu1, res_tpm[i,]$mu2) <= 1) &
        sapply(1:nrow(res_tpm), function(i) max(res_tpm[i,]$mu1, res_tpm[i,]$mu2) > 1)  &
            sapply(1:nrow(res_tpm), function(i) res_tpm[i,]$lambda1 <= cond_labmda)
return(condition)

}
