########### import and process data

mydir <- "/mnt/gtklab01/ahjung/bivalent-tmp/GSE75748/"

bulk_cell <- read.csv(paste0(mydir,"GSE75748_bulk_cell_type_ec.csv"),row.names=1)
bulk_time <- read.csv(paste0(mydir,"GSE75748_bulk_time_course_ec.csv"),row.names=1)
sc_cell <- read.csv(paste0(mydir,"GSE75748_sc_cell_type_ec.csv"),row.names=1)
sc_time <- read.csv(paste0(mydir,"GSE75748_sc_time_course_ec.csv"),row.names=1)

extract_coldata_sub <- function(df,i){
    df_split <- strsplit(colnames(df)[i],"_|[.]")[[1]]

    coldata <- data.frame("cell"=df_split[1],
                          "exp"=df_split[2],
                          "idx"=df_split[3])
    return(coldata)}

extract_coldata <- function(df) {
    coldata <- extract_coldata_sub(df,1)
    invisible(sapply(2:ncol(df),
                     function(i) coldata <<- rbind(coldata,
                                                   extract_coldata_sub(df,i))))
    return(coldata)}

sc_cell_coldata <- extract_coldata(sc_cell)
sc_time_coldata <- extract_coldata(sc_time)
bulk_time_coldata <- extract_coldata(bulk_time)
bulk_cell_coldata <- data.frame("cell"=strsplit(colnames(bulk_cell)[1],"_")[[1]][1],
                              "rep"=strsplit(colnames(bulk_cell)[1],"_")[[1]][2])
invisible(sapply(2:ncol(bulk_cell),
       function(i) bulk_cell_coldata <<- rbind(bulk_cell_coldata,
                                               data.frame("cell"=strsplit(colnames(bulk_cell)[i],"_")[[1]][1],
                                                          "rep"=strsplit(colnames(bulk_cell)[i],"_")[[1]][2]))))

