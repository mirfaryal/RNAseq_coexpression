# Prepare to Run ----

tic <- as.integer(as.POSIXct(Sys.time()))
#numeric timestamp based on system time

library("tidyverse")
library("recount")
#load required packages
#to run package recount for MacOS use command line script brew install libpng and Rscript install.packages(latticeExtra)


# Functions for Counts, CPM, RPKM, and TPM from recount2 rse ----

#Get counts data_frame from recount2 rse
get_counts <- function(rse){
  libtype_factor <- rep(1, length(colData(rse)$paired_end))
  libtype_factor[colData(rse)$paired_end] <- 2
  
  round(sweep(assays(rse)$counts, 2, libtype_factor*colData(rse)$avg_read_length, "/")) %>%
    as.data.frame()
}

#Calculate cpm from recount2 rse
get_cpm <- function(tis_gene_count) {
  sweep(tis_gene_count, 2, c(colSums(tis_gene_count)/(10**6)), "/") %>%
    apply(., 2, round, digits = 5) %>%
    as.data.frame()
}

#Calculate tpm from recount2 rse 
get_tpm <- function(tis_gene_count){
  rpk <- sweep(tis_gene_count, 1, c((bplength/1000)), "/")
  sweep(rpk, 2, c((colSums(rpk)/(10**6))), "/") %>%
    apply(., 2, round, digits = 5) %>%
    as.data.frame()
}

#Calculate rpkm from recount2 rse
get_rpkm <- function(tis_gene_count){
  sweep(get_cpm(tis_gene_count), 1, c((bplength/1000)), "/") %>%
    apply(., 2, round, digits = 5) %>%
    as.data.frame()
}


# Make Some Directories if They Don't Exist ----

dirname <- "./rse"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}
#make directory named rse if it does not exist

dirname <- "./counts"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}
#make directory named counts if it does not exist

dirname <- "./metadata"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}
#make directory named metadata if it does not exist

dirname <- "./rowranges"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}
#make directory named rowranges if it does not exist

dirname <- "./rpkm"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}
#make directory named rpkm if it does not exist

dirname <- "./cpm"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}
#make directory named cpm if it does not exist

dirname <- "./tpm"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}
#make directory named tpm if it does not exist


# Get Selected Projects ----

selected_projects_sample_tissue <- read_delim("../data/gtex_project_sample_run_tissue.txt",
                                              delim = "\t",
                                              col_names = T)
#read txt file with projects, tissues, samples, and run (each sample has unique run number) 
#into selected_projects_sample_tissue
#for MacOS use selected_projects_sample_tissue <- read_delim("./data/gtex_project_sample_run_tissue.txt", delim = "\t", col_names = T)

selected_projects <- selected_projects_sample_tissue %>%
  pull(project) %>%
  unique()
#make unique projects = selected_projects vector (will be length one - gtex all one project)


# Download rse_gene and write row and column metadata ----

SRP <- "SRP012682" #assign project to object SRP (might it be more efficient to assign object value selected_project)
print(SRP)
#options(timeout = 250) 
#if timeout issues
download_study(SRP, type = "rse-gene") #download data from recount
file.rename(paste0(SRP,"/rse_gene.Rdata"),
            paste0("./rse/", SRP, "_rse_gene.Rdata")) #rename rse_gene for given SRP
#download and rename rse_gene for given SRP

load(file.path("./rse", paste0(SRP, "_rse_gene.Rdata")))
#load the downloaded data files

temp = rse_gene %>%
  colData() %>%
  as.data.frame() %>%
  write_delim(paste0("./metadata/", paste(SRP, "metadata.txt", sep = "_")),
              delim = "\t",
              col_names = T)
#keep all metadata for rse, file = "SRP_metadata.txt"
#ERROR: Error in stream_delim_(df, file, ..., bom = bom, quote_escape = quote_escape,  : Don't know how to handle vector of type list in column 'characteristics'.

rse_gene %>%
  rowRanges() %>%
  as.data.frame() %>%
  write_delim(paste0("./rowranges/", paste(SRP, "rowranges.txt", sep = "_")),
              delim = "\t",
              col_names = T)
#keep rowRanges for rse (length of gene, etc), file = "SRP_rowranges.txt"

bplength <- rse_gene %>%
  rowRanges() %>%
  as_tibble() %>%
  pull(bp_length)
#extract bp_length from rowranges


# Get gene_counts Per Tissue and Normalize Counts Matrix ----

gene_counts <- rse_gene %>%
  get_counts()
#convert rse coverage to actual gene counts with get_counts()

tissues <- selected_projects_sample_tissue %>%
  filter(project == SRP) %>%
  pull(sharq_beta_tissue) %>%
  unique()
#pull unique tissues from each indiviual project 

for(tis in tissues) {
  runs <- selected_projects_sample_tissue %>%
    filter(project == SRP & sharq_beta_tissue == tis) %>%
    pull(run)
  #get samples corresponding to project and tissue
  
  gene_counts %>% 
    select(intersect(runs, colnames(gene_counts))) %>% 
    rownames_to_column(var = "gene") %>%
    write_delim(paste0("./counts/", paste(SRP, tis, "counts.txt", sep = "_")),
                delim = "\t",
                col_names = T)
  #get gene_counts for given project and tissue
  
  tissue_gene_count <- read_delim(paste0("./counts/", paste(SRP, tis, "counts.txt", sep = "_")),
                                  delim = "\t",
                                  col_names = T,
                                  col_types = cols(.default = "d", gene = "c")) %>%
    column_to_rownames("gene")
  #tissue_gene_count = argument for get_cpm, get_tpm, and get_rpkm
  
  get_cpm(tissue_gene_count) %>%
    rownames_to_column(var = "gene") %>%
    write_delim(paste0("./cpm/", paste(SRP, tis, "cpm.txt", sep = "_")), 
                delim = "\t", 
                col_names = T)
  #get cpm, write to file
  
  get_rpkm(tissue_gene_count) %>%
    rownames_to_column(var = "gene") %>%
    write_delim(paste0("./rpkm/", paste(SRP, tis, "rpkm.txt", sep = "_")), 
                delim = "\t", 
                col_names = T)
  #get rpkm, write to file
  
  get_tpm(tissue_gene_count) %>%
    rownames_to_column(var = "gene") %>%
    write_delim(paste0("./tpm/", paste(SRP, tis, "tpm.txt", sep = "_")), 
                delim = "\t", 
                col_names = T)
  #get tpm, write to file
}
#filter for the samples (runs) we actualely want from each project (using project and desired tissue from that project)

toc <- as.integer(as.POSIXct(Sys.time()))
#numerical timestamp based on system time
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
#print runtime
