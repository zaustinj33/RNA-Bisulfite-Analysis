## script for generating filtered overlap and union files for meRanCall output files
"Input: 2 biological replicate files labelled *_[0-9]_m5c* stored in <home>/RNA_BS/m5C_siteLists/
        2 RNAseq gene counts files with the same beginning text as the meRanCall files
        optional 2 RNAseq transcript files 
        RNAseq files stored in <home>/RNAseq/gene_lists

 Output: 2 RNAseq-filtered m5C counts files
          One overlap file
          One union file
          All files contain: coverage, count, methylation level, and seq context info

The user only needs to specify home_dir and project. As long as files are in the
right locations, and files and columns are labelled correctly this will work.
There are no checkpoints to protect the user from these problems, so please hastag 
and write functions until all tests have passed.
"
# set directories #
# parent directory
home_dir = "~/OneDrive/Xie_lab/Bioinformatics/"
project = "NA_project"

# location of 5mC files (filtered, overlap, and union files saves here)
setwd(paste0(home_dir,"RNA_BS/", project,"/m5C_siteLists/"))

# location of RNAseq count files
RNAseq_dir=paste0(home_dir,"RNAseq/", project,"/gene_lists/")
dir.exists(RNAseq_dir)

# desired save directory for figures
save_dir=paste0(home_dir,"RNA_BS/", project,"/Figures/")
dir.exists(save_dir)

## --- read in raw data files --- ##
RNABS_files <- list.files(pattern = "txt$")
print(RNABS_files[1])

for(i in RNABS_files) {
  #print(i)
  assign(paste0("rep_",match(i,RNABS_files)),i)
}
rm(i)

# Import RNAseq data, convert in to memory-friendly format for filtering

get_genes_list <- function (rep) {
  name <- gsub("_BS_.*$","",rep)
  
  ## if tx counts exist, import gene and tx counts (from HISAT count)
  # note: need a better identifier if tx files don't exist
  if (grepl("Xiaoran",name) == TRUE){
    name1 <- gsub("BS","",name1)
    name2 <- gsub("1RNA","2RNA",name1)
    
    gene_name1 <- paste0(name1,"seq_genome_final.genes.results")
    gene_name2 <- paste0(name2,"seq_genome_final.genes.results")
    tx_name1 <- paste0(name1,"seq_genome_final.isoforms.results")
    
    genome_final.genes1 <- read.delim(paste0("~/Dropbox/FA_project/RNAseq/",gene_name1))
    genome_final.genes1 <- genome_final.genes1[genome_final.genes1$TPM > 0,]
    genome_final.genes2 <- read.delim(paste0("~/Dropbox/FA_project/RNAseq/",gene_name2))
    genome_final.genes2 <- genome_final.genes2[genome_final.genes2$TPM > 0,]
    genes_tmp <- merge(genome_final.genes1, genome_final.genes2, by = "gene_id", all = T)
    
    
    genome_final.isoforms1 <- read.delim(paste0("~/Dropbox/FA_project/RNAseq/",tx_name1))
    genome_final.isoforms2 <- read.delim(paste0("~/Dropbox/FA_project/RNAseq/",tx_name2))
    genome_final.isoforms1 <- genome_final.isoforms1[genome_final.isoforms1$TPM > 0,]
    genome_final.isoforms2 <- genome_final.isoforms2[genome_final.isoforms2$TPM > 0,]
    isoforms_tmp <- merge(genome_final.isoforms1, genome_final.isoforms2, by = "tx_id", all = T)
    
    # gene list
    genes_list <- data.frame(genes_tmp[1])
    
    #transcript list
    tx_list <- data.frame(isoforms_tmp[1], isoforms_tmp[2])
  }
  
  else {
    ## if STAR gene counts
    cat("reading STAR RNAseq counts unique to file \n")
    gene_name1 <- paste0(name,"_1_STARoutReadsPerGene.out.tab")
    gene_name2 <- paste0(name,"_2_STARoutReadsPerGene.out.tab")
    
    genome_final.genes1 <- read.delim(paste0(RNAseq_dir,gene_name1), header=F,skip=4)
    genome_final.genes1 <- as.data.frame(subset(genome_final.genes1, rowSums(genome_final.genes1[,2:4]) > 0)[,1])
    colnames(genome_final.genes1) <- c("gene_id")
    genome_final.genes2 <- read.delim(paste0(RNAseq_dir,gene_name2), header=F,skip=4)
    genome_final.genes2 <- as.data.frame(subset(genome_final.genes2, rowSums(genome_final.genes2[,2:4]) > 0)[,1])
    colnames(genome_final.genes2) <- c("gene_id")
    
    genes_list <- merge(genome_final.genes1, genome_final.genes2, by = "gene_id", all = T)
    #tx-ome_exp <- read.delim(#need to use kallisto to get tx counts)
    #ID_list <- subset(tx_gene_names, tx_gene_names$gene_id %in% genome_exp[,1])
    
    # since no tx file exists, convert gene names to tx
    tx_table <- read.csv(paste0(home_dir,"Mouse_tx_gene_symbol_annotations.csv"))
    genes_list <- merge(tx_table, genes_list, by="gene_id")
    genes_list <- genes_list[!duplicated(genes_list$tx_name),]
    tx_list <- as.data.frame(genes_list$tx_name)
    colnames(tx_list) = c("transcript_id")
  }
  return(tx_list)
}

#test
#ID_genes <- get_genes_list(rep_3)


# Subsetting function 
# filters by coverage, counts, and methyl rate, then by RNAseq data
# for clarity this was already done by meRanCall, but is not always complete
RNAseq_filter <- function (rep){
  name <- gsub("_m5c.*$","",rep)
  raw <- read.delim(rep)
  # by coverage, C-count, and methyl rate
  final_tmp <- raw[raw$C_count >= 3 & raw$cov >= 20 & raw$methRate >= 0.1,]
  
  ## by RNAseq, hashtag if no RNAseq files present
  genes_list <- get_genes_list(name)
  #first filter by genes
  final <- merge(final_tmp, genes_list, by.x = "X.SeqID", by.y = "transcript_id")
  #second filter by transcripts
  #final <- merge(final_tmp, genes_list[[1]], by.x = "gene_id", by.y = "gene_id")
  
  #remove na's 
  final <- final[complete.cases(final),]
  final$group <- paste(final$X.SeqID, final$refPos)
  write.csv(final, paste0(name, "_RNAseqFiltered.csv"))
  return(final)
}

#test
#RNAseq_filter(rep_3) -> subrep1
#RNAseq_filter(rep_4) -> subrep2

# Merging replicates 1+2, creating union and overlap csv files 
# overlap
overlap <- function(rep1_file, rep2_file) {
  rep1_name <- gsub("_m5c.*$","",rep1_file)
  rep2_name <- gsub("_m5c.*$","",rep2_file)
  name <- gsub("_[0-9]_m5c.*$","",rep1_file)
  
  rep1 <- read.csv(paste0(rep1_name, "_RNAseqFiltered.csv"))
  rep2 <- read.csv(paste0(rep2_name, "_RNAseqFiltered.csv"))
  
  # Create relevant columns for overlap downstream analysis
  rep1_tmp <- data.frame(rep1$group, rep1$cov, rep1$C_count, rep1$seqContext)
  colnames(rep1_tmp) <- c("group", "rep1_cov", "rep1_C_count", "rep1_seqContext")
  #calculate methyl level
  rep1_tmp$methylLevel <- rep1_tmp$rep1_C_count/rep1_tmp$rep1_cov
  colnames(rep1_tmp) <- c("group", paste0(name,"_1_cov"), paste0(name,"_1_C_count"),
                          paste0(name,"_2_seqContext"), paste0(name,"_2_methylLevel"))
  
  rep2_tmp <- data.frame(rep2$group, rep2$cov, rep2$C_count, rep2$seqContext)
  colnames(rep2_tmp) <- c("group", "rep2_cov", "rep2_C_count", "rep2_seqContext")
  rep2_tmp$methylLevel <- rep2_tmp$rep2_C_count/rep2_tmp$rep2_cov
  colnames(rep2_tmp) <- c("group", paste0(name,"_2_cov"), paste0(name,"_2_C_count"),
                          paste0(name,"_2_seqContext"),paste0(name,"_2_methylLevel"))
  
  overlap <- merge(rep1_tmp, rep2_tmp, by = "group")
  #write.csv(overlap, paste(name, "_overlap.csv",sep=""))
  
  return(overlap)
}
# test
#overlap(rep_3, rep_4) -> sample_overlap

# union
union <- function(rep1, rep2, name){
  rep1_tmp <- data.frame(rep1$group, rep1$cov, rep1$C_count)
  colnames(rep1_tmp) <- c("group", "rep1_cov", "rep1_C_count")
  rep1_tmp$methylLevel <- rep1_tmp$rep1_C_count/rep1_tmp$rep1_cov
  
  rep2_tmp <- data.frame(rep2$group, rep2$cov, rep2$C_count)
  colnames(rep2_tmp) <- c("group", "rep2_cov", "rep2_C_count")
  rep2_tmp$methylLevel <- rep2_tmp$rep2_C_count/rep2_tmp$rep2_cov
  
  unionSites <- merge(rep1_tmp, rep2_tmp, by = "group", all = T)
  unionSites[is.na(unionSites)] = "/"
  
  write.csv(unionSites, paste0(name, "_union.csv"))
  
  return(unionSites)
}

# test
#union(subrep1, subrep2, rep_3) -> sample_union


## Full run function
Create_overlap_files <- function (rep1, rep2) {
  #x <- RNAseq_filter(rep1)
  #y <- RNAseq_filter(rep2)
  final <- overlap(rep1, rep2)
  return(final)
}

test <- Create_overlap_files(rep_1, rep_2)
