## Venn diagram to compare m5C replicates and overlap

# set directories if necessary #
# parent directory
home_dir = "~/OneDrive/Xie_lab/Bioinformatics/"
project = "NA_project"

# location of 5mC files (filtered, overlap, and union files saves here)
setwd(paste0(home_dir,"RNA_BS/", project,"/m5C_siteLists/"))

# desired save directory for figures
save_dir=paste0(home_dir,"RNA_BS/", project,"/Figures/")
dir.exists(save_dir)

FilteredFiles <- list.files(pattern = "_RNAseqFiltered.csv$")
print(RNABS_files[1])

for(i in FilteredFiles) {
  #print(i)
  assign(paste0("rep_",match(i,FilteredFiles)),i)
}
rm(i)


replicate_venn <- function (rep1, rep2) {
  name <- gsub("_[0-9]_RNAseqFiltered.*$","",rep_1)
  rep1 <- read.csv(rep1)
  rep2 <- read.csv(rep2)
  overlap <- read.csv(paste0(name,"_overlap.csv"))
  
  overlap <- nrow(overlap)
  rep1_genes_red <- nrow(subset(rep1, !(group %in% rep2$group)))
  rep2_genes_blue <-nrow(subset(rep2, !(group %in% rep1$group)))
  
  ## if rep2 has more than rep1 use inverted=TRUE
  grid.newpage()  
  if ((rep1_genes_red) > (rep2_genes_blue)){
    venn.plot <- draw.pairwise.venn((rep1_genes_red+overlap), (rep2_genes_blue+overlap), 
                                    overlap, lty=0,fill=c('red','blue'),inverted=TRUE)
  }
  else {
    venn.plot <- draw.pairwise.venn((rep1_genes_red+overlap), (rep2_genes_blue+overlap), 
                                    overlap, lty=0,fill=c('red','blue'),inverted=FALSE)
  }
  
  
  grid.draw(venn.plot)
  #dev.off()
}


replicate_venn(rep_1,rep_2)
