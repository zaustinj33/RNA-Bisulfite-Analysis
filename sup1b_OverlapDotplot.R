## Dotplot showing similarity of methylation level between two replicates for 
# all transcripts

home_dir = "~/OneDrive/Xie_lab/Bioinformatics/"
project = "NA_project"

setwd(paste0(home_dir,"RNA_BS/", project,"/m5C_siteLists/"))
save_dir=paste0(home_dir,"RNA_BS/", project,"/Figures/")
dir.exists(save_dir)

overlap_files <- list.files(pattern = "_overlap.csv$")


for(k in overlap_files) {
  #print(i)
  assign(paste0("overlap_",match(k, overlap_files)),k)
}
rm(k)

## dotplot ##
replicate_dotplot <- function (overlap) {
  name <- gsub("_overlap.csv$","",overlap_1)
  overlap <- read.csv(overlap)

  
  # only grab methylLevel columns
  methylDF <- overlap[ , grepl( "methylLevel" , names(overlap))]
  
  pearson <- cor.test(methylDF[,1], methylDF[,2], method="pearson")
  ggplot(overlap, aes(x = methylDF[,1], y = methylDF[,2])) + 
    geom_point(alpha=.5, color="black") + theme_bw() + 
    theme(legend.title = element_blank()) + 
    theme(panel.grid=element_blank()) +
    theme(panel.border=element_blank())+
    theme(axis.line = element_line(size=0.5, colour = "black")) +
    xlab("m5c level of rep1") + 
    ylab("m5c level of rep2") +
    stat_cor(method="pearson")
  #png(paste0(name, "_dotplot.png"), device=NULL)
  #return(pearson)
}

replicate_dotplot(overlap_1)

