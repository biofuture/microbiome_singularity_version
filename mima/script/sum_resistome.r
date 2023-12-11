#################################################################
# Function: calculate total resistance gene and generate violin plots  
# Call: Rscript sum_resistome_violin.R -m meta.txt  -i column_name_of_sampleID -g Group -a normalize_16s.type.tab.txt -b normalize_cellnumber.type.tab.txt -o out_dir
# R packages used: "optparse","ggpubr"
# R version: R-3.5.1
# Last update: 2019-04-05, HuiMin Zheng, zhenghuimin91@126.com
# Copyright: ZHOU HongWei
#################################################################

## install necessary libraries
doInstall<-"FALSE" # Change to FALSE if you don't want packages installed.
toInstall<-c("optparse","ggpubr")
if(doInstall){install.packages(toInstall, repos = "http://cran.us.r-project.org")}
toLib<-toInstall
lapply(toLib, library, character.only = TRUE)

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
  make_option(c("-m","--meta"), type="character", help="Input the metadata file [required]."),
  make_option(c("-i", "--sid"), type="character", help="Input the column name of SampleID[required]."),
  make_option(c("-g", "--group"), type="character", help="Input the column name of Group[required]."),
  make_option(c("-a", "--r16s"), type="character", help="Input the resistome.normalize_16s.type.tab.txt[required]."),
  make_option(c("-b","--cellnum"),type="character",help="Input the resistome.normalize_cellnumber.type.tab.txt[required]."),
  make_option(c("-o", "--out_dir"), type="character", help="Output directory[required].")
)


opt_parser <- OptionParser(option_list=option_list) 
opts <- parse_args(opt_parser, args=args)

# check the parse
if (is.null(opts$meta)|is.null(opts$sid)|is.null(opts$group)|is.null(opts$r16s)|is.null(opts$cellnum)|is.null(opts$out_dir)){
  print_help(opt_parser)
  stop("Please make sure all arguments required are supplied", call.=FALSE)
}


##main 
meta <- read.table(opts$meta,header = T,sep = "\t",comment.char = "",check.names = F)
r_16s <- read.table(opts$r16s,header = T,sep = "\t",comment.char = "",skip = 1,row.names = 1)
r_cn <- read.table(opts$cellnum,header = T,sep = "\t",comment.char = "",skip = 1,row.names = 1)

resistome_sum <- function(resistometable,colname_sum,meta,colname_sid){
  t <- data.frame(t(resistometable))
  t[,colname_sum] <- apply(t,1,sum)
  t <- subset(t,select=colname_sum)
  dta <- merge(meta,t,by.x = colname_sid,by.y = "row.names")
  return(dta)
}

dta <- resistome_sum(r_16s,"sum_16s",meta,opts$sid)
dta <- resistome_sum(r_cn,"sum_cn",dta,opts$sid)

dir.create(opts$out_dir)
library(ggpubr)
comp <- combn(levels(as.factor(dta[,opts$group])),2,list)
ggviolin(dta, x = opts$group, y = "sum_16s", fill = opts$group,
         palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = comp)
ggsave(paste(opts$out_dir,"/sum16s_",opts$group,".pdf",sep =""))
 
ggviolin(dta, x = opts$group, y = "sum_cn", fill = opts$group,
         palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = comp)
ggsave(paste(opts$out_dir,"/sumcellnumb_",opts$group,".pdf",sep =""))

