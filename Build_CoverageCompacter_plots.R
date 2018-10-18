
library(tidyverse)


inpath_m="N:/Documents/MEGA/home_share/merged_bam_analysis/CoverageCompacterOutput"
inpath_s="N:/Documents/MEGA/home_share/WGS_bam_single_analysis/CoverageCompacterOutput"
files_s <- dir(inpath_s, pattern = ".bed_")

# MAKE TEMPLATE FOR DF
df <- tibble("Sample"=character(), "chr"=character(), "start"=integer(), "end"=integer(), "size"=double(),
             "firstCoveredBase"=integer(), "lastCoveredBase"=integer(), "meanCoverage"=double(), "NBasesCovered"=integer(),
             "DepthSum"=integer(), "coverageFraction"=double(), sep='\t')

# READ IN FILES
for (i in seq_along(files_s)) { 
  # get the string for the chrom column
  sample <- unlist(strsplit(unlist(strsplit(files_s[[i]], "bed_"))[2], '.txt'))[1]
  print(sample)
  # read in file and add new col, then append to df
  df <- read_delim(paste0(inpath_s,'/',files_s[[i]]), '\t') %>% mutate(Sample=sample) %>% rbind(df)
}

# sort df by chrom according to factor
chrFactor <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M')
df <- df %>% arrange(Chrom) %>% mutate(Chrom = factor(Chrom, chrFactor))
df$Chrom <- as.character(df$Chrom)
