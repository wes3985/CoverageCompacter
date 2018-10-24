
library(tidyverse)


inpath_m="N:/Documents/MEGA/home_share/merged_bam_analysis/CoverageCompacterOutput"
inpath_s30="N:/Documents/MEGA/home_share/WGS_bam_single_analysis/CoverageCompacterOutput30Kb"
inpath_s="N:/Documents/MEGA/home_share/WGS_bam_single_analysis/CoverageCompacterOutput"
files_s <- dir(inpath_s, pattern = ".bed_")
files_m <- dir(inpath_m, pattern = ".bed_")
files_s30 <- dir(inpath_s30, pattern = ".bed_")

# MAKE TEMPLATE FOR DF
df <- tibble("Sample"=character(), "chr"=character(), "start"=integer(), "end"=integer(), "size"=double(),
             "firstCoveredBase"=integer(), "lastCoveredBase"=integer(), "meanCoverage"=double(), "NBasesCovered"=integer(),
             "DepthSum"=integer(), "coverageFraction"=double(), sep='\t')

# READ IN FILES_s
for (i in seq_along(files_s)) { 
  # get the string for the chrom column
  sample <- unlist(strsplit(unlist(strsplit(files_s[[i]], "bed_"))[2], '.txt'))[1]
  print(sample)
  # read in file and add new col, then append to df
  df <- read_delim(paste0(inpath_s,'/',files_s[[i]]), '\t', na = c("", "NA", "None")) %>% mutate(Sample=sample) %>% rbind(df)
}
# READ IN FILES_m
for (i in seq_along(files_m)) { 
  # get the string for the chrom column
  sample <- unlist(strsplit(unlist(strsplit(files_m30[[i]], "bed_"))[2], '.txt'))[1]
  print(sample)
  # read in file and add new col, then append to df
  df <- read_delim(paste0(inpath_m,'/',files_m[[i]]), '\t', na = c("", "NA", "None")) %>% mutate(Sample=sample) %>% rbind(df)
}

# replace None for NA (if required)
df <- na_if(df, 'None')

# sort df by chrom according to factor
chrFactor <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M')
df <- df %>% arrange(chr) %>% mutate(chr = factor(chr, chrFactor))
df$Chrom <- as.character(df$chr)

# Add col for tissue type
df <- df %>% mutate(Tissue = case_when(grepl("HE",Sample) ~ "HE",
                                       grepl("BE",Sample) ~ "BE",
                                       grepl("AE",Sample) ~ "AE",
                                       grepl("PBMC",Sample) ~ "C",
                                       grepl("FC69",Sample) ~ "C",
                                       grepl("HC",Sample) ~ "HC"))

# calculate relative size for the bins
df <- df %>% group_by(chr) %>% mutate(RelSize= size/((max(end)+1)-min(start)))

# make a binary column for covered vc non-covered bins
df <- df %>% mutate(Covered = case_when(firstCoveredBase=="None" | is.na(firstCoveredBase) ~ "non-covered",
                                        firstCoveredBase!="None" | is.na(firstCoveredBase) ~ "covered"))

# Remove the X,Y,M chrom because it has almost complete coverage
df <- subset(df, !chr %in% c("chrY", "chrX", "chrM"))

# Collect top level coverage data into new tibble
SampleCovSummary <- tibble("Sample"=character(), "Tissue"=character(), "covFrac"=double(), "projectedCovFrac"=double(), 
                           "totalDepthSum"=integer(), "CovMeanCov"=double(), "MeanCov"=double(), "CovRatio"=double(), "alpha"=double())
ChromCovSummary <- tibble("Sample"=character(), "Tissue"=character(), "chr"=character(), "covFrac"=double(), "projectedCovFrac"=double(), 
                     "totalDepthSum"=integer(), "CovMeanCov"=double(), "MeanCov"=double(), "CovRatio"=double(), "alpha"=double())
# (1) Percent coverage by sample
for (s in unique(df$Sample)){
  dfs <- subset(df, Sample==s);
  print(s)
  genomeSize=sum(as.numeric(dfs$size))    # sum of genome size
  covbins=subset(dfs, NBasesCovered!=0)
  cov=sum(as.numeric(covbins$size))
  noCov=subset(dfs, NBasesCovered==0)  
  noCov=sum(as.numeric(noCov$size))
  covFrac=(sum(as.numeric(covbins$NBasesCovered))/sum(as.numeric(dfs$size)))    # actual coverage
  ProjectedCovFrac=(cov/(cov+noCov))      # projected coverage with deeper sequencing
  totalDepthSum=(sum(as.numeric(dfs$DepthSum)))    # Total bases sequenced for sample
  CovMeanCov=totalDepthSum/sum(as.numeric(covbins$NBasesCovered))                  # mean coverage of covered regions (total bases sequencd / total bases covered)
  MeanCov=totalDepthSum/genomeSize                                                 # mean coverage across whole sample
  CovRatio=CovMeanCov/MeanCov                                                   # indicates how reads are distributed across the genome, 
  alpha=covFrac * totalDepthSum
  SampleCovSummary = add_row(SampleCovSummary, Sample=s, Tissue=dfs$Tissue[1], covFrac=covFrac, projectedCovFrac=ProjectedCovFrac, 
          totalDepthSum=totalDepthSum, CovMeanCov=CovMeanCov, MeanCov=MeanCov, CovRatio=CovRatio,  alpha=alpha)

}

SampleCovSummary30kb=SampleCovSummary

ggplot() + geom_point(data=SampleCovSummary30kb, aes(x=CovRatio, y=covFrac, colour=Tissue)) + 
  xlim(0,50) + geom_smooth(data=SampleCovSummary30kb, aes(x=CovRatio, y=covFrac)) + 
  geom_point(data=SampleCovSummary, aes(x=CovRatio, y=covFrac, colour=Tissue)) + 
  geom_smooth(data=SampleCovSummary, aes(x=CovRatio, y=covFrac)) +
  ggtitle("Coverage achieved at ~1x as a function of coverage ratio")


ggplot() + 
  geom_point(data=SampleCovSummary30kb, aes(x=covFrac, y=projectedCovFrac, colour=Tissue)) +
  geom_smooth(data=SampleCovSummary30kb, aes(x=covFrac, y=projectedCovFrac, colour='red')) + 
  geom_point(data=SampleCovSummary, aes(x=covFrac, y=projectedCovFrac, colour=Tissue)) +
  geom_smooth(data=SampleCovSummary, aes(x=covFrac, y=projectedCovFrac)) +
  ggtitle("Projected coverage of samples based on ultra-low sequencing") + 
  xlim(0,0.6) + ylim(0,1)


SampleCovSummary %>% ggplot(aes(x=covFrac, y=projectedCovFrac, colour=Tissue)) + 
  geom_point() + ggtitle("Projected coverage of samples based on ultra-low sequencing") + xlim(0,0.6) + ylim(0,1)
SampleCovSummary %>% ggplot(aes(x=CovRatio, y=projectedCovFrac, colour=Tissue)) + geom_point(aes(colour=Tissue)) + xlim(0,50)
SampleCovSummary %>% ggplot(aes(x=CovRatio, y=covFrac)) + geom_point(aes(colour=Tissue)) + xlim(0,50) + ylim(0,0.75) +
  ggtitle("Coverage achieved at ~1x as a function of coverage ratio")  +
  geom_smooth()
  geom_line(data=new_x, aes(new_x$x, y=predict(m, newdata=new_x$x)), colour='red')

SampleCovSummary %>% ggplot(aes(x=MeanCov, y=covFrac, colour=Tissue)) + geom_point()
SampleCovSummary %>% ggplot(aes(x=MeanCov, y=projectedCovFrac, colour=Tissue)) + geom_point()
#SampleCovSummary %>% ggplot(aes(x=alpha, y=projectedCovFrac, colour=Tissue)) + geom_point()
#SampleCovSummary %>% ggplot(aes(x=covFrac, y=CovRatio, colour=Tissue)) + geom_point()
#SampleCovSummary %>% ggplot(aes(x=MeanCov, y=projectedCovFrac, colour=Tissue)) + geom_point()
#SampleCovSummary %>% ggplot(aes(x=MeanCov, y=covFrac, colour=Tissue)) + geom_point()


#### BUILDING A MODEL OF THE CovRatio vs covFrac
library(splines)
library(modelr)
model_matrix(SampleCovSummary, covFrac ~ ns(CovRatio, 2))
m<-nls(data=SampleCovSummary, covFrac ~ a*CovRatio/(b+CovRatio),start=list(a=1,b=1))
#m<-nls(data=SampleCovSummary, covFrac ~ a + b * I(CovRatio^z), start=list(a=0.1,b=0.1,z=1))
cor(SampleCovSummary$covFrac, predict(m))
new_x = data.frame(x=seq(from=1, to=50, length.out=27))

# (2) Percent coverage by merged tissue
dfH <- subset(df, Tissue=="HE" )
dfB <- subset(df, Tissue=="BE" )
dfA <- subset(df, Tissue=="AE" )
dfC <- subset(df, Tissue=="C" )

sum(as.numeric(dfH$NBasesCovered))

dfH %>%  ggplot( aes(x=RelSize, colour=Covered)) + geom_density() + xlim(0,0.0005)
