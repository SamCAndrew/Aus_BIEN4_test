# prep full bien data for target group and genus data
# by cory on 6.16.17

library(dplyr)
library(data.table)
library(doParallel)
registerDoParallel(4)

bigfile='/Users/ctg/Documents/Data/nw_occurrences_for_cory_6_15_17.csv'
bigfile.sample <- read.csv(bigfile, stringsAsFactors=FALSE, header=T, nrows=20)
bigfile.colclass <- sapply(bigfile.sample,class)

bigfile.raw <- tbl_df(read.csv(bigfile,
                               stringsAsFactors=FALSE, header=T,
                               colClasses=bigfile.colclass, comment.char=""))
start.time=proc.time()
bigfile.raw=fread(bigfile)
proc.time-start.time
object.size(bigfile.raw)/1e6 


keep=!is.na(bigfile.raw$longitude)

bigfile.raw=bigfile.raw[keep,]
keep=!is.na(bigfile.raw$latitude)
bigfile.raw=data.frame(bigfile.raw[keep,])
write.csv(bigfile.raw,file='/Users/ctg/Dropbox/Projects/BIEN/Occurrence_Data/All_BIEN_6_15_17.csv',row.names=FALSE)

bigfile.raw=data.frame(fread('/Users/ctg/Dropbox/Projects/BIEN/Occurrence_Data/All_BIEN_6_15_17.csv'))
genus=unique(bigfile.raw[,c('scrubbed_family','scrubbed_genus')])
family=unique(bigfile.raw[,c('scrubbed_family')])

#=============================================================================#write out family files
outdir.family='/Users/ctg/Dropbox/Projects/BIEN/Occurrence_Data/BIEN_Family_Pres/'
if(!file.exists(outdir.family)) dir.create(outdir.family)

foreach(ii = 1:length(family)) %dopar% {
  out=subset(bigfile.raw,scrubbed_family==family[ii])[, c('longitude','latitude')]
  write.csv(out,file=paste0(outdir.family,family[ii],'.csv'),row.names=FALSE)
}

#=============================================================================#write out genus files
outdir.genus='/Users/ctg/Dropbox/Projects/BIEN/Occurrence_Data/BIEN_Genus_Pres/'
if(!file.exists(outdir.genus)) dir.create(outdir.genus)

foreach(ii = 1:nrow(genus)) %dopar% {
  out=subset(bigfile.raw,scrubbed_family==genus[ii,1] & scrubbed_genus==genus[ii,2] )[,c('longitude','latitude')]
  write.csv(out,file=paste0(outdir.genus,paste0(genus[ii,], collapse='_'),'.csv'),row.names=FALSE)
}

#=============================================================================# Make lookup table of family - genus - species
ofgs.all=bigfile.raw[,c(1:4,7)]
ofgs=unique(ofgs.all)
names(ofgs)=c("order","family","genus","species","growth_form")
ofgs$species=gsub(' ','_',ofgs$species)
write.csv(ofgs,file='/Users/ctg/Dropbox/Projects/BIEN/Occurrence_Data/OFGS_Lookup_Table.csv',row.names=F)
write.csv(ofgs,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow/inst/extdata/OFGS_Lookup_Table.csv',row.names=F)

#=============================================================================# check the size number of samples per family and genus to see if there are any that will lead to trouble with target groups