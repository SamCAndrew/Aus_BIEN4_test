#library(colorout)
if(Sys.info()['user']=='saman') {
  # this installs my development versions; you can just install the package once
  devtools::install('Packages_3_12_18/cmsdm')
	devtools::install('Packages_3_12_18/BIENWorkflow')
	devtools::install('Packages_3_12_18/trinaryMaps')  
	#devtools::install( '/Users/ctg/Dropbox/Projects/BIEN/Private_RBIEN/todoBIEN')
}
library(cmsdm)
library(BIENWorkflow)
library(trinaryMaps)
#library(todoBIEN)
