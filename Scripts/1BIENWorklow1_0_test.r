#== 1. packages
library(cmsdm)
library(BIENWorkflow)
library(trinaryMaps)
## Make the custom background files automatically copy in sdmdirectories
myInputDir="D:/Documents_2/BIEN4_NCEAS_Test"
myModelingDir="D:/Documents_2/BIEN4_NCEAS_Test/BIEN_v1"
if(!file.exists(myModelingDir)) dir.create(myModelingDir)
#====================================================================
#== 2. Create directories
#====================================================================

runName='sam_demo_180321'

# system specific
baseDir=paste0(myModelingDir,'/', runName,'_inputs')
mySpeciesDir=paste0(myInputDir,'/Pres_demo')
# mySpeciesDir=paste0(myInputDir,'/Pres') # for all species
myCustomBackgroundDir=paste0(myInputDir,'/BIEN_Family_Pres')
myCustomBackgroundTable= system.file("extdata/OFGS_Lookup_Table.csv", package='BIENWorkflow')
myEnvDir=paste0(myInputDir,'/Env_1_21_18')
# each of these should be stored as a single stack
myOtherEnvDir=paste0(myInputDir,'/Future_Env') 
#myOtherEnvDir=paste0(myInputDir,'/Future_Env_Demo') # for testing

# myOtherEnvDir=list(system.file("extdata/OtherEnv/", package='BIENWorkflow'))
# for using offsets; not yet done in BIEN
myRangePolyDir=NULL
mySamplingModelDir=NULL
mySamplingDataDir=NULL
myOffsetDirs=NULL

# system independent
otherDirs=list(
  rangePolyDir=NULL,
  rdistDir=paste0(baseDir,'/rDists'),
  domains=paste0(baseDir,'/domains'))
sortDirNames=c('Points','Bounding_Box','Convex_Hull','SDM','SDM10_20')
offsetDataDirNames=c('Expert_Maps')
samplingDataDirNames=NULL

allBaseDirs=sdmDirectories(baseDir,
                           warn=FALSE,
                           mySpeciesDir=mySpeciesDir,
                           myEnvDir=myEnvDir,
                           myOffsetDirs=myOffsetDirs,
                           mySamplingDataDir=mySamplingDataDir,
                           mySamplingModelDir=mySamplingModelDir,
                           myOtherEnvDir=myOtherEnvDir,
                           myCustomBackgroundDir=myCustomBackgroundDir,
                           myCustomBackgroundTable=myCustomBackgroundTable,
                           offsetDataDirNames=offsetDataDirNames,
                           samplingDataDirNames=samplingDataDirNames,
                           overwrite=FALSE,
                           sortDirNames=sortDirNames,
                           otherDirs=otherDirs)

#====================================================================
#== 3. Prep environmental layers
#====================================================================
source(paste0(myInputDir,'/2BIENWorklow1_0_test_Make_Data.r'))
save(allBaseDirs,file=paste0(baseDir,'/allBaseDirs.rdata'))

#====================================================================
#== 4. Sort species
#====================================================================
#--------------------------------------------------------------------
### 4.a sort species by algorithm
#--------------------------------------------------------------------
sortDone=checkSortDone(allBaseDirs, deleteSorted=FALSE)
pointsProj=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
sortSpeciesBulk(allBaseDirs,pointsProj,myprj=myprj,doThin=TRUE,thinCutoff=20,verbose=TRUE,overwrite=FALSE,doClean=TRUE)

#--------------------------------------------------------------------
### 4.b specify which species to run
#--------------------------------------------------------------------
speciesList=list.files(file.path(allBaseDirs$speciesDir,'SDM'),full.names=T)
done=checkSDMDone(allBaseDirs,speciesList) # check whether any already done
speciesList=done$notRun
str(done)

#====================================================================
#== 5. Model settings
#====================================================================
modelSettings=list()
modelSettings$samplingSettings=c('noBias','targetBG') 
modelSettings$expertSettings=list(prob=1e-6,rate=0,skew=1e-6,shift=0)
modelSettings$predictorSettings=NULL
modelSettings$regularizationSettings=NULL
modelSettings$priorSettings=NULL
modelSettings$backgroundSettings=NULL
modelSettings$algorithmSettings=c('maxnet')
modelSettings$formulaMaker="sdmGenericFormulas"
modelSettings$maxTime=500 # time for an individual call to glmnet()
modelSettings$samplingModel='none'

source(system.file("extdata/Demo_Workflow/Common_BIENWorkflow_Settings_v1.r", package='BIENWorkflow'))
if(Sys.info()['user']=='ctg') {
	source( '~/Dropbox/Projects/BIEN/BIENWorkflow/inst/extdata/Demo_Workflow/Common_BIENWorkflow_Settings_v1.r')
}

toDo$maxBGFit=3e4
toDo$otherEnvProjections=TRUE
toDo$plotOtherEnvPred=TRUE
toDo$makeBiasedBackground=TRUE
toDo$openFig=TRUE
toDo$whichFutures='all' #c('he') # put anything here you could grep from the scenario name
toDo$whichFuturesToPlot=c('he')
toDo$plotBinaryShp=F
toSave$shapeFile=F
toDo$maxTrimDomainTime=60
# this could give you trouble; you might locate it with Sys.which('gdal_polygonize.py') or set toSave$shapefile=F; toDo$$plotBinaryShp=F
toDo$pypath="/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdal_polygonize.py"
#on PC you probably should set toDo$evalFork=F
toDo$evalFork=F
#--------------------------------------------------------------------
### 6. Run Models
#--------------------------------------------------------------------
j=1; (speciesCSV=speciesList[j])

for (j in 1:length(speciesList)){
  print(j)
  out=sdmBIEN(speciesCSV=speciesList[j],
              allBaseDirs,
              modelSettings=modelSettings,
              toDo=toDo,
              toSave=toSave,
              toOverwrite=toOverwrite)
  gc()
}
