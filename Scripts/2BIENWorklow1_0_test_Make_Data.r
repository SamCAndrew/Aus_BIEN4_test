
#!! NOTE: Users shouldn't need to modify this file; it's sourced from file 1

#====================================================================
#====================================================================
# Env
#====================================================================
#====================================================================

# !!!NOTE: this may change from workflow to workflow (maybe different file type)
#myprj = "+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs=0,0,0" for New world
myprj<-"+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#env.f=paste0(allBaseDirs$envDir,'/bio_sp8.tif')
env.f=paste0(allBaseDirs$envDir,'/Aus_bien_v1_env.grd')
if(!file.exists(env.f)){
	print('make your env file, dummy!')
#   original.files=list.files(allBaseDirs$envDir,'.tif',full.names=T)
#   env=stack(original.files)
#   projection(env)=myprj
#   #bb=crop(env,c(-5e6,5e6,-7.5e6,3e6))
#   bb=env
#   writeRaster(bb,file=paste0(allBaseDirs$envDir,'/bio_sp8.tif'),overwrite=T)
#   file.remove(original.files)
}
env=stack(env.f)
print('loaded env raster')

## Other Env layers
#if(length(list.files(allBaseDirs$otherEnvDir))==0){
#  original.files=list.files(myOtherEnvDir,full.names=T)
#  file.copy(original.files,paste0(allBaseDirs$otherEnvDir,'/', basename(original.files)))
#   otherEnv=lapply(original.files, function(x){
#     name=file_path_sans_ext(basename(x))
#     tmp1=stack(x)
#     #tmp1=tmp1[[c(1,3,4,2,5,6)]] # 1,2,12,15,20, aridity
#     projection(tmp1)=myprj # Will have to be sure this is true before hand!!
#     tmp2=crop(tmp1,env)
#     writeRaster(tmp2,file=paste0(allBaseDirs$otherEnvDir,'/',name,'.tif'), overwrite=T)
#   })
#}
#print('prepped other env rasters')

## shapefile for plotting
world.shp=readOGR(system.file("extdata/TM_WORLD_BORDERS_SIMPL-0.3/TM_WORLD_BORDERS_SIMPL-0.3.shp",package='BIENWorkflow'),'TM_WORLD_BORDERS_SIMPL-0.3',verbose=F)
#world.shp=crop(world.shp,c(-140,-30,-60,45))
world.shp=spTransform(world.shp,projection(env))
print('loaded world shapefile')


## aggregate env layer for subsampling very large species

coarse.grid.file=allBaseDirs$coarseGridFile
if(!file.exists(coarse.grid.file)){
	env.template=stack(env.f)[[1]]
	coarse.grid=aggregate(env.template,fact=2)
	# plot(coarse.grid)
	writeRaster(coarse.grid,file=coarse.grid.file,overwrite=T)
	print('wrote coarse grid file')
}


## raster template for buffering domains (used in trimDomain)
template.file=paste0(allBaseDirs$miscDir,'/templateRaster.tif')
if(!file.exists(template.file)){
	template=!is.na(env[[1]])
	template[template==0]=NA
	template[template==1]=0
	writeRaster(template,template.file,options = c("COMPRESS=LZW", "PREDICTOR=2"),datatype = "INT2S", overwrite = TRUE )
	# plot(template)
	print('wrote env template')
}
