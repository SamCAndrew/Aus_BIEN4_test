
# NOTE: Users don't need to run this; i already ran it to create the environmental layers you're using

#====================================================================
# prep environmental layers (some of these are made below)
# get climate
cl.f=list.files('/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Patrick_Test_1_8/Env3',full.names=T)
cl=stack(cl.f)
myprj = "+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs=0,0,0"
projection(cl)=myprj
# get soil
so.f=list.files('/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soils_10km_clean',full.names=T)
so=stack(so.f)
so=projectRaster(so,cl)
out=cropcrop(cl,so)
env=stack(out[[1]],out[[2]])
names(env)=c('aridity','bio1','bio12','bio15','bio2','bio20','bedrockDepth','bulkDensity1_4','clay1_4','ph1_4','silt1_4')
writeRaster(env,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Env_1_21_18/bien_v1_env.grd',overwrite=T)

# env.f=paste0(allBaseDirs$envDir,'/bien_v1_env.grd')
# if(!file.exists(env.f)){
#   original.files=list.files(allBaseDirs$envDir,'.tif',full.names=T)
#   env=stack(original.files)
#   projection(env)=myprj
#   #bb=crop(env,c(-5e6,5e6,-7.5e6,3e6))
#   bb=env
#   writeRaster(bb,file=paste0(allBaseDirs$envDir,'/bio_sp8.tif'),overwrite=T)
#   file.remove(original.files)
# }
# env=stack(env.f)
# print('loaded env raster')


#====================================================================
#====================================================================
# prep future env
	# convert to grd
	# project and crop
	# stack with present soil
outDir='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/BIEN4_NCEAS_Test/Future_Env'
pres.env=stack('/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Env_1_21_18/bien_v1_env.grd')
soil=pres.env[[c('bedrockDepth','bulkDensity1_4','clay1_4','ph1_4','silt1_4')]]
fut.f=list.files('/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Future_Env_Raw_3_1',full.names=T)
for(i in seq_along(fut.f)){
	f.r=stack(fut.f[i])
	f.r1=f.r[[c(6,1,3,4,2,5)]]
	names(f.r1)=c('aridity','bio1','bio12','bio15','bio2','bio20')
	tmp=cropcrop(f.r1,soil)
	f.all=stack(tmp[[1]],tmp[[2]])
	writeRaster(f.all,file= paste0(outDir,'/',tools::file_path_sans_ext(basename(fut.f[i])),'.grd'),overwrite=T)
}


#====================================================================
#====================================================================
#====================================================================
#====================================================================
# TINKERING/EXPLORATION BELOW


#====================================================================
# organize the other env data

myOtherEnvDir='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/OtherEnv'

e50=stack(paste0(myOtherEnvDir,'/cc8550.tif'))

a50=raster('/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/BIENWorkflow1-0_tests/accum_aridity_cc8550_10m.tif')
a50.1=projectRaster(a50,e50)

e=stack(e50,a50.1)

writeRaster(e,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/BIENWorkflow1-0_tests/Env3_OtherEnv/cc8550.tif')

e70=stack(paste0(myOtherEnvDir,'/cc8570.tif'))

a70=raster('/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/BIENWorkflow1-0_tests/accum_aridity_cc8570_10m.tif')
a70.1=projectRaster(a70,e70)
e1=stack(e70,a70.1)

writeRaster(e1,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/BIENWorkflow1-0_tests/Env3_OtherEnv/cc8570.tif')


# read original 

e2=stack(list.files('/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/BIENWorkflow1-0_tests/Env3',full.names=T))
e2.1=e2[[c(2,4,3,5,6,1)]]
writeRaster(e2.1,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/BIENWorkflow1-0_tests/Env3_clean/present.tif')

plot(e1)
plot(e2.1)
plot(e)

#----------------------------------------------------------------------
# explore soils
s.f=list.files('/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soils',full.names=T)
s.r=stack(s.f)

aa=raster('/Users/ctg/Documents/Worlclim/V2/wc2.0_10m_prec/wc2.0_10m_prec_01.tif')
plot(aa)
aa

be=stack('/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/BIENWorkflow1-0_tests/Env3_clean/present.tif')
projection(be)="+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs=0,0,0"

aa=raster('/Users/ctg/Documents/Data/Soils/BDRICM_M_250m_ll.tif')
aa1=crop(aa,c(-180,-30,-60,80))
aa2=raster::aggregate(aa1,fact=3)
writeRaster(aa2,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids/BDRICM_M_ll.tif')
aa3=projectRaster(aa2,crs="+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs=0,0,0")
writeRaster(aa3,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids_reproj/BDRICM_M_ll.tif')
gc()
aa4=cropcrop(aa3,be)
aa5=raster::resample(aa4[[1]],be)
writeRaster(aa5,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids_Clean/BDRICM_M_ll.tif')


aa=raster('/Users/ctg/Documents/Data/Soils/PHIHOX_M_sl3_250m_ll.tif')
aa1=crop(aa,c(-180,-30,-60,80))
aa2=raster::aggregate(aa1,fact=3)
writeRaster(aa2,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids/PHIHOX_M_sl3_ll.tif')
aa3=projectRaster(aa2,crs="+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs=0,0,0")
writeRaster(aa3,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids_reproj/PHIHOX_M_sl3_ll.tif',overwrite=T)
gc()
aa4=cropcrop(aa3,be)
aa5=raster::resample(aa4[[1]],be)
writeRaster(aa5,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids_Clean/PHIHOX_M_sl3_ll.tif')


aa=raster('/Users/ctg/Documents/Data/Soils/SLTPPT_M_sl3_250m_ll.tif')
aa1=crop(aa,c(-180,-30,-60,80))
aa2=raster::aggregate(aa1,fact=3)
writeRaster(aa2,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids/SLTPPT_M_sl3_ll.tif')
aa3=projectRaster(aa2,crs="+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs=0,0,0")
writeRaster(aa3,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids_reproj/SLTPPT_M_sl3_ll.tif')
gc()
aa4=cropcrop(aa3,be)
aa5=raster::resample(aa4[[1]],be)
writeRaster(aa5,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids_Clean/SLTPPT_M_sl3_ll.tif')


aa=raster('/Users/ctg/Documents/Data/Soils/TAXOUSDA_250m_ll.tif')
aa1=crop(aa,c(-180,-30,-60,80))
aa2=raster::aggregate(aa1,fact=3)
writeRaster(aa2,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids/TAXOUSDA_ll.tif')
aa3=projectRaster(aa2[[1]],crs="+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs=0,0,0")
writeRaster(aa3,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids_reproj/TAXOUSDA_ll.tif')
gc()
aa4=cropcrop(aa3,be)
aa5=raster::resample(aa4[[1]],be)
writeRaster(aa5,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids_Clean/TAXOUSDA_ll.tif')

aa=raster('/Users/ctg/Documents/Data/Soils/af_CLYPPT_T__M_sd3_250m.tif')
aa1=crop(aa,c(-180,-30,-60,80))
aa2=raster::aggregate(aa1,fact=3)
writeRaster(aa2,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids/af_CLYPPT_T__M_sd3.tif')
aa3=projectRaster(aa2[[1]],crs="+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs=0,0,0")
writeRaster(aa3,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids_reproj/af_CLYPPT_T__M_sd3.tif')
gc()
aa4=cropcrop(aa3,be)
aa5=raster::resample(aa4[[1]],be)
writeRaster(aa5,file='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids_Clean/af_CLYPPT_T__M_sd3.tif')

pdf('~/Desktop/soil.pdf',w=8,h=8)
soil.r=stack(list.files('/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soil_Grids_Clean/',full.names=T))
names(soil.r)=c('depth.to.bedrock','soilPh.at..15','silt.at..15','class')
plot(soil.r)
dev.off()

#==========================================================================
#== use 10km layers
figDir='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Figures'
soilDir='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soils_10km'

# plot all soils
soil.f=list.files(soilDir,full.names=T)
soil.r=stack(soil.f)
soil.r=crop(soil.r,c(-180,-30,-60,80))
pdf(paste0(figDir,'/soil_all.pdf'),w=8,h=11)
par(mfrow=c(3,3))
for(i in 1:nlayers(soil.r)) plot(soil.r[[i]],main=names(soil.r)[i])
dev.off()

# write out just 5 summaries
outDir='/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soils_10km_clean'

name=c('BLDFIE','CLYPPT','PHIHOX','SLTPPT')
for( i in seq_along(name)){
	soil.f=list.files(soilDir,full.names=T,pattern=name[i])[1:4]
	soil.r=stack(soil.f)
	soil.r=crop(soil.r,c(-180,-30,-60,80))
	calc(soil.r,mean,file=paste0(outDir,'/',name[i],'1_4_10km.tif'),overwrite=T)
}

bedrock=raster('/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/Soils_10km/BDRICM_M_10km_ll.tif')
bedrock=crop(bedrock,c(-180,-30,-60,80))
writeRaster(bedrock,paste0(outDir,'/BDRICM_M_10km.tif'),overwrite=T)

# check correlations among soils
s.r=stack(list.files(outDir,full.names=T))
cors=cor(values(s.r),use='complete.obs')

fig.f=paste0(figDir,'/soil_layers_10km.pdf')
pdf(fig.f,w=10,h=10)
	par(mfrow=c(2,3))
	for(i in 1:nlayers(s.r)) plot(s.r[[i]],main=names(s.r)[i])
	corrplot::corrplot(cors,order = "AOE", addCoef.col = "grey",number.cex=.6)
dev.off()
system(paste0('open ',fig.f))

# check correlations among soils and climate
cl=stack('/Users/ctg/Dropbox/Projects/BIEN/BIENWorkflow_Input/BIENWorkflow1-0_tests/Env3_clean/present.tif')
myprj = "+proj=laea +lat_0=15 +lon_0=-80 +x_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs=0,0,0"
projection(cl)=myprj
so=projectRaster(s.r,cl)
out=cropcrop(cl,so)
env=stack(out[[1]],out[[2]])

cors=cor(values(env),use='complete.obs')
fig.f=paste0(figDir,'/soil_clim_cor.pdf')
pdf(fig.f,w=6,h=6)
	corrplot::corrplot(cors,order = "AOE", addCoef.col = "grey",number.cex=.6)
dev.off()
system(paste0('open ',fig.f))

