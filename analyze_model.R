#################################
######Analysis of results and creation of rasters for upload to GEE
#################################
#only run if you have not before
#renv::init()

libs=c(
  "dplyr",
  "tidyr",
  "raster",
  "coda",
  "rjags")
lapply(libs, require, character.only=T)



#file locations and names
mdatwd <- "data/"
mname <- "peninsulaDec2019" #model name for file naming

#download the results if you did not create them in fit_model.R:
download.file('https://storage.googleapis.com/data-sharing-gmoncrieff/peninsulaDec2019_modeloutput.Rdata', destfile = paste0(mdatwd, mname, "_modeloutput.Rdata"))
download.file('https://storage.googleapis.com/data-sharing-gmoncrieff/peninsulaDec2019_envdata.Rdata', destfile = paste0(mdatwd, mname, "_envdata.Rdata"))
download.file('https://storage.googleapis.com/data-sharing-gmoncrieff/peninsulaDec2019_inputdata_small.Rdata', destfile = paste0(mdatwd, mname, "_inputdata_small.Rdata"))

#load results
foutput <- paste0(mdatwd, mname, "_modeloutput.Rdata")
envdata <- paste0(mdatwd,mname,"_envdata.Rdata")

#load model results
#load env data
load(foutput)
load(envdata)

#check model fitting - select a parameter
par <- m[,"phi"]
plot(par)

#extract parameter summaries
res <- as.data.frame(summary(m)$statistics)

#summarise non spatial pars
phi <- res %>%
  mutate(pname = rownames(res)) %>%
  filter(pname =='phi')

sig <- res %>%
  mutate(pname = rownames(res)) %>%
  filter(pname =='sigma')

#summarise spatial pars
alphas <- res %>%
  mutate(par = rownames(res)) %>%
  filter(str_detect(par, "alpha")) %>%
  filter(str_detect(par, "\\[")) %>%
  mutate(parstr = str_sub(par,1,5)) %>%
  mutate(parnum = str_remove(par,"alpha")) %>%
  mutate(parnum = str_sub(parnum,2,-2)) %>%
  select(c("Mean","SD","parstr","parnum"))

lambdas <- res %>%
  mutate(par = rownames(res)) %>%
  filter(str_detect(par, "lambda")) %>%
  filter(str_detect(par, "\\[")) %>%
  filter(str_detect(par, "beta",negate=TRUE)) %>%
  mutate(parstr = str_sub(par,1,6)) %>%
  mutate(parnum = str_remove(par,"lambda")) %>%
  mutate(parnum = str_sub(parnum,2,-2)) %>%
  select(c("Mean","SD","parstr","parnum"))

gammas <- res %>%
  mutate(par = rownames(res)) %>%
  filter(str_detect(par, "gamma")) %>%
  filter(str_detect(par, "\\[")) %>%
  filter(str_detect(par, "beta",negate=TRUE)) %>%
  mutate(parstr = str_sub(par,1,5)) %>%
  mutate(parnum = str_remove(par,"gamma")) %>%
  mutate(parnum = str_sub(parnum,2,-2)) %>%
  select(c("Mean","SD","parstr","parnum"))

As <- res %>%
  mutate(par = rownames(res)) %>%
  filter(str_detect(par, "A")) %>%
  filter(str_detect(par, "\\[")) %>%
  filter(str_detect(par, "beta",negate=TRUE)) %>%
  mutate(parstr = str_sub(par,1,2)) %>%
  mutate(parnum = str_remove(par,"A")) %>%
  mutate(parnum = str_sub(parnum,2,-2)) %>%
  select(c("Mean","SD","parstr","parnum"))

######################################
####convert to rasters
#####################################

env = env %>%
  separate(UIJ,c("lon","lat"),sep="_")

#rasterise
gammaM <- data.frame(lon = env$lon, lat = env$lat, gammas$Mean)
gammaSD <- data.frame(lon = env$lon, lat = env$lat, gammas$SD)
alphaM <- data.frame(lon = env$lon, lat = env$lat, alphas$Mean)
alphaSD <- data.frame(lon = env$lon, lat = env$lat, alphas$SD)
lambdaM <- data.frame(lon = env$lon, lat = env$lat, lambdas$Mean)
lambdaSD <- data.frame(lon = env$lon, lat = env$lat, lambdas$SD)
AM <- data.frame(lon = env$lon, lat = env$lat, As$Mean)
ASD <- data.frame(lon = env$lon, lat = env$lat, As$SD)

gammaMras<-rasterFromXYZ(gammaM,res=c(0.002607436,0.002607436),crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',digits=0.3)
gammaSDras<-rasterFromXYZ(gammaSD,res=c(0.002607436,0.002607436),crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',digits=0.3)
alphaMras<-rasterFromXYZ(alphaM,res=c(0.002607436,0.002607436),crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',digits=0.3)
alphaSDras<-rasterFromXYZ(alphaSD,res=c(0.002607436,0.002607436),crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',digits=0.3)
lambdaMras<-rasterFromXYZ(lambdaM,res=c(0.002607436,0.002607436),crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',digits=0.3)
lambdaSDras<-rasterFromXYZ(lambdaSD,res=c(0.002607436,0.002607436),crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',digits=0.3)
AMras<-rasterFromXYZ(AM,res=c(0.002607436,0.002607436),crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',digits=0.3)
ASDras<-rasterFromXYZ(ASD,res=c(0.002607436,0.002607436),crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',digits=0.3)

#write to disk
writeRaster(gammaMras,paste0(mdatwd, "/EarthEngine/gammaM.tif"))
writeRaster(gammaSDras,paste0(mdatwd, "/EarthEngine/gammaSD.tif"))
writeRaster(lambdaMras,paste0(mdatwd, "/EarthEngine/lambdaM.tif"))
writeRaster(lambdaSDras,paste0(mdatwd, "/EarthEngine/lambdaSD.tif"))
writeRaster(alphaMras,paste0(mdatwd, "/EarthEngine/alphaM.tif"))
writeRaster(alphaSDras,paste0(mdatwd, "/EarthEngine/alphaSD.tif"))
writeRaster(AMras,paste0(mdatwd, "/EarthEngine/AM.tif"))
writeRaster(ASDras,paste0(mdatwd, "/EarthEngine/ASD.tif"))
