######
##data prep and model fitting for Slingsby et al 2020
######

####################
##setup
####################
renv::init()

### Load libraries
libs=c(
  "doParallel",
  "rasterVis",
  "rgdal",
  "reshape2",
  "sp",
  "knitr",
  "rmarkdown",
  "ggplot2",
  "tidyr",
  "dplyr",
  "minpack.lm",
  "maptools",
  "lubridate",
  "rjags",
  "dclone",
  "raster")
lapply(libs, require, character.only=T)

#file locations and names
mdatwd <- "data/"
mname <- "peninsulaDec2019" #model name for file naming

# Calculate the number of cores
no_cores <- detectCores() - 1

# # Initiate cluster
cl <- makeCluster(no_cores, type = "FORK")
registerDoParallel(cl)

###########################################################
###Load data
###########################################################

#spatial covariates and vegetation age
load(paste0(mdatwd,"inputData.RData"))
#modis ndvi and quality flag from MODIS/006/MYD13Q1 and MODIS/006/MOD13Q1
NDVI <- stack(paste0(mdatwd,"NDVI"))
QA <- stack(paste0(mdatwd,"QA"))

#mask input data to where we have fire data
rfi <- mask(rfi, vegtype)
NDVI <- mask(NDVI, rfi[[1]]) #note use of rfi - already cropped by LC
QA <- mask(QA, rfi[[1]])


###########################################################
###Convert NDVI, QA and age data to df
###########################################################

#NDVI
dat <- as.data.frame(rasterToPoints(NDVI))
datXY <- dat[,1:2]
dat <- as.data.frame(t(dat[,3:ncol(dat)])) #transpose data.frame
dat$Date <- as.Date(rownames(dat), format = "X%Y.%m.%d")
dat <- melt(dat, id=c("Date"))
datXY$variable <- paste0("V", rownames(datXY))
dat <- merge(dat, datXY)
names(dat)[which(names(dat)=="value")] <- "NDVI"

#QA data
qdat <- as.data.frame(rasterToPoints(QA))
qdatXY <- qdat[,1:2]
qdat <- as.data.frame(t(qdat[,3:ncol(qdat)])) #transpose data.frame
qdat$Date <- as.Date(rownames(qdat), format = "X%Y.%m.%d")
qdat <- melt(qdat, id=c("Date"))
qdatXY$variable <- paste0("V", rownames(qdatXY))
qdat <- merge(qdat, qdatXY)
names(qdat)[which(names(qdat)=="value")] <- "QA"

#Age (rfi) data
adat <- as.data.frame(rasterToPoints(rfi))
adatXY <- adat[,1:2]
adat <- as.data.frame(t(adat[,3:ncol(adat)])) #transpose data.frame
adat$Date <- as.Date(rownames(adat), format = "X%Y.%m.%d")
adat <- melt(adat, id=c("Date"))
adatXY$variable <- paste0("V", rownames(adatXY))
adat <- merge(adat, adatXY)
names(adat)[which(names(adat)=="value")] <- "Age"

#Merge the three
cdat <- merge(qdat, dat, by = c("Date", "x", "y"))
cdat <- merge(cdat, adat, by = c("Date", "x", "y"))
cdat <- cdat[,c("Date", "x", "y", "Age", "NDVI", "QA")] #drop unwanted columns
cdat$UI <- paste(cdat$x, cdat$y, sep = "_") #add unique identifier not rounded
cdat$UIJ <- paste(round(cdat$x, 3), round(cdat$y, 3), sep = "_") #add unique identifier rounded to coords with 3 decimal places (that matches covariates)

###########################################################
###Filter by a few NB criteria and trim covariates to match
###########################################################


#drop NDVI values < 0 or bad quality (QA !=0) and pixels (UI) with =<300 data points
cdat <- cdat %>% filter(NDVI > 0) %>% group_by(UIJ) %>% filter(QA == 0) %>% filter(n() >= 300)

#trim covariates and temporal data to match
cov$UIJ <- paste(round(cov$x, 3), round(cov$y, 3), sep = "_") #replace unique identifier rounded to coords with 3 decimal places (that matches temporal data)
cov$UI <- paste(cov$x, cov$y, sep = "_") #unique identifier not rounded

cov <- cov %>% filter(complete.cases(cov)) %>% 
  filter(UIJ %in% cdat$UIJ)
cdat <- filter(cdat, UIJ %in% cov$UIJ)

#drop unwanted levels
cov <- droplevels(cov)
cdat <- droplevels(cdat)

###########################################################
###Quick data check and clear unwanted objects from memory
###########################################################

rm(list = c("NDVI", "QA", "rfi", "adat", "adatXY", "dat", "datXY", "libs", "qdat", "qdatXY")); gc()

###########################################################
###Add columns for month of fire and age in years
###########################################################

#create new variable indicating month of fire
cdat <- cdat %>% mutate(firemonth=month(Date-Age))

#create age column in years
cdat$DA <- cdat$Age/365.25

###########################################################
###Have a look at the raw data
###########################################################

# P <- ggplot(data = cdat, aes(x = DA, y = NDVI, colour = UI)) +
#     geom_line()
# #+ geom_vline(aes(xintercept = as.numeric(Fire)))


###########################################################
### Create dummy variables for veg type and select and scale environmental data
###########################################################

#create dummy vars for veg type
cov$vegnut <- recode_factor(cov$vegtype, "Granite Fynbos" = "High", "Sand Fynbos" = "High", "Sandstone Fynbos" = "Low", "Shale Fynbos" = "High", "Shale Renosterveld" = "High", "Western Strandveld"  = "High")
cov$vegnut <- as.factor(cov$vegnut)
tveg <- as.numeric(cov$vegnut) - 1
dummies <- model.matrix(~as.factor(tveg))
dummies <- dummies[,-1]

#select vars and scale data
envars <- c("slope", "aspect", "tpi", "prec1", "prec7", "tmax1", "tmin7") #"dem", 
scaled <- scale(as.matrix(cov[,envars]))
env <- as.data.frame(cbind(intercept=1, scaled, tveg))
env$UI <- cov$UI
env$UIJ <- cov$UIJ
env <- cbind(env,dummies)

#have a look at collinearity
# require(GGally)
# ggpairs(as.data.frame(env_full))

#save the scaling parameters to convert fitted coefficients back to metric units later
beta.mu=c(intercept=0,attr(scaled,"scaled:center")) #,rep(0,ncol(tveg)))
beta.sd=c(intercept=1,attr(scaled,"scaled:scale")) #,rep(1,ncol(tveg)))
rm(scaled)  #drop the scaled data

###########################################################
###Format data for JAGS
###########################################################

tdat <- cdat

# #set up cross-validation dataset?
# holdout <- 0.00
# set.seed(111)
# s <- sort(sample(unique(cov$UI),round(length(unique(cov$UI))*(1-holdout)))); length(s)
# #NA sites in cv set
# tdat$NDIN=tdat$NDVI
# tdat$NDIN[!(tdat$UI%in%s)]=NA; gc() 

tdat$NDIN <- tdat$NDVI
#if we want to predict NDVI for date beyond 2014-05-31
#tdat$NDIN[(tdat$Date>as.Date("2014-05-31"))]=NA; gc() 

#other wise only send data up to 2014-05-31 to jags
tdat <- tdat[(tdat$Date<=as.Date("2014-05-31")),]; gc() 

#create new id that goes from 1 to nGrid (to order env and tdat in the same way)
env$jag_id <- as.integer(as.factor(env$UIJ))
jtab <- data.frame(UIJ=env$UIJ,jag_id=env$jag_id, stringsAsFactors=F)
tdat <- left_join(tdat, jtab, by='UIJ')

#arrange temporal and env data into same order
drop.cols <- c('jag_id', 'UI','UIJ')
env <- env[order(env$jag_id),] 
save(env,file=paste(mdatwd,mname,"_envdata.Rdata",sep="")) #save env for analysing results
env <- env %>% dplyr::select(-one_of(drop.cols))
env <- as.matrix(env)
tdat <- tdat[order(tdat$jag_id),]

#final check
if(length(unique(tdat$jag_id)) != nrow(env))  print("sites not matching between spatial and temporal data!")

#for alter analysis
save(tdat,file=paste(mdatwd,mname,"_inputdata_small.Rdata",sep="")) #save env for analysing results
###########################################################
###Prep JAGS inputs
###########################################################

#get counts
nGrid=length(unique(tdat$jag_id))       ;nGrid
nTime=length(unique(tdat$Date))          ;nTime
nBeta=ncol(env)                          ;nBeta

#write data object
data=list(
  age=tdat$DA,
  ndvi=tdat$NDIN,
  id=tdat$jag_id,
  firemonth=tdat$firemonth,
  nObs=nrow(tdat),
  env=env,
  nGrid=nGrid,
  nBeta=nBeta
)

#function to generate initial values
gen.inits=function(nGrid,nBeta) { list(
  ## spatial terms
  alpha=runif(nGrid,0.1,0.5),
  gamma=runif(nGrid,0.1,.9),
  A=runif(nGrid,0.1,.9),
  lambda=runif(nGrid,0.2,1),
  ## spatial means
  alpha.mu=runif(1,0.1,0.2),
  ## priors  
  gamma.beta=runif(nBeta,0,1),
  gamma.tau=runif(1,1,5),
  alpha.tau=runif(1,1,5),
  lambda.beta=runif(nBeta,0,2),
  lambda.tau=runif(1,0,2),
  A.beta=runif(nBeta,0,1),
  A.tau=runif(1,1,5),
  tau=runif(1,0,2)
)
}

#list of parameters to monitor (save)
params=c("phi","gamma.beta","gamma.sigma","A.beta","A.sigma","alpha","gamma","lambda","A",
         "alpha.mu","alpha.sigma","lambda.beta","lambda.sigma","sigma")

###########################################################
###Save all data into Rdata object for model fitting
###########################################################

save.image(file=paste(mdatwd,mname,"_inputdata.Rdata",sep="")) 

rm(list = ls()[-which(ls() %in% c("mdat", "mname", "data", "params", "cl", "mdatwd", "gen.inits"))])
gc()

###########################################################
###Run JAGS
###########################################################

foutput=paste0(mdatwd, mname, "_modeloutput.Rdata")


m <- jags.parfit(cl = cl, #runs chains in parallel with library(dclone)
                 data = data, 
                 params = params, 
                 model = paste0(mdatwd,"Model.R"), 
                 inits = gen.inits(data$nGrid,data$nBeta), 
                 n.chains = 3,
                 n.adapt=1000,n.update=1000,
                 thin = 2, n.iter = 2000
)


save(m,file=foutput)

