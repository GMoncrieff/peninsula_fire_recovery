###########################################################
###Script to generate plots for Slingsby, Moncriedff and Wilson 2020
###########################################################


#only run if you have not before
#renv::init()

libs=c("dplyr",
"tidyr",
"ggplot2",
"reshape2",
"raster",
"coda",
"rjags",
"tictoc",
"readxl",
"scales",
"sf",
"cowplot",
"stringr")
lapply(libs, require, character.only=T)

#file locations and names
mdatwd <- "data/"
mname <- "peninsulaDec2019" #model name for file naming

###########################################################
###Get exceedance rasters (GEE output) and plot
###NOTE: further plotting based on JAGS outputs below
###Hash out this section if you don't have GEE outputs yet
###########################################################

# Get data and wrangle for plotting
exceed <- stack("data/exceed_below.tif", "data/exceed_above.tif")
names(exceed) <- c("below", "above")
exceed <- projectRaster(exceed, crs = CRS("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +units=m +no_defs"))
edat <- as.data.frame(rasterToPoints(exceed, spatial = F))
edat <- melt(edat, id = c("x", "y"))
edat <- fortify(edat)

# Get a pretty coastline for plotting
coast <- st_read("Data/coastline") #, layer = "coastline")
coast <- fortify(coast)

# Plot

g <- ggplot() +
  geom_sf(data = coast, fill = "skyblue1") +
  geom_tile(data = edat, aes(x, y, fill = value)) + 
  scale_fill_gradient(low = "#FFF5F0", high = "#67000D", na.value = "transparent") +
  facet_wrap(~variable) +
  theme_void() +
  labs(fill="Deviance") + 
  annotate("rect", xmin = 2045607, xmax = 2054169, ymin = -4045661, ymax = -4040601, fill = "transparent", colour = "grey30") + 
  annotate("text", label = "Silvermine", x = 2050000, y = -4047000, colour = "grey30") +
  annotate("rect", xmin = 2037356, xmax = 2042970, ymin = -4037517, ymax = -4031284, fill = "transparent", colour = "grey30") + 
  annotate("text", label = "Karbonkelberg", x = 2040000, y = -4038500, colour = "grey30") +
  annotate("rect", xmin = 2050000, xmax = 2052500, ymin = -4072500, ymax = -4069250, fill = "transparent", colour = "grey30") + 
  annotate("text", label = "Cape of Good Hope", x = 2050000, y = -4074000, colour = "grey30") +
  annotate("rect", xmin = 2053500, xmax = 2056500, ymin = -4061000, ymax = -4058500, fill = "transparent", colour = "grey30") + 
  annotate("text", label = "Miller's Point", x = 2055500, y = -4062500, colour = "grey30")

ggsave(filename = "figures/exceedmap.png", plot = g, device = NULL, path = NULL, scale = 1, width = 18, height = 18, units = "cm", dpi = 300, limitsize = TRUE)

###########################################################
###Get model data and prep for model prediction and plotting
###########################################################

#download the results if you did not create them in fit_model.R:
download.file('https://storage.googleapis.com/data-sharing-gmoncrieff/peninsulaDec2019_modeloutput.Rdata', destfile = paste0(mdatwd, mname, "_modeloutput.Rdata"))
download.file('https://storage.googleapis.com/data-sharing-gmoncrieff/peninsulaDec2019_envdata.Rdata', destfile = paste0(mdatwd, mname, "_envdata.Rdata"))
download.file('https://storage.googleapis.com/data-sharing-gmoncrieff/peninsulaDec2019_inputdata_small.Rdata', destfile = paste0(mdatwd, mname, "_inputdata_small.Rdata"))

#load results
foutput <- paste0(mdatwd, mname, "_modeloutput.Rdata")
envdata <- paste0(mdatwd,mname,"_envdata.Rdata")
inputdata <- paste0(mdatwd,mname,"_inputdata_small.Rdata")

#load model results
#load env data
load(foutput)
load(envdata)
load(inputdata)

#add columns for results
tdat$mean <- NA
tdat$upper <- NA
tdat$lower <- NA
tdat$lq <- NA
tdat$uq <- NA

#all results
res <- as.data.frame(summary(m)$statistics)

#firemonth
phi <- res %>%
  mutate(pname = rownames(res)) %>%
  filter(pname =='phi')

#model sd
sig <- res %>%
  mutate(pname = rownames(res)) %>%
  filter(pname =='sigma')

#extract spatial pars
alphas <- res %>%
  mutate(par = rownames(res)) %>%
  filter(str_detect(par, "alpha")) %>%
  filter(str_detect(par, "\\[")) %>%
  mutate(parstr = str_sub(par,1,5)) %>%
  mutate(parnum = str_remove(par,"alpha")) %>%
  mutate(parnum = str_sub(parnum,2,-2)) %>%
  dplyr::select(c("Mean","SD","parstr","parnum"))

lambdas <- res %>%
  mutate(par = rownames(res)) %>%
  filter(str_detect(par, "lambda")) %>%
  filter(str_detect(par, "\\[")) %>%
  filter(str_detect(par, "beta",negate=TRUE)) %>%
  mutate(parstr = str_sub(par,1,6)) %>%
  mutate(parnum = str_remove(par,"lambda")) %>%
  mutate(parnum = str_sub(parnum,2,-2)) %>%
  dplyr::select(c("Mean","SD","parstr","parnum"))

gammas <- res %>%
  mutate(par = rownames(res)) %>%
  filter(str_detect(par, "gamma")) %>%
  filter(str_detect(par, "\\[")) %>%
  filter(str_detect(par, "beta",negate=TRUE)) %>%
  mutate(parstr = str_sub(par,1,5)) %>%
  mutate(parnum = str_remove(par,"gamma")) %>%
  mutate(parnum = str_sub(parnum,2,-2)) %>%
  dplyr::select(c("Mean","SD","parstr","parnum"))

As <- res %>%
  mutate(par = rownames(res)) %>%
  filter(str_detect(par, "A")) %>%
  filter(str_detect(par, "\\[")) %>%
  filter(str_detect(par, "beta",negate=TRUE)) %>%
  mutate(parstr = str_sub(par,1,1)) %>%
  mutate(parnum = str_remove(par,"A")) %>%
  mutate(parnum = str_sub(parnum,2,-2)) %>%
  dplyr::select(c("Mean","SD","parstr","parnum"))


#to get a raster of ids:
env = env %>%
  separate(UIJ,c("lon","lat"),sep="_")
jag_id <- data.frame(lon = env$lon, lat = env$lat, gammas$parnum)
jag_id_ras<-rasterFromXYZ(jag_id,res=c(0.002607436,0.002607436),crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',digits=0.3)

###########################################################
###Plot maps of parameters and boxplot of regressino coefficients
###########################################################

###########################################################
###Plot maps of parameters and regression coefficients
###########################################################

# Parameter coefficients

gcols <- grep("gamma.beta", colnames(m[[1]]))
gamma.coef <- data.frame(variable = rep(sort(rep(colnames(m[[1]])[gcols],1000)),3),
                         value = as.vector(sapply(m, FUN = function(x){x[,gcols]})),
                         stringsAsFactors = F)

lcols <- grep("lambda.beta", colnames(m[[1]]))
lambda.coef <- data.frame(variable = rep(sort(rep(colnames(m[[1]])[lcols],1000)),3),
                          value = as.vector(sapply(m, FUN = function(x){x[,lcols]})),
                          stringsAsFactors = F)

Acols <- grep("A.beta", colnames(m[[1]]))
A.coef <- data.frame(variable = rep(sort(rep(colnames(m[[1]])[Acols],1000)),3),
                     value = as.vector(sapply(m, FUN = function(x){x[,Acols]})),
                     stringsAsFactors = F)

betas <- bind_rows(lambda.coef, gamma.coef, A.coef)
betas$covariate <- str_sub(betas$variable, -7, -1)
betas$variable <- str_sub(betas$variable, 1, -9)
betas <- filter(betas, !covariate == "beta[1]")
betas$covariate <- recode(betas$covariate, 
                          'beta[2]' = "slope", 
                          'beta[3]' = "aspect",
                          'beta[4]' = "TPI",
                          'beta[5]' = "precip_Jan",
                          'beta[6]' = "precip_July",
                          'beta[7]' = "tmax_Jan",
                          'beta[8]' = "tmin_July",
                          'beta[9]' = "soiltype")

b <- ggplot(betas) +
  geom_boxplot(aes(x = covariate, y = value)) +
  facet_wrap(.~variable, scales = "free_x") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_hline(aes(yintercept = 0), colour = "gray50") +
  coord_flip()

# Maps

dat <- bind_rows(lambdas, gammas, alphas, As)
dat <- left_join(dat, jag_id, by = c("parnum" = "gammas.parnum"))
dat$lon <- as.numeric(as.character(dat$lon))
dat$lat <- as.numeric(as.character(dat$lat))
dat <- fortify(dat)

###Get coastline for pretty plotting

coast <- st_read("Data/coastline") #, layer = "coastline")
coast <- st_transform(coast, '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
coast <- fortify(coast)

###Plot

# Means

l <- ggplot() +
  geom_sf(data = coast, fill = "skyblue1") +
  geom_tile(data = filter(dat, parstr == "lambda"), aes(lon, lat, fill = log(Mean))) + 
  scale_fill_gradient(low = "#FFF5F0", high = "#67000D", na.value = "transparent") +
  theme_void() +
  labs(fill="ln(lambda)")

g <- ggplot() +
  geom_sf(data = coast, fill = "skyblue1") +
  geom_tile(data = filter(dat, parstr == "gamma"), aes(lon, lat, fill = log(Mean))) + 
  scale_fill_gradient(low = "steelblue4", high = "white", na.value = "transparent") +
  theme_void() +
  labs(fill="ln(gamma)")

a <- ggplot() +
  geom_sf(data = coast, fill = "skyblue1") +
  geom_tile(data = filter(dat, parstr == "alpha"), aes(lon, lat, fill = Mean)) + 
  scale_fill_gradient(low = "#EDFAED", high = "#228B22", na.value = "transparent") +
  theme_void() +
  labs(fill="alpha")

A <- ggplot() +
  geom_sf(data = coast, fill = "skyblue1") +
  geom_tile(data = filter(dat, parstr == "A"), aes(lon, lat, fill = Mean)) + 
  scale_fill_gradient(low = "#FAF3E3", high = "#B8860B", na.value = "transparent") +
  theme_void() +
  labs(fill="Alpha")

parmeans <- ggdraw() +
  draw_plot(l + theme(legend.position=c(.25,.25), legend.key.size=unit(.5, "cm")), x = 0, y = 0, width = .25, height = 1) +
  draw_plot(g + theme(legend.position=c(.25,.25), legend.key.size=unit(.5, "cm")), x = .25, y = 0, width = .25, height = 1) +
  draw_plot(A + theme(legend.position=c(.25,.25), legend.key.size=unit(.5, "cm")), x = .5, y = 0, width = .25, height = 1) +
  draw_plot(a + theme(legend.position=c(.25,.25), legend.key.size=unit(.5, "cm")), x = .75, y = 0, width = .25, height = 1)

ggsave(filename = "figures/parametermap_means.png", plot = parmeans, device = NULL, path = NULL, scale = 1, width = 26, height = 14, units = "cm", dpi = 300, limitsize = TRUE)

pars <- ggdraw() +
  draw_plot(A + theme(legend.position=c(.25,.25), legend.key.size=unit(.5, "cm")), x = 0, y = .4, width = .33, height = .6) +
  draw_plot(g + theme(legend.position=c(.25,.25), legend.key.size=unit(.5, "cm")), x = .33, y = .4, width = .33, height = .6) +
  draw_plot(l + theme(legend.position=c(.25,.25), legend.key.size=unit(.5, "cm")), x = .66, y = .4, width = .33, height = .6) +
  draw_plot(b + theme(legend.position=c(.25,.25), legend.key.size=unit(.5, "cm")), x = 0, y = 0, width = .9, height = .4) +
  draw_plot_label(label = c("(a)", "(b)"), size = 15, x = c(0.025, 0.025), y = c(.965, .4), fontface = "bold")

ggsave(filename = "parametermap.png", plot = pars, device = NULL, path = NULL, scale = 1, width = 16, height = 20, units = "cm", dpi = 300, limitsize = TRUE)


# Standard Deviations

l <- ggplot() +
  geom_sf(data = coast, fill = "skyblue1") +
  geom_tile(data = filter(dat, parstr == "lambda"), aes(lon, lat, fill = log(SD))) + 
  scale_fill_gradient(low = "#FFF5F0", high = "#67000D", na.value = "transparent") +
  theme_void() +
  labs(fill="ln(lambda)")

g <- ggplot() +
  geom_sf(data = coast, fill = "skyblue1") +
  geom_tile(data = filter(dat, parstr == "gamma"), aes(lon, lat, fill = log(SD))) + 
  scale_fill_gradient(low = "steelblue4", high = "white", na.value = "transparent") +
  theme_void() +
  labs(fill="ln(gamma)")

a <- ggplot() +
  geom_sf(data = coast, fill = "skyblue1") +
  geom_tile(data = filter(dat, parstr == "alpha"), aes(lon, lat, fill = SD)) + 
  scale_fill_gradient(low = "#EDFAED", high = "#228B22", na.value = "transparent") +
  theme_void() +
  labs(fill="alpha")

A <- ggplot() +
  geom_sf(data = coast, fill = "skyblue1") +
  geom_tile(data = filter(dat, parstr == "A"), aes(lon, lat, fill = SD)) + 
  scale_fill_gradient(low = "#FAF3E3", high = "#B8860B", na.value = "transparent") +
  theme_void() +
  labs(fill="Alpha")

parsd <- ggdraw() +
  draw_plot(l + theme(legend.position=c(.25,.25), legend.key.size=unit(.5, "cm")), x = 0, y = 0, width = .25, height = 1) +
  draw_plot(g + theme(legend.position=c(.25,.25), legend.key.size=unit(.5, "cm")), x = .25, y = 0, width = .25, height = 1) +
  draw_plot(A + theme(legend.position=c(.25,.25), legend.key.size=unit(.5, "cm")), x = .5, y = 0, width = .25, height = 1) +
  draw_plot(a + theme(legend.position=c(.25,.25), legend.key.size=unit(.5, "cm")), x = .75, y = 0, width = .25, height = 1)

ggsave(filename = "figures/parametermap_SD.png", plot = parsd, device = NULL, path = NULL, scale = 1, width = 26, height = 14, units = "cm", dpi = 300, limitsize = TRUE)


###########################################################
###Plot observed NDVI time-series with model predictions for points of interest
###########################################################

# Get points of interest and extract "id" by intersecting with "jag_id_ras" raster
pts <- data.frame(
  Site = c("Miller's Point: Alien clearing", "Cape of Good Hope: Drought", 
  "Silvermine: Development", "Silvermine: Drought",
  "Silvermine: Aliens", "Karbonkelberg: Fire"),
  Latitude = c(-34.225279, -34.316068, -34.113911,
               -34.115129, -34.03837, -34.039321),
  Longitude = c(18.463551, 18.428294, 18.39334,
                18.397745, 18.37375, 18.334125)
  )
  
coordinates(pts) <- ~ Longitude + Latitude
ids <- extract(jag_id_ras, pts)

###NOTE: if you dont want to do this for all pixels it will take a long time!!!###

# Filter on the the jag_id column of tdat and env
tdat <- tdat %>% filter(jag_id %in% ids)
env <- env %>% filter(jag_id %in% ids)

#set number of samples
nsamp <- 1000

#sample global parameters
phi_par <- rnorm(nsamp, phi$Mean,phi$SD)
sigma_par <- rnorm(nsamp, sig$Mean,sig$SD)

#make an empty df to which we can write results
new_dat = tdat[FALSE,]

# This loop models ndvi using the estimated parameters and covariate data for the points of interest
##NB model was fitted on data up to 2014-05-31##

##time
tic()

for (j in 1:nrow(env)){
  #pixel_id
  pid <- env$jag_id[j]
  
  tdat_temp <- tdat %>% filter(jag_id == pid)
  
  #get pixel pars
  alpha_m <- alphas$Mean[which(alphas$parnum==pid)]
  alpha_sd <- alphas$SD[which(alphas$parnum==pid)]
  
  lambda_m <- lambdas$Mean[which(lambdas$parnum==pid)]
  lambda_sd <- lambdas$SD[which(lambdas$parnum==pid)]
  
  A_m <- As$Mean[which(As$parnum==pid)]
  A_sd <- As$SD[which(As$parnum==pid)]
  
  gamma_m <- gammas$Mean[which(gammas$parnum==pid)]
  gamma_sd <- gammas$SD[which(gammas$parnum==pid)]
  
  #sample
  alpha <- rnorm(nsamp, alpha_m,alpha_sd)
  lambda <- rnorm(nsamp, lambda_m,lambda_sd)
  gamma <- rnorm(nsamp, gamma_m,gamma_sd)
  A <- rnorm(nsamp, A_m,A_sd)
  
  for (i in 1:nrow(tdat_temp)){
    #sample loop
    mu <- numeric(nsamp)
    ndvi <- numeric(nsamp)
    
    for (k in 1:nsamp){
      mu[k] <- alpha[k]+gamma[k]-gamma[k]*exp(-(tdat_temp$DA[i]/lambda[k]))+
        sin((phi_par[k]+((tdat_temp$firemonth[i]-1)*3.141593/6))+6.283185*tdat_temp$DA[i])*A[k]
      
      ndvi[k] <- rnorm(1,mu[k], sigma_par[k])
    }
    #summarize samples
    upper <- quantile(ndvi,probs=0.975)
    uq <- quantile(ndvi,probs=0.75) 
    mean <- quantile(ndvi,probs=0.5)
    lq <- quantile(ndvi,probs=0.25)
    lower <- quantile(ndvi,probs=0.025)
    #write to df
    tdat_temp$mean[i] <- mean
    tdat_temp$upper[i] <- upper
    tdat_temp$lower[i] <- lower
    tdat_temp$uq[i] <- uq
    tdat_temp$lq[i] <- lq
  }
  
#output the final results data frame
new_dat <- bind_rows(new_dat,tdat_temp)
}

toc()
  
### Plot

#Fix names
nms <- data.frame(jag_id = ids, Name = pts$Site[1:6])
new_dat <- merge(new_dat, nms)

#cape point
Pcp <- ggplot(data=new_dat, aes(x=Date,y=NDVI)) +
  #geom_point() +
  geom_line(color="blue") +
  geom_ribbon(aes(ymin=lower,ymax=upper,alpha=0.1))+
  geom_ribbon(aes(ymin=lq,ymax=uq,alpha=0.1))+
  scale_x_date(date_breaks = "1 year",
               labels=date_format("%Y"),
               limits = as.Date(c('2013-01-01','2017-06-20'))) +
  scale_y_continuous(limits=c(0,1)) +
  facet_wrap(~Name) +
  xlab("Date") +
  ylab("NDVI") +
  theme_bw() +
  theme(legend.position="none") +
  geom_vline(xintercept = as.Date("2014-05-31")) +
  annotate("text", label = "Fit", x = as.Date("2013-03-15"), y = 0.1) +
  annotate("text", label = "Forecast", x = as.Date("2015-06-01"), y = 0.1)

ggsave(filename = "figures/postfire_curves_points_of_interest.png", plot = Pcp, device = NULL, path = NULL, scale = 1, width = 18, height = 12, units = "cm", dpi = 300, limitsize = TRUE)

    
    
    
    
    



