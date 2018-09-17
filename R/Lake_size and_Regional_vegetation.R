library(raster)
library(sp)
library(disqover)
library(reshape2)
library(ggplot2)
library(dplyr)
library(neotoma)
library(maps)
library(RColorBrewer)
library(fields)
library(zipfR)
library(rioja)
###################################################################################################################################
#load pollen data
mywd <- '~/git_upload/K/' 
setwd(mywd)
data.loc <-'data/'
plot.loc <-'figures/'

load(paste(data.loc,'elicitation_neus_certainty_median_101_sites_only_abies_new_species.RData',sep=''))
y1 <- y[,!colnames(y)%in%'Chestnut']
#reveals wants and age estimate in column1 (or at least no pollen data)
y2 <- cbind(1:nrow(y1),y1)


taxa <- colnames(y1)

###################################################################################################################################
####################################################################################################
#read REVEALS estimates
####################################################################################################
svs = read.csv(paste(data.loc,'svs_LC6K.csv',sep=''), sep=',', header=TRUE, stringsAsFactors=FALSE)

svs[which(is.na(svs$sv)), 'sv'] = 0.01
svs_agg = aggregate(sv ~ taxon, svs, median)
svs_agg$taxon[svs_agg$taxon=='Alder'] <- 'Other hardwood'
svs_agg$taxon[svs_agg$taxon=='Fir'] <- 'Other conifer'
svs_agg$taxon[svs_agg$taxon=='Larch'] <- 'Tamarack'
svs_agg <- svs_agg[svs_agg$taxon%in%taxa,]
svs_agg <- svs_agg[order(svs_agg$taxon),]

ppes <- readRDS(paste(data.loc,'PPEs_agg.RDS',sep=''))
ppes$taxon <- as.character(ppes$taxon)
ppes$taxon[ppes$taxon=='Alder'] <- 'Other hardwood'
ppes$taxon[ppes$taxon=='Fir'] <- 'Other conifer'
ppes$taxon[ppes$taxon=='Larch'] <- 'Tamarack'
ppes <- ppes[ppes$taxon%in%taxa,]
ppes <- ppes[order(ppes$taxon),]

params <- cbind(ppes,svs_agg$sv)
colnames(params)[3:4] <- c('ppe.error','fallspeed') 
#################################################################################################
#test effects of lake radii
rad.lake <- c(5.5,55,550,5500)
recon.lake.size <-
lapply(rad.lake,function(radius){
a <- REVEALSinR(pollen = y2,#"~/Reveals_NEUS/data/reveals_test.csv",
                params = params,#"~/Reveals_NEUS/data/reveals_input_params.csv",
                dwm        = 'GPM neutral',
                tBasin     = "lake",
                dBasin     = 2*radius, # diameter!
                regionCutoff = 100000,#radius,
                n      = 1000)

recon.reveals <- a[,grep('mean',colnames(a))]
})




lake.radius.diff <- recon.lake.size[[1]] -recon.lake.size[[4]]
max(abs(lake.radius.diff))

################################################################################################################################
#test effects of regional vegetation
rad.veg <- c(50000,100000,200000,300000,400000)
recon.tot <-
  lapply(rad.veg,function(radius){
    a <- REVEALSinR(pollen = y2,#"~/Reveals_NEUS/data/reveals_test.csv",
                    params = params,#"~/Reveals_NEUS/data/reveals_input_params.csv",
                    dwm        = 'GPM neutral',
                    tBasin     = "lake",
                    dBasin     = 2*550, # diameter!
                    regionCutoff = radius,
                    n      = 1000)
    
    recon.reveals <- a[,grep('mean',colnames(a))]
  })

recon.diff <- recon.tot[[1]] -recon.tot[[4]]
max(abs(recon.diff))

################################################################################################################################
#functions to estimate K



#--------------------------------------------------------------------------------------------------------------]
#function to estimate species dependent depositional coefficient
Ki <- function(b, R, zmax=400000){
  ul1 <- b*(zmax-R)^(1/8)
  ul2 <- b*(zmax+R)^(1/8)
  ul3 <- b*(2*R)^(1/8)
  gamma_Ki <- Igamma(8, ul1, lower=TRUE) -  
    Igamma(8, ul2, lower=TRUE) +  
    Igamma(8, ul3, lower=TRUE)
  
  return(4*pi*R/b^8*gamma_Ki)
}

K_calc <- function(fall_speed,R,zmax,n,u,c){
  
  #calculate parameter b later on used to estimate deposition coeffcient 
  b <- 4/sqrt(pi) * fall_speed/(n*u*c)  
  
  #species specific depositional coefficient
  K_species <- 
    sapply(b,function(x){
      Ki(x[[1]],R=R,zmax = zmax)
    })
  
  #eq (5) in Sugita (2007a)
}

n=0.25
u=3
c=0.12
gamma <- n/2


###############################################################################################################
#estimate K
area.h <- c(0.01,1,100,10000)
area.m <- area.h * 10000

radius <- sqrt(area.m/pi)

K_lake_size <-
  sapply(radius,function(x){
    K_calc(fall_speed = params$fallspeed, R = x,zmax = 100000,n=n,u=u,c=c)
  })

rownames(K_lake_size) <- params$taxon

K_lake_size_stand_oak <- t(K_lake_size)/K_lake_size['Oak',]
K_lake_size_mean <- t(K_lake_size)/colMeans(K_lake_size)


matplot(log10(radius),K_lake_size_stand_oak,pch = 1:12,xlim=c(0,4),col=rep(1:6,2))
legend('topleft',pch = 1:12,legend=params$taxon,col=rep(1:6,2))

#################################################################################################################
#different code to estimate K
b <- 4/sqrt(pi) * params$fallspeed/(n*u*c)

K_integration <- 
  sapply(b, function(coef){
    sapply(radius,function(radius){
      b <- coef[1]
      integrand <- function(z)   {b * gamma*z^(gamma-1) * exp(-b*z^gamma)}# * z^2}
      integrate(f=integrand,lower=radius,upper = 100000)$value
    })
  })

colnames(K_integration) <- params$taxon

K_integration_stand_oak <- K_integration/K_integration[,'Oak']

matplot(log10(radius),K_integration_stand_oak,pch = 1:12,xlim=c(0,4))
legend('topleft',pch = 1:12,legend=params$taxon,col=rep(1:6,2))
##################################################################################################################
#load STEPPS parameters
catchment.radii <-readRDS(paste(data.loc,'catchment_radii.RDS',sep=''))
stepps.radii.plk <- catchment.radii[[2]]
stepps.radii.plk  <- stepps.radii.plk[stepps.radii.plk$taxon!='Chestnut',]
radius_text <- paste(c(5.5,55,550,5500),'m lake radius')


pdf(paste(plot.loc,'K_vs_catchment_radius.pdf'),height= 12, width = 7)
par(mfrow=c(2,1))
plot(stepps.radii.plk[,6],K_lake_size_stand_oak[1,],pch = 1:12,
     ylim=c(0,max(K_lake_size_stand_oak)),xlab='',ylab='',xlim=c(100,550))
for(i in 1:4){
  points(stepps.radii.plk[,6],K_lake_size_stand_oak[i,],col=i,pch = 1:12)
}
legend('topleft',legend = c(as.character(stepps.radii.plk$taxon),radius_text),
       pch = c(1:12,rep(15,4)),cex=0.75,col=c(rep(1,12),1:4))
mtext(side=1,line=2.2,'STEPPS 70% catchment radius [km]',font=2)
mtext(side=2,line=2.2,'K',font=2)


plot(stepps.radii.plk[,6],K_lake_size_mean[1,],pch = 1:12,
     ylim=c(0,max(K_lake_size_mean)),xlab='',ylab='',xlim=c(100,550))
for(i in 1:4){
  points(stepps.radii.plk[,6],K_lake_size_mean[i,],col=i,pch = 1:12)
}
legend('topleft',legend = c(as.character(stepps.radii.plk$taxon),radius_text),
       pch = c(1:12,rep(15,4)),cex=0.75,col=c(rep(1,12),1:4))
mtext(side=1,line=2.2,'STEPPS 70% catchment radius [km]',font=2)
mtext(side=2,line=2.2,'K',font=2)
dev.off()

#################################################################################################################
#make vegetation reconstruction
plot(K_lake_size[,3],K_integration[3,])# always nice to call indexes...
cor(K_lake_size[,3],K_integration[3,])

weights <- params$ppe*K_lake_size[,3]

weighted_pollen <- y1/weights
veg <- round(weighted_pollen/rowSums(weighted_pollen),3)
#################################################################################################################
# extract K from Martin Theurkauf's code
source('REVEASL_K.R')
files.reveals <- list.files('~/disqover/R')
sapply(files.reveals[grep('.R',files.reveals)],function(x) source(paste('~/disqover/R/',x,sep='')))

rad.lake <- c(5.5,55,550,5500)
K <-
  sapply(rad.lake,function(radius){
    a <- REVEALSinR_K(pollen = y2[1,],#"~/Reveals_NEUS/data/reveals_test.csv",
                    params = params,#"~/Reveals_NEUS/data/reveals_input_params.csv",
                    dwm        = 'GPM neutral',
                    tBasin     = "lake",
                    dBasin     = 2*radius, # diameter!
                    regionCutoff = 100000,#radius,
                    n      = 1000)
    
    a[[2]]
  })

rownames(K) <- ppes$taxon
colnames(K) <- rad.lake

saveRDS(K,paste(data.loc,'K.RDS',sep=''))

