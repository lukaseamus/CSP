#####################################################################
##### Project: Climate-driven shifts in kelp forest composition #####
#####          reduce carbon sequestration potential            #####
##### Script purpose: Estimation of carbon sequestration        #####
##### Author: Luka Seamus Wright                                #####
#####################################################################

#### 1.    Input ####
#### 1.1   Temperature ####
E <- read.csv("~/Desktop/Plymouth University/Dissertation/Publication/Data/Environmental.csv")
fieldtemp <- mean(E$temp) # 13.49867°C was the mean temperature during decomposition

#### 1.2   Light ####
# gross primary production
ogpp <- function(x) {
          if(x < 37){
            -3.476632e-05 * x + 1.264358e-03
          } else {
            0
          }
        } # L. ochroleuca carbon assimilation (g C g-1 h-1)
ogpp <- Vectorize(ogpp) # allow function to take vectors as input

hgpp <- function(x) {
          if(x < 108){
            -1.649451e-05 * x + 1.766794e-03
            } else {
              0
            }
        } # L. hyperborea carbon assimilation (g C g-1 h-1)
hgpp <- Vectorize(hgpp) # allow function to take vectors as input

dgpp <- function(x) {6.797230e-06 * x + 1.266816e-03} # L. digitata carbon assimilation (g C g-1 h-1)
# this function does not need to be constrained because it's slope is positive

# photosynthetically active radiation
L4PAR <- function(x) {exp(-0.1520991 * x + 5.0476945)} # annual relationship of PAR with depth at L4
L4spPAR <- function(x) {exp(-0.1636840 * x + 5.7091730)} # spring relationship of PAR with depth at L4
L4suPAR <- function(x) {exp(-0.1496887 * x + 5.8214088)} # summer relationship of PAR with depth at L4
L4auPAR <- function(x) {exp(-0.1447491 * x + 4.6215877)} # autumn relationship of PAR with depth at L4
L4wiPAR <- function(x) {exp(-0.1484400 * x + 3.5215956)} # winter relationship of PAR with depth at L4

# photosynthesis-irradiance
PI <- function(x) {1.8 * tanh(0.006 * x / 1.8)} # relationship between photosynthesis and PAR (doi 10.3354/ab00515)

# seabed slope
WH <- cbind(-4.144978, 50.363629) # West Hoe coordinates
L4 <-  c(-4.189722222, 50.22297222) # L4 coordinates (doi 10.1002/ecm.1366)
require(geosphere)
dist <- distm(WH, L4, fun = distGeo)
sl <- 48/dist # seabed slope = -0.003006109

# daylight hours
dh <- 12 # annual average daylight hours
spdh <- 13.672233333333333
sudh <- 15.577766666666667
audh <- 10.711133333333333
widh <- 8.877766666666667

# Effect of PAR on GPP
# While 0.842857142857143 m d-1 has been found to be the mean detrital velocity 
# for L. hyp. in Norway (doi 10.1007/s00442-018-4121-7), we find that in order to 
# reach the sink site where L. dig. and L. hyp. have been found (doi 10.1002/ecm.1366),
# detritus must travel there in about 100 days at a constant decomposition rate.
# The nearest kelp forest is Rame Head 
RH <- c(-4.221083, 50.31164) # Rame Head coordinates (doi 10.1002/ecm.1366)
v <- sqrt(distm(RH, L4, fun = distGeo)^2 + 48^2)/50
# detritus needs to travel at a velocity of 202.2655 m d-1 to reach the sink site on time
# assuming it needs to get there in 50 days 

# daytime irradiance during experiment
i <- mean(E[E$daytime == "day",]$lux)*0.019 # 0.54 μmol photons m-2 s-1
# 0.019 is the conversion factor from lux to μmol

PAR <- function(time, slope, period) {
       m <- (1 - PI(i)/PI(50.4)) * # 50.4 μmol photons m-2 s-1 was the lab PAR
            PI(L4PAR(time*v/sqrt(slope^2 + 1)*slope))/PI(50.4) 
       sp <- (1 - PI(i)/PI(50.4)) *
             PI(L4spPAR(time*v/sqrt(slope^2 + 1)*slope))/PI(50.4)
       su <- (1 - PI(i)/PI(50.4)) *
             PI(L4suPAR(time*v/sqrt(slope^2 + 1)*slope))/PI(50.4)
       au <- (1 - PI(i)/PI(50.4)) *
             PI(L4auPAR(time*v/sqrt(slope^2 + 1)*slope))/PI(50.4)
       wi <- (1 - PI(i)/PI(50.4)) *
             PI(L4wiPAR(time*v/sqrt(slope^2 + 1)*slope))/PI(50.4)
      
       if (period == "year") {
         m
       }
       else if (period == "spring") {
         sp
       }
       else if (period == "summer") {
         su
       }
       else if (period == "autumn") {
         au
       }
       else if (period == "winter") {
         wi
       }
} # change in photosynthesis due to PAR (m-1 depth) as a function of seabed slope and number of days

#### 1.3   Decomposition ####
# increase in decomposition with temperature = 0.00051 proportion °C-1 (doi 10.21203/rs.3.rs-38503/v1)
od <- function(time) {1 - time * 0.0153889922394093} # L. ochroleuca decomposition (proportion left d-1)
hd <- function(time) {1 - time * 0.00603222583027008} # L. hyperborea decomposition (proportion left d-1)
dd <- function(time) {1 - time * 0.00926997162947394} # L. digitata decomposition (proportion left d-1)

# constrained decomposition (proportion left >= 0)
od0 <- function(time) {
  deco <- 1 - time * 0.0153889922394093
  if(time < 65){
    deco
  } else {
    0
  }
} # L. ochroleuca decomposition (proportion left d-1)
od0 <- Vectorize(od0)

hd0 <- function(time) {
  deco <- 1 - time * 0.00603222583027008
  if(time < 166){
    deco
  } else {
    0
  }
} # L. hyperborea decomposition (proportion left d-1)
hd0 <- Vectorize(hd0)

dd0 <- function(time) {
  deco <- 1 - time * 0.00926997162947394
  if(time < 108){
    deco
  } else {
    0
  }
} # L. digitata decomposition (proportion left d-1)
dd0 <- Vectorize(dd0)

#### 1.4   Export ####
c <- read.csv("~/Desktop/Plymouth University/Dissertation/Publication/Data/Constants.csv")
# Biomass export (g plant-1 yr-1), carbon export (g C plant-1 yr-1), 2017 plant density (plants m-2),
# product 95% confidence interval (g C m-2 yr-1) constants for all species and periods

#### 2.    Estimation functions ####
#### 2.1   Carbon export function ####
# variables: density = plants m-2
#            species = "ochroleuca", "hyperborea" or "digitata"
#            period = annual or seasons 

CE <- function (density, species, period) {
  d <- density * c$carbon[1]
  h <- density * c$carbon[2]
  o <- density * c$carbon[3]
  
  dsp <- density * c$carbon[4]
  hsp <- density * c$carbon[5]
  osp <- density * c$carbon[6]
  
  dsu <- density * c$carbon[7]
  hsu <- density * c$carbon[8]
  osu <- density * c$carbon[9]
  
  dau <- density * c$carbon[10]
  hau <- density * c$carbon[11]
  oau <- density * c$carbon[12]

  dwi <- density * c$carbon[13]
  hwi <- density * c$carbon[14]
  owi <- density * c$carbon[15]
  
  if (species == "ochroleuca" && period == "year") {
    o
  }
  else if (species == "hyperborea" && period == "year") {
    h
  }
  else if (species == "digitata" && period == "year") {
    d
  }
  else if (species == "ochroleuca" && period == "spring") {
    osp
  }
  else if (species == "hyperborea" && period == "spring") {
    hsp
  }
  else if (species == "digitata" && period == "spring") {
    dsp
  }
  else if (species == "ochroleuca" && period == "summer") {
    osu
  }
  else if (species == "hyperborea" && period == "summer") {
    hsu
  }
  else if (species == "digitata" && period == "summer") {
    dsu
  }
  else if (species == "ochroleuca" && period == "autumn") {
    oau
  }
  else if (species == "hyperborea" && period == "autumn") {
    hau
  }
  else if (species == "digitata" && period == "autumn") {
    dau
  }
  else if (species == "ochroleuca" && period == "winter") {
    owi
  }
  else if (species == "hyperborea" && period == "winter") {
    hwi
  }
  else if (species == "digitata" && period == "winter") {
    dwi
  }
  else {
    stop("Not a valid species.")
  }
}

#### 2.2   Carbon sequestraion potential function ####
# variables: density = plants m-2
#            species = "ochroleuca", "hyperborea" or "digitata"
#            period = annual or seasons 
#            time = vector of detrital age
#            temperature = temperature correction of decomposition 

CSP <- function (density, species, period, time) {
  
  d <- density * c$carbon[1] * dd(time)
  h <- density * c$carbon[2] * hd(time)
  o <- density * c$carbon[3] * od(time)

  dsp <- density * c$carbon[4] * dd(time)
  hsp <- density * c$carbon[5] * hd(time)
  osp <- density * c$carbon[6] * od(time)
  
  dsu <- density * c$carbon[7] * dd(time)
  hsu <- density * c$carbon[8] * hd(time)
  osu <- density * c$carbon[9] * od(time)
  
  dau <- density * c$carbon[10] * dd(time)
  hau <- density * c$carbon[11] * hd(time)
  oau <- density * c$carbon[12] * od(time)
  
  dwi <- density * c$carbon[13] * dd(time)
  hwi <- density * c$carbon[14] * hd(time)
  owi <- density * c$carbon[15] * od(time)

  if (species == "ochroleuca" && period == "year") {
    o
  }
  else if (species == "hyperborea" && period == "year") {
    h
  }
  else if (species == "digitata" && period == "year") {
    d
  }
  else if (species == "ochroleuca" && period == "spring") {
    osp
  }
  else if (species == "hyperborea" && period == "spring") {
    hsp
  }
  else if (species == "digitata" && period == "spring") {
    dsp
  }
  else if (species == "ochroleuca" && period == "summer") {
    osu
  }
  else if (species == "hyperborea" && period == "summer") {
    hsu
  }
  else if (species == "digitata" && period == "summer") {
    dsu
  }
  else if (species == "ochroleuca" && period == "autumn") {
    oau
  }
  else if (species == "hyperborea" && period == "autumn") {
    hau
  }
  else if (species == "digitata" && period == "autumn") {
    dau
  }
  else if (species == "ochroleuca" && period == "winter") {
    owi
  }
  else if (species == "hyperborea" && period == "winter") {
    hwi
  }
  else if (species == "digitata" && period == "winter") {
    dwi
  }
  else {
    stop("Not a valid species.")
  }
}

#### 2.3   Carbon assimilation function ####
# variables: species = "ochroleuca", "hyperborea" or "digitata"
#            period = annual or seasons 
#            time = vector of detrital age
#            slope = seabed slope

CA <- function (species, period, time, slope) {
  
  time <- 1:time
  
  d <- sum(c$density[1] * c$biomass[1] * dd0(time) * # exported biomass left after n days
           dh * dgpp(time) * PAR(time, slope, "year")) # GPP at n days
         
  h <- sum(c$density[2] * c$biomass[2] * hd0(time) *
           dh * hgpp(time) * PAR(time, slope, "year"))
  o <- sum(c$density[3] * c$biomass[3] * od0(time) *
           dh * ogpp(time) * PAR(time, slope, "year"))

  dsp <- sum(c$density[4] * c$biomass[4] * dd0(time) *
             spdh * dgpp(time) * PAR(time, slope, "spring"))
  hsp <- sum(c$density[5] * c$biomass[5] * hd0(time) *
             spdh * hgpp(time) * PAR(time, slope, "spring"))
  osp <- sum(c$density[6] * c$biomass[6] * od0(time) *
             spdh * ogpp(time) * PAR(time, slope, "spring"))

  dsu <- sum(c$density[7] * c$biomass[7] * dd0(time) *
             sudh * dgpp(time) * PAR(time, slope, "summer"))
  hsu <- sum(c$density[8] * c$biomass[8] * hd0(time) *
             sudh * hgpp(time) * PAR(time, slope, "summer"))
  osu <- sum(c$density[9] * c$biomass[9] * od0(time) *
             sudh * ogpp(time) * PAR(time, slope, "summer"))

  dau <- sum(c$density[10] * c$biomass[10] * dd0(time) *
             audh * dgpp(time) * PAR(time, slope, "autumn"))
  hau <- sum(c$density[11] * c$biomass[11] * hd0(time) *
             audh * hgpp(time) * PAR(time, slope, "autumn"))
  oau <- sum(c$density[12] * c$biomass[12] * od0(time) *
             audh * ogpp(time) * PAR(time, slope, "autumn"))

  dwi <- sum(c$density[13] * c$biomass[13] * dd0(time) *
             widh * dgpp(time) * PAR(time, slope, "winter"))
  hwi <- sum(c$density[14] * c$biomass[14] * hd0(time) *
             widh * hgpp(time) * PAR(time, slope, "winter"))
  owi <- sum(c$density[15] * c$biomass[15] * od0(time) *
             widh * ogpp(time) * PAR(time, slope, "winter"))
  
  if (species == "ochroleuca" && period == "year") {
    o
  }
  else if (species == "ochroleuca" && period == "spring") {
    osp
  }
  else if (species == "ochroleuca" && period == "summer") {
    osu
  }
  else if (species == "ochroleuca" && period == "autumn") {
    oau
  }
  else if (species == "ochroleuca" && period == "winter") {
    owi
  }
  else if (species == "hyperborea" && period == "year") {
    h
  }
  else if (species == "hyperborea" && period == "spring") {
    hsp
  }
  else if (species == "hyperborea" && period == "summer") {
    hsu
  }
  else if (species == "hyperborea" && period == "autumn") {
    hau
  }
  else if (species == "hyperborea" && period == "winter") {
    hwi
  }
  else if (species == "digitata" && period == "year") {
    d
  }
  else if (species == "digitata" && period == "spring") {
    dsp
  }
  else if (species == "digitata" && period == "summer") {
    dsu
  }
  else if (species == "digitata" && period == "autumn") {
    dau
  }
  else if (species == "digitata" && period == "winter") {
    dwi
  }
  else {
    stop("Not a valid input.")
  }
}

CA <- Vectorize(CA)

#### 3.    Estimations ####
#### 3.1   Present CSP ####
#### 3.1.1 Species ####
dCSP <- CSP(c$density[1], "digitata", "year", 0:250)
hCSP <- CSP(c$density[2], "hyperborea", "year", 0:250)
oCSP <- CSP(c$density[3], "ochroleuca", "year", 0:250)

# build dataframe
pseq <- data.frame(time = 0:250,
                   CSP = c(dCSP, hCSP, oCSP),
                   lo = c(dCSP-c$cCI[1], hCSP-c$cCI[2], oCSP-c$cCI[3]),
                   hi = c(dCSP+c$cCI[1], hCSP+c$cCI[2], oCSP+c$cCI[3]),
                   species = c(rep("d", 251), rep("h", 251), rep("o", 251)))

# calculate equation for each species
# L. digitata
c$density[1] * c$carbon[1] # y-intercept = 89.84651
c$density[1] * c$carbon[1] - c$density[1] * c$carbon[1] * dd(1) # slope = 0.8328746
# y = -0.83x + 89.85

# L. hyperborea
c$density[2] * c$carbon[2] # y-intercept = 211.1804
c$density[2] * c$carbon[2] - c$density[2] * c$carbon[2] * hd(1) # slope = 1.273888
# y = -1.27x + 211.18

# L. ochroleuca
c$density[3] * c$carbon[3] # y-intercept = 127.2262
c$density[3] * c$carbon[3] - c$density[3] * c$carbon[3] * od(1) # slope = 1.957883
# y = -1.96x + 127.23

# calculate x intercept ± 95% CI for each species
xint <- function(x1, x2, y1, y2){
  -(x2*y1-x1*y2)/(x2-x1) / (y2-y1)/(x2-x1)
}

# L. digitata
with(pseq, xint(time[1], time[2], CSP[1], CSP[2])) # mean = 107.8752 days
with(pseq, xint(time[1], time[2], hi[1], hi[2])) # upper = 132.3463 days
with(pseq, xint(time[1], time[2], lo[1], lo[2])) # lower = 83.40409 days

# L. hyperborea
with(pseq, xint(time[252], time[253], CSP[252], CSP[253])) # mean = 165.7763 days
with(pseq, xint(time[252], time[253], hi[252], hi[253])) # upper = 207.8915 days
with(pseq, xint(time[252], time[253], lo[252], lo[253])) # lower = 123.661 days

# L. ochroleuca
with(pseq, xint(time[503], time[504], CSP[503], CSP[504])) # mean = 64.98151 days
with(pseq, xint(time[503], time[504], hi[503], hi[504])) # upper = 86.59629 days
with(pseq, xint(time[503], time[504], lo[503], lo[504])) # lower = 43.36674 days

#### 3.1.2 Seasons ####
spdCSP <- CSP(c$density[4], "digitata", "spring", 0:400)
sphCSP <- CSP(c$density[5], "hyperborea", "spring", 0:400)
spoCSP <- CSP(c$density[6], "ochroleuca", "spring", 0:400)
sudCSP <- CSP(c$density[7], "digitata", "summer", 0:400)
suhCSP <- CSP(c$density[8], "hyperborea", "summer", 0:400)
suoCSP <- CSP(c$density[9], "ochroleuca", "summer", 0:400)
audCSP <- CSP(c$density[10], "digitata", "autumn", 0:400)
auhCSP <- CSP(c$density[11], "hyperborea", "autumn", 0:400)
auoCSP <- CSP(c$density[12], "ochroleuca", "autumn", 0:400)
widCSP <- CSP(c$density[13], "digitata", "winter", 0:400)
wihCSP <- CSP(c$density[14], "hyperborea", "winter", 0:400)
wioCSP <- CSP(c$density[15], "ochroleuca", "winter", 0:400)

# build dataframe
pseqs <- data.frame(time = 0:400,
                    CSP = c(spdCSP, sphCSP, spoCSP, sudCSP,
                            suhCSP, suoCSP, audCSP, auhCSP,
                            auoCSP, widCSP, wihCSP, wioCSP),
                    lo = c(spdCSP-c$cCI[4], sphCSP-c$cCI[5], spoCSP-c$cCI[6],
                           sudCSP-c$cCI[7], suhCSP-c$cCI[8], suoCSP-c$cCI[9],
                           audCSP-c$cCI[10], auhCSP-c$cCI[11],auoCSP-c$cCI[12],
                           widCSP-c$cCI[13], wihCSP-c$cCI[14], wioCSP-c$cCI[15]),
                    hi = c(spdCSP+c$cCI[4], sphCSP+c$cCI[5], spoCSP+c$cCI[6],
                           sudCSP+c$cCI[7], suhCSP+c$cCI[8], suoCSP+c$cCI[9],
                           audCSP+c$cCI[10], auhCSP+c$cCI[11],auoCSP+c$cCI[12],
                           widCSP+c$cCI[13], wihCSP+c$cCI[14], wioCSP+c$cCI[15]),
                    species = rep(c(rep("d", 401), rep("h", 401),
                                    rep("o", 401)),4),
                    season = c(rep("Spring", 1203), rep("Summer", 1203),
                               rep("Autumn", 1203), rep("Winter", 1203)))


# Spring
with(pseqs, xint(time[1], time[2], CSP[1], CSP[2])) # 107.8752 
with(pseqs, xint(time[402], time[403], CSP[402], CSP[403])) # 165.7763 
with(pseqs, xint(time[803], time[804], CSP[803], CSP[804])) # 64.98151 

# Summer
with(pseqs, xint(time[1204], time[1205], CSP[1204], CSP[1205])) # 107.8752 
with(pseqs, xint(time[1605], time[1606], CSP[1605], CSP[1606])) # 165.7763 
with(pseqs, xint(time[2006], time[2007], CSP[2006], CSP[2007])) # 64.98151 

# Detrital longevity is the same across seasons

#### 3.1   Present CA ####
#### 3.1.1 Species ####
dCA <- CA("digitata", "year", 0:250, sl)
hCA <- CA("hyperborea", "year", 0:250, sl)
oCA <- CA("ochroleuca", "year", 0:250, sl)
dCA[1] <- 0
hCA[1] <- 0
oCA[1] <- 0

# build dataframe
pass <- data.frame(time = 0:250,
                   CA = c(dCA, hCA, oCA),
                   lo = c(dCA-c$cCI[1], hCA-c$cCI[2], oCA-c$cCI[3]),
                   hi = c(dCA+c$cCI[1], hCA+c$cCI[2], oCA+c$cCI[3]),
                   species = c(rep("d", 251), rep("h", 251), rep("o", 251)))


#### 3.1.2 Seasons ####
spdCA <- CA("digitata", "spring", 0:400, sl)
sphCA <- CA("hyperborea", "spring", 0:400, sl)
spoCA <- CA("ochroleuca", "spring", 0:400, sl)
sudCA <- CA("digitata", "summer", 0:400, sl)
suhCA <- CA("hyperborea", "summer", 0:400, sl)
suoCA <- CA("ochroleuca", "summer", 0:400, sl)
audCA <- CA("digitata", "autumn", 0:400, sl)
auhCA <- CA("hyperborea", "autumn", 0:400, sl)
auoCA <- CA("ochroleuca", "autumn", 0:400, sl)
widCA <- CA("digitata", "winter", 0:400, sl)
wihCA <- CA("hyperborea", "winter", 0:400, sl)
wioCA <- CA("ochroleuca", "winter", 0:400, sl)
spdCA[1] <- 0
sphCA[1] <- 0
spoCA[1] <- 0
sudCA[1] <- 0
suhCA[1] <- 0
suoCA[1] <- 0
audCA[1] <- 0
auhCA[1] <- 0
auoCA[1] <- 0
widCA[1] <- 0
wihCA[1] <- 0
wioCA[1] <- 0

# build dataframe
passs <- data.frame(time = 0:400,
                    CA = c(spdCA, sphCA, spoCA, sudCA,
                            suhCA, suoCA, audCA, auhCA,
                            auoCA, widCA, wihCA, wioCA),
                    lo = c(spdCA-c$cCI[4], sphCA-c$cCI[5], spoCA-c$cCI[6],
                           sudCA-c$cCI[7], suhCA-c$cCI[8], suoCA-c$cCI[9],
                           audCA-c$cCI[10], auhCA-c$cCI[11],auoCA-c$cCI[12],
                           widCA-c$cCI[13], wihCA-c$cCI[14], wioCA-c$cCI[15]),
                    hi = c(spdCA+c$cCI[4], sphCA+c$cCI[5], spoCA+c$cCI[6],
                           sudCA+c$cCI[7],suhCA+c$cCI[8], suoCA+c$cCI[9],
                           audCA+c$cCI[10], auhCA+c$cCI[11],auoCA+c$cCI[12],
                           widCA+c$cCI[13], wihCA+c$cCI[14], wioCA+c$cCI[15]),
                    species = rep(c(rep("d", 401), rep("h", 401),
                                    rep("o", 401)),4),
                    season = c(rep("Spring", 1203), rep("Summer", 1203),
                               rep("Autumn", 1203), rep("Winter", 1203)))


#### 3.2   Past, present and future CSP ####
# contribution of each species to the kelp forest is altered according
# to temperature predicted by future climate scenarios

#### 3.2.1 Load future temperature data ####
require(raster)
# Data source: https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html (doi 10.1029/2002jd002670)
# HadISST sea surface temperature data between 1870 and 2019
had.sst <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/HadISST_sst.nc")
# Data source: https://www.bio-oracle.org/downloads-to-email.php (doi 10.1111/geb.12693)
# Bio-ORACLE minimum sea surface temperature data for 2050 (RCP8.5)
min2050.85 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2050AOGCM.RCP85.Surface.Temperature.Lt.min.asc.BOv2_1.asc")
# Bio-ORACLE mean sea surface temperature data for 2050 (RCP8.5)
mean2050.85 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2050AOGCM.RCP85.Surface.Temperature.Mean.asc.BOv2_1.asc") 
# Bio-ORACLE maximum sea surface temperature data for 2050 (RCP8.5)
max2050.85 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2050AOGCM.RCP85.Surface.Temperature.Lt.max.asc.BOv2_1.asc")
# Bio-ORACLE minimum sea surface temperature data for 2100 (RCP8.5)
min2100.85 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2100AOGCM.RCP85.Surface.Temperature.Lt.min.asc.BOv2_1.asc") 
# Bio-ORACLE mean sea surface temperature data for 2100 (RCP8.5)
mean2100.85 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2100AOGCM.RCP85.Surface.Temperature.Mean.asc.BOv2_1.asc")
# Bio-ORACLE maximum sea surface temperature data for 2100 (RCP8.5)
max2100.85 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2100AOGCM.RCP85.Surface.Temperature.Lt.max.asc.BOv2_1.asc") 
# Bio-ORACLE minimum sea surface temperature data for 2050 (RCP6.0)
min2050.60 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2050AOGCM.RCP60.Surface.Temperature.Lt.min.asc.BOv2_1.asc")
# Bio-ORACLE mean sea surface temperature data for 2050 (RCP6.0)
mean2050.60 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2050AOGCM.RCP60.Surface.Temperature.Mean.asc.BOv2_1.asc") 
# Bio-ORACLE maximum sea surface temperature data for 2050 (RCP6.0)
max2050.60 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2050AOGCM.RCP60.Surface.Temperature.Lt.max.asc.BOv2_1.asc")
# Bio-ORACLE minimum sea surface temperature data for 2100 (RCP6.0)
min2100.60 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2100AOGCM.RCP60.Surface.Temperature.Lt.min.asc.BOv2_1.asc") 
# Bio-ORACLE mean sea surface temperature data for 2100 (RCP6.0)
mean2100.60 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2100AOGCM.RCP60.Surface.Temperature.Mean.asc.BOv2_1.asc")
# Bio-ORACLE maximum sea surface temperature data for 2100 (RCP6.0)
max2100.60 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2100AOGCM.RCP60.Surface.Temperature.Lt.max.asc.BOv2_1.asc") 
# Bio-ORACLE minimum sea surface temperature data for 2050 (RCP2.6)
min2050.26 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2050AOGCM.RCP26.Surface.Temperature.Lt.min.asc.BOv2_1.asc")
# Bio-ORACLE mean sea surface temperature data for 2050 (RCP2.6)
mean2050.26 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2050AOGCM.RCP26.Surface.Temperature.Mean.asc.BOv2_1.asc") 
# Bio-ORACLE maximum sea surface temperature data for 2050 (RCP2.6)
max2050.26 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2050AOGCM.RCP26.Surface.Temperature.Lt.max.asc.BOv2_1.asc")
# Bio-ORACLE minimum sea surface temperature data for 2100 (RCP2.6)
min2100.26 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2100AOGCM.RCP26.Surface.Temperature.Lt.min.asc.BOv2_1.asc") 
# Bio-ORACLE mean sea surface temperature data for 2100 (RCP2.6)
mean2100.26 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2100AOGCM.RCP26.Surface.Temperature.Mean.asc.BOv2_1.asc")
# Bio-ORACLE maximum sea surface temperature data for 2100 (RCP2.6)
max2100.26 <- brick("~/Desktop/Plymouth University/Dissertation/Publication/Data/2100AOGCM.RCP26.Surface.Temperature.Lt.max.asc.BOv2_1.asc") 

# HadISST February data for West Hoe
Feb <- had.sst[[seq(2, 1790, by = 12)]] # extract February for each year
min1 <- extract(Feb, WH, method = "bilinear") # extract data
min1 <- t(min1) # transpose data

# HadISST August data for West Hoe
Aug <- had.sst[[seq(8, 1796, by = 12)]] # extract August for each year
max1 <- extract(Aug, WH, method = "bilinear") # extract data
max1 <- t(max1) # transpose data

# HadISST mean data for West Hoe
M <- had.sst[[1:1800]] # extract each month and year
mean1 <- extract(M, WH, method = "bilinear") # extract data
mean1 <- data.frame(mean = t(mean1),
                    month = 1:12)
mean1$names <- rownames(mean1)
mean1$year <- substr(mean1$names, start = 2, stop = 5)
mean1 <- aggregate(mean ~ year, mean1, mean) # average for each year

## RCP8.5
# Bio-ORACLE 2050 February data for West Hoe
min2.85 <- extract(min2050.85, WH, method = "bilinear")

# Bio-ORACLE 2050 August data for West Hoe
max2.85 <- extract(max2050.85, WH, method = "bilinear")

# Bio-ORACLE 2050 mean data for West Hoe
mean2.85 <- extract(mean2050.85, WH, method = "bilinear")

# Bio-ORACLE 2100 February data for West Hoe
min3.85 <- extract(min2100.85, WH, method = "bilinear")

# Bio-ORACLE 2100 August data for West Hoe
max3.85 <- extract(max2100.85, WH, method = "bilinear")

# Bio-ORACLE 2100 mean data for West Hoe
mean3.85 <- extract(mean2100.85, WH, method = "bilinear")

## RCP2.6
# Bio-ORACLE 2050 February data for West Hoe
min2.26 <- extract(min2050.26, WH, method = "bilinear")

# Bio-ORACLE 2050 August data for West Hoe
max2.26 <- extract(max2050.26, WH, method = "bilinear")

# Bio-ORACLE 2050 mean data for West Hoe
mean2.26 <- extract(mean2050.26, WH, method = "bilinear")

# Bio-ORACLE 2100 February data for West Hoe
min3.26 <- extract(min2100.26, WH, method = "bilinear")

# Bio-ORACLE 2100 August data for West Hoe
max3.26 <- extract(max2100.26, WH, method = "bilinear")

# Bio-ORACLE 2100 mean data for West Hoe
mean3.26 <- extract(mean2100.26, WH, method = "bilinear")

## RCP6.0
# Bio-ORACLE 2050 February data for West Hoe
min2.60 <- extract(min2050.60, WH, method = "bilinear")

# Bio-ORACLE 2050 August data for West Hoe
max2.60 <- extract(max2050.60, WH, method = "bilinear")

# Bio-ORACLE 2050 mean data for West Hoe
mean2.60 <- extract(mean2050.60, WH, method = "bilinear")

# Bio-ORACLE 2100 February data for West Hoe
min3.60 <- extract(min2100.60, WH, method = "bilinear")

# Bio-ORACLE 2100 August data for West Hoe
max3.60 <- extract(max2100.60, WH, method = "bilinear")

# Bio-ORACLE 2100 mean data for West Hoe
mean3.60 <- extract(mean2100.60, WH, method = "bilinear")

# Combine into one dataframe
sst.85 <- data.frame(year = rep(c(1870:2019, 2050, 2100), 3),
                     temp = c(rbind(min1, min2.85, min3.85), rbind(max1, max2.85, max3.85),
                              rbind(mean1, mean2.85, mean3.85)[, 2]), 
                     group = c(rep("Minimum", 152), rep("Maximum", 152), 
                               rep("Mean", 152)),
                     row.names = c(1:456))

sst.26 <- data.frame(year = rep(c(1870:2019, 2050, 2100), 3),
                     temp = c(rbind(min1, min2.26, min3.26), rbind(max1, max2.26, max3.26),
                              rbind(mean1, mean2.26, mean3.26)[, 2]), 
                     group = c(rep("Minimum", 152), rep("Maximum", 152), 
                               rep("Mean", 152)),
                     row.names = c(1:456))

sst.60 <- data.frame(year = rep(c(1870:2019, 2050, 2100), 3),
                     temp = c(rbind(min1, min2.60, min3.60), rbind(max1, max2.60, max3.60),
                              rbind(mean1, mean2.60, mean3.60)[, 2]), 
                     group = c(rep("Minimum", 152), rep("Maximum", 152), 
                               rep("Mean", 152)),
                     row.names = c(1:456))

# Remove pre-1900 data
sst.85 <- data.frame(sst.85[-c(1:30, 153:182, 305:334),], row.names = NULL)
sst.26 <- data.frame(sst.26[-c(1:30, 153:182, 305:334),], row.names = NULL)
sst.60 <- data.frame(sst.60[-c(1:30, 153:182, 305:334),], row.names = NULL)

#### 3.2.2 Visualise temperature data for West Hoe ####
# Combine dataframes
sst <- data.frame(rbind(sst.26, sst.60, sst.85),
                  scenario = c(rep("RCP2.6", 366), rep("RCP6.0", 366), rep("RCP8.5", 366)))

# Define custom theme
require(ggplot2)
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .5, .2, .2),"cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.background = element_blank(),
                 legend.text.align = 0,
                 legend.title = element_blank(),
                 text = element_text(family = "Helvetica Neue"))

sstp <- ggplot(sst, aes(year, temp)) +
          annotate("segment", x = Inf, y = Inf, xend = -Inf, yend = Inf, lwd = 1) +
          geom_point(aes(colour = group), shape = 16, alpha = 0.4, show.legend = F) +
          geom_smooth(aes(fill = group, colour = group),
                     method = "loess", span = 1, size = 0.5, alpha = 0.5,
                     show.legend = F) +
          scale_fill_manual(values = c("#f5a54a", "#898b8e", "#6ea4be")) +
          scale_colour_manual(values = c("#f5a54a", "#898b8e", "#6ea4be")) +
          annotate("segment", x = Inf, y = 20, xend = -Inf, yend = 20, lwd = 0.5) + 
          # upper tolerance for cold temperate kelps (predicted from experimentals)
          annotate("segment", x = Inf, y = 9.53, xend = -Inf, yend = 9.53, lwd = 0.5) + 
          # lower tolerance for warm temperate kelp (predicted to be 10°C from experiments 
          # but arrival in Plymouth Sound suggests 9.53°C)
          ylab("Sea surface temperature (°C)") +
          xlab("Year") +
          scale_x_continuous(expand = c(0,0), breaks = seq(1900, 2100, by = 100)) +
          scale_y_continuous(expand = c(0,0), breaks = seq(6, 22, by = 4)) +
          coord_cartesian(ylim = c(6, 22), xlim = c(1900, 2100)) +
          facet_grid(~scenario) +
          theme(legend.position = c(0.78, 0.9),
                strip.background = element_blank(),
                strip.text = element_text(size = 15, hjust = 0),
                panel.spacing = unit(1.5, "cm")) +
          mytheme

sstp # dimensions 4 x 8.5 in

#### 3.2.3 Extract fitted values as seen in plot ####
## RCP8.5
min.sst <- sst.85[1:122,]
min.mod <- loess(data = min.sst, temp ~ year, span = 1)
min.fit <- data.frame(year = 1900:2100,
                      min.temp = predict(min.mod, newdata = 1900:2100),
                      group = rep("Minimum", 201))

max.sst <- sst.85[123:244,]
max.mod <- loess(data = max.sst, temp ~ year, span = 1)
max.fit <- data.frame(year = 1900:2100,
                      max.temp = predict(max.mod, newdata = 1900:2100),
                      group = rep("Maximum", 201))

mean.sst <- sst.85[245:366,]
mean.mod <- loess(data = mean.sst, temp ~ year, span = 1)
mean.fit <- data.frame(year = 1900:2100,
                       mean.temp = predict(mean.mod, newdata = 1900:2100),
                       group = rep("Mean", 201))

sst.fit.85 <- cbind(min.fit, max.fit, mean.fit)[,c(1,2,5,8)]

## RCP2.6
min.sst <- sst.26[1:122,]
min.mod <- loess(data = min.sst, temp ~ year, span = 1)
min.fit <- data.frame(year = 1900:2100,
                      min.temp = predict(min.mod, newdata = 1900:2100),
                      group = rep("Minimum", 201))

max.sst <- sst.26[123:244,]
max.mod <- loess(data = max.sst, temp ~ year, span = 1)
max.fit <- data.frame(year = 1900:2100,
                      max.temp = predict(max.mod, newdata = 1900:2100),
                      group = rep("Maximum", 201))

mean.sst <- sst.26[245:366,]
mean.mod <- loess(data = mean.sst, temp ~ year, span = 1)
mean.fit <- data.frame(year = 1900:2100,
                       mean.temp = predict(mean.mod, newdata = 1900:2100),
                       group = rep("Mean", 201))

sst.fit.26 <- cbind(min.fit, max.fit, mean.fit)[,c(1,2,5,8)]

## RCP6.0
min.sst <- sst.60[1:122,]
min.mod <- loess(data = min.sst, temp ~ year, span = 1)
min.fit <- data.frame(year = 1900:2100,
                      min.temp = predict(min.mod, newdata = 1900:2100),
                      group = rep("Minimum", 201))

max.sst <- sst.60[123:244,]
max.mod <- loess(data = max.sst, temp ~ year, span = 1)
max.fit <- data.frame(year = 1900:2100,
                      max.temp = predict(max.mod, newdata = 1900:2100),
                      group = rep("Maximum", 201))

mean.sst <- sst.60[245:366,]
mean.mod <- loess(data = mean.sst, temp ~ year, span = 1)
mean.fit <- data.frame(year = 1900:2100,
                       mean.temp = predict(mean.mod, newdata = 1900:2100),
                       group = rep("Mean", 201))

sst.fit.60 <- cbind(min.fit, max.fit, mean.fit)[,c(1,2,5,8)]


#### 3.2.4 Calculate contribution of each species ####

# thermal limits of the different Laminaria spp.:
# L. hyp. and L. dig. have 20°C as their upper limit for sporophyte growth
# and 18°C and 17°C as their upper limit for gametophyte fertility.
# Because we are talking about the trailing edge of these species, the
# ceasing of sporophyte growth is taken as the upper limit: 20.1°C.
# L. och. has 5°C as the lower limit for sporophyte growth
# and 10°C as its lower limit for gametophyte fertility.
# Because we are talking about the leading edge of this species, gametophyte
# fertility is taken as the lower limit: 10°C. Since we have accurate data
# on the arrival of L. och. in Plymouth Sound in 1946 (Parke, 1948), this value 
# can be slightly adjusted to the minimum temperature in the region in 1945, 
# one year before it was found
sst.fit.85$min.temp[46] # 9.523089°C
sst.fit.60$min.temp[46] # 9.528703°C
sst.fit.26$min.temp[46] # 9.525496°C

## RCP8.5
# calculate proportion of 2016 density remaining for each year
up.lim <- 20.1 - sst.fit.85$max.temp
lo.lim <- sst.fit.85$min.temp - sst.fit.85$min.temp[46]
# standardising all values by the 2016 value gives us a proportion we 
# can multiply the 2016 density by
up.prop <- up.lim/up.lim[117] 
lo.prop <- lo.lim/lo.lim[117]

# calculate density contribution across years based on calculated proportion
sst.fit.85$dp <- c$density[1]*up.prop
sst.fit.85$hp <- c$density[2]*up.prop
sst.fit.85$op <- c$density[3]*lo.prop

# the kelp forest at West Hoe is split into a 16-m mixed band of L. hyperborea and L. ochroleuca
# and a 4-m band of L. digitata
# L. ochroleuca cannot occupy the L. digitata zone so it can only take up
# the space occupied by L. hyperborea once this species disappears (doi: 10.1111/1365-2435.12977)
# Using 2016 density data, the maximum density for L. ochroleuca is estimates as 12.09863 m-2

sst.fit.85$op[181:201] <- 12.09863 # L. ochroleuca densities maxes out at 12.09863 m-1

## RCP2.6
# calculate proportion of 2016 density remaining for each year
up.lim <- 20.1 - sst.fit.26$max.temp
lo.lim <- sst.fit.26$min.temp - sst.fit.26$min.temp[46]
# standardising all values by the 2016 value gives us a proportion we 
# can multiply the 2016 density by
up.prop <- up.lim/up.lim[117] 
lo.prop <- lo.lim/lo.lim[117]

# calculate density contribution across years based on calculated proportion
sst.fit.26$dp <- c$density[1]*up.prop
sst.fit.26$hp <- c$density[2]*up.prop
sst.fit.26$op <- c$density[3]*lo.prop

## RCP6.0
# calculate proportion of 2016 density remaining for each year
up.lim <- 20.1 - sst.fit.60$max.temp
lo.lim <- sst.fit.60$min.temp - sst.fit.60$min.temp[46]
# standardising all values by the 2016 value gives us a proportion we 
# can multiply the 2016 density by
up.prop <- up.lim/up.lim[117] 
lo.prop <- lo.lim/lo.lim[117]

# calculate density contribution across years based on calculated proportion
sst.fit.60$dp <- c$density[1]*up.prop
sst.fit.60$hp <- c$density[2]*up.prop
sst.fit.60$op <- c$density[3]*lo.prop

#### 3.2.5 Species ####
## RCP8.5
# 50 days is chosen because detritus of all species is still present at that stage
sst.fit.85$dCSP <- with(sst.fit.85, CSP(dp, "digitata", "year", 50))
sst.fit.85$hCSP <- with(sst.fit.85, CSP(hp, "hyperborea", "year", 50))
sst.fit.85$oCSP <- with(sst.fit.85, CSP(op, "ochroleuca", "year", 50))

# build dataframe
fseq.85 <- with(sst.fit.85,
                data.frame(time = 1900:2100,
                           CSP = c(dCSP, hCSP, oCSP),
                           lo = c(dCSP-c$cCI[1], hCSP-c$cCI[2], oCSP-c$cCI[3]),
                           hi = c(dCSP+c$cCI[1], hCSP+c$cCI[2], oCSP+c$cCI[3]),
                           species = c(rep("d", 201), rep("h", 201), rep("o", 201))))

## RCP2.6
# 50 days is chosen because detritus of all species is still present at that stage
sst.fit.26$dCSP <- with(sst.fit.26, CSP(dp, "digitata", "year", 50))
sst.fit.26$hCSP <- with(sst.fit.26, CSP(hp, "hyperborea", "year", 50))
sst.fit.26$oCSP <- with(sst.fit.26, CSP(op, "ochroleuca", "year", 50))

# build dataframe
fseq.26 <- with(sst.fit.26,
                data.frame(time = 1900:2100,
                           CSP = c(dCSP, hCSP, oCSP),
                           lo = c(dCSP-c$cCI[1], hCSP-c$cCI[2], oCSP-c$cCI[3]),
                           hi = c(dCSP+c$cCI[1], hCSP+c$cCI[2], oCSP+c$cCI[3]),
                           species = c(rep("d", 201), rep("h", 201), rep("o", 201))))

## RCP6.0
# 50 days is chosen because detritus of all species is still present at that stage
sst.fit.60$dCSP <- with(sst.fit.60, CSP(dp, "digitata", "year", 50))
sst.fit.60$hCSP <- with(sst.fit.60, CSP(hp, "hyperborea", "year", 50))
sst.fit.60$oCSP <- with(sst.fit.60, CSP(op, "ochroleuca", "year", 50))

# build dataframe
fseq.60 <- with(sst.fit.60,
                data.frame(time = 1900:2100,
                           CSP = c(dCSP, hCSP, oCSP),
                           lo = c(dCSP-c$cCI[1], hCSP-c$cCI[2], oCSP-c$cCI[3]),
                           hi = c(dCSP+c$cCI[1], hCSP+c$cCI[2], oCSP+c$cCI[3]),
                           species = c(rep("d", 201), rep("h", 201), rep("o", 201))))



#### 3.2.6 Total ####
## RCP8.5
# set L. digitata and L. hyperborea CSP to 0 from 2090 onward
sst.fit.85$dCSP[191:201] <- 0 
sst.fit.85$hCSP[191:201] <- 0 
# set L. ochroleuca CSP to 0 before 1945
sst.fit.85$oCSP[1:45] <- 0

sst.fit.85$CSP <- with(sst.fit.85, dCSP + hCSP + oCSP)

# build dataframe
fseqt.85 <- with(sst.fit.85, 
                 data.frame(time = 1900:2100,
                            CSP = CSP))

fseqt.85$lo <-  fseqt.85$CSP - sqrt(c$cCI[1]^2 + c$cCI[2]^2 + c$cCI[3]^2)
fseqt.85$hi <-  fseqt.85$CSP + sqrt(c$cCI[1]^2 + c$cCI[2]^2 + c$cCI[3]^2)

# estimated linear rate
with(fseqt.85, (CSP[1] - CSP[201])/200) # decrease of 0.9251831 (g C m-2 yr-1) per year
with(fseqt.85, (CSP[1] - CSP[201])/200/CSP[1]*100) # decrease of 0.330093% yr-1
with(fseqt.85, (CSP[1] - CSP[117])/116/CSP[1]*100) # decrease of 0.1699618% yr-1
with(fseqt.85, (CSP[117] - CSP[201])/84/CSP[117]*100) # decrease of 0.6865921% yr-1

## RCP2.6
# set L. ochroleuca CSP to 0 before 1945 and after 2094
sst.fit.26$oCSP[c(1:45, 196:201)] <- 0

sst.fit.26$CSP <- with(sst.fit.26, dCSP + hCSP + oCSP)

# build dataframe
fseqt.26 <- with(sst.fit.26, 
                 data.frame(time = 1900:2100,
                            CSP = CSP))

fseqt.26$lo <-  fseqt.26$CSP - sqrt(c$cCI[1]^2 + c$cCI[2]^2 + c$cCI[3]^2)
fseqt.26$hi <-  fseqt.26$CSP + sqrt(c$cCI[1]^2 + c$cCI[2]^2 + c$cCI[3]^2)

# estimated linear rate
with(fseqt.26, (CSP[1] - CSP[201])/200) # decrease of 0.6063763 (g C m-2 yr-1) per year
with(fseqt.26, (CSP[1] - CSP[201])/200/CSP[1]*100) # decrease of 0.2230919% yr-1

## RCP6.0
# set L. ochroleuca CSP to 0 before 1945
sst.fit.60$oCSP[1:45] <- 0

sst.fit.60$CSP <- with(sst.fit.60, dCSP + hCSP + oCSP)

# build dataframe
fseqt.60 <- with(sst.fit.60, 
                 data.frame(time = 1900:2100,
                            CSP = CSP))

fseqt.60$lo <-  fseqt.60$CSP - sqrt(c$cCI[1]^2 + c$cCI[2]^2 + c$cCI[3]^2)
fseqt.60$hi <-  fseqt.60$CSP + sqrt(c$cCI[1]^2 + c$cCI[2]^2 + c$cCI[3]^2)

# estimated linear rate
with(fseqt.60, (CSP[1] - CSP[201])/200) # decrease of 0.6866185 (g C m-2 yr-1) per year
with(fseqt.60, (CSP[1] - CSP[201])/200/CSP[1]*100) # decrease of 0.2514098% yr-1

# comparison between scenarios
(with(fseqt.85, CSP[117] - CSP[201]) - with(fseqt.60, CSP[117] - CSP[201]))/with(fseqt.85, CSP[117] - CSP[201])
(with(fseqt.85, CSP[117] - CSP[201]) - with(fseqt.26, CSP[117] - CSP[201]))/with(fseqt.85, CSP[117] - CSP[201])

#### 3.3   Past, present and future CE ####
#### 3.3.1 Species ####
## RCP8.5
sst.fit.85$dCE <- CE(sst.fit.85$dp, "digitata", "year")
sst.fit.85$hCE <- CE(sst.fit.85$hp, "hyperborea", "year")
sst.fit.85$oCE <- CE(sst.fit.85$op, "ochroleuca", "year")

## RCP6.0
sst.fit.60$dCE <- CE(sst.fit.60$dp, "digitata", "year")
sst.fit.60$hCE <- CE(sst.fit.60$hp, "hyperborea", "year")
sst.fit.60$oCE <- CE(sst.fit.60$op, "ochroleuca", "year")

## RCP2.6
sst.fit.26$dCE <- CE(sst.fit.26$dp, "digitata", "year")
sst.fit.26$hCE <- CE(sst.fit.26$hp, "hyperborea", "year")
sst.fit.26$oCE <- CE(sst.fit.26$op, "ochroleuca", "year")

# build dataframe
fexp <- data.frame(time = 1900:2100,
                   CE = c(with(sst.fit.26, c(dCE, hCE, oCE)),
                          with(sst.fit.60, c(dCE, hCE, oCE)),
                          with(sst.fit.85, c(dCE, hCE, oCE))),
                   lo = c(with(sst.fit.26, c(dCE-c$cCI[1], hCE-c$cCI[2], oCE-c$cCI[3])),
                          with(sst.fit.60, c(dCE-c$cCI[1], hCE-c$cCI[2], oCE-c$cCI[3])),
                          with(sst.fit.85, c(dCE-c$cCI[1], hCE-c$cCI[2], oCE-c$cCI[3]))),
                   hi = c(with(sst.fit.26, c(dCE+c$cCI[1], hCE+c$cCI[2], oCE+c$cCI[3])),
                          with(sst.fit.60, c(dCE+c$cCI[1], hCE+c$cCI[2], oCE+c$cCI[3])),
                          with(sst.fit.85, c(dCE+c$cCI[1], hCE+c$cCI[2], oCE+c$cCI[3]))),
                   scenario = c(rep("RCP2.6", 603), rep("RCP6.0", 603), rep("RCP8.5", 603)),
                   species = rep(c(rep("d", 201), rep("h", 201), rep("o", 201)), 3))

#### 3.3.2 Total ####
## RCP8.5
sst.fit.85$dCE[191:201] <- 0
sst.fit.85$hCE[191:201] <- 0
sst.fit.85$oCE[1:45] <- 0
sst.fit.85$CE <- with(sst.fit.85, dCE + hCE + oCE)

## RCP6.0
sst.fit.60$oCE[1:45] <- 0
sst.fit.60$CE <- with(sst.fit.60, dCE + hCE + oCE)

## RCP2.6
sst.fit.26$oCE[c(1:45, 196:201)] <- 0
sst.fit.26$CE <- with(sst.fit.26, dCE + hCE + oCE)

# build dataframe
fexpt <- data.frame(time = 1900:2100,
                    CE = c(sst.fit.26$CE, sst.fit.60$CE, sst.fit.85$CE),
                    lo = c(sst.fit.26$CE - sqrt(c$cCI[1]^2 + c$cCI[2]^2 + c$cCI[3]^2),
                           sst.fit.60$CE - sqrt(c$cCI[1]^2 + c$cCI[2]^2 + c$cCI[3]^2),
                           sst.fit.85$CE - sqrt(c$cCI[1]^2 + c$cCI[2]^2 + c$cCI[3]^2)),
                    hi = c(sst.fit.26$CE + sqrt(c$cCI[1]^2 + c$cCI[2]^2 + c$cCI[3]^2),
                           sst.fit.60$CE + sqrt(c$cCI[1]^2 + c$cCI[2]^2 + c$cCI[3]^2),
                           sst.fit.85$CE + sqrt(c$cCI[1]^2 + c$cCI[2]^2 + c$cCI[3]^2)),
                    scenario = c(rep("RCP2.6", 201), rep("RCP6.0", 201), rep("RCP8.5", 201)))

# estimated linear rate
with(fexpt, (CE[403] - CE[603])/200) # decrease of 0.09020656 (g C m-2 yr-1) per year
with(fexpt, (CE[403] - CE[603])/200/CE[403]*100) # decrease of 0.0209222% yr-1

#### 4.    Visualisation ####
#### 4.1   Present CSP ####
#### 4.1.1 Species ####
pcs <- ggplot(pseq, aes(time, CSP)) +
          geom_line(aes(colour = species), size = 0.5) +
          geom_ribbon(aes(ymin = lo, ymax = hi, fill = species), alpha = 0.5) +
          scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                              labels = c(expression(italic("L. digitata")*"        y = –0.83x + 89.85"),
                                         expression(italic("L. hyperborea")*"  y = –1.27x + 211.18"),
                                         expression(italic("L. ochroleuca")*"  y = –1.96x + 127.23")),
                              guide = guide_legend()) +
          scale_fill_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")*"        y = –0.83x + 89.85"),
                                       expression(italic("L. hyperborea")*"  y = –1.27x + 211.18"),
                                       expression(italic("L. ochroleuca")*"  y = –1.96x + 127.23")),
                            guide = guide_legend()) +
          xlab("Detrital age (d)") +
          ylab(expression("Carbon sequestration (g C m"^-2*" yr"^-1*")")) +
          scale_x_continuous(expand = c(0,0)) +
          scale_y_continuous(expand = c(0,0)) +
          coord_cartesian(ylim = c(0, 400), xlim = c(0, 250)) +
          theme(legend.position = c(0.5, 0.905)) +
          mytheme

pcs

#### 4.1.2 Seasons ####
# Relevel seasons to avoid alphabetical order
pseqs$season <- factor(pseqs$season, levels = c("Spring", "Summer", "Autumn", "Winter"))

pcss <- ggplot(pseqs, aes(time, CSP)) +
          annotate("segment", x = Inf, y = Inf, xend = -Inf, yend = Inf, lwd = 1) +
          geom_line(aes(colour = species), size = 0.5) +
          geom_ribbon(aes(ymin = lo, ymax = hi, fill = species), alpha = 0.5) +
          scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                              labels = c(expression(italic("L. digitata")),
                                         expression(italic("L. hyperborea")),
                                         expression(italic("L. ochroleuca"))),
                              guide = guide_legend()) +
          scale_fill_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")),
                                       expression(italic("L. hyperborea")),
                                       expression(italic("L. ochroleuca"))),
                            guide = guide_legend()) +
          xlab("Detrital age (d)") +
          ylab(expression("Carbon sequestration (g C m"^-2*")")) +
          scale_x_continuous(expand = c(0,0), breaks = seq(0, 300, by = 100)) +
          scale_y_continuous(expand = c(0,0)) +
          coord_cartesian(ylim = c(0, 200), xlim = c(0, 300)) +
          facet_grid(~season) +
          theme(legend.position = c(0.91, 0.86),
                strip.background = element_blank(),
                strip.text = element_text(size = 15, hjust = 0),
                panel.spacing = unit(.8, "cm")) +
          mytheme

pcss # save as 4 x 8.5 in

#### 4.2   Present cumulative CA ####
#### 4.2.1 Species ####
pca <- ggplot(pass, aes(time, CA)) +
          geom_line(aes(colour = species), size = 0.5) +
          geom_ribbon(aes(ymin = lo, ymax = hi, fill = species), alpha = 0.5) +
          scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                              labels = c(expression(italic("L. digitata")),
                                         expression(italic("L. hyperborea")),
                                         expression(italic("L. ochroleuca"))),
                              guide = guide_legend()) +
          scale_fill_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")),
                                       expression(italic("L. hyperborea")),
                                       expression(italic("L. ochroleuca"))),
                            guide = guide_legend()) +
          xlab("Detrital age (d)") +
          ylab(expression("Carbon assimilation (g C m"^-2*" yr"^-1*")")) +
          scale_x_continuous(expand = c(0,0)) +
          scale_y_continuous(expand = c(0,0)) +
          coord_cartesian(ylim = c(0, 500), xlim = c(0, 100)) +
          theme(legend.position = c(0.1, 0.905)) +
          mytheme

pca # save as 4 x 8.5 in

#### 4.2.2 Seasons ####
# Relevel seasons to avoid alphabetical order
passs$season <- factor(passs$season, levels = c("Spring", "Summer", "Autumn", "Winter"))

pcas <- ggplot(passs, aes(time, CA)) +
          annotate("segment", x = Inf, y = Inf, xend = -Inf, yend = Inf, lwd = 1) +
          geom_line(aes(colour = species), size = 0.5) +
          geom_ribbon(aes(ymin = lo, ymax = hi, fill = species), alpha = 0.5) +
          scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                              labels = c(expression(italic("L. digitata")),
                                         expression(italic("L. hyperborea")),
                                         expression(italic("L. ochroleuca"))),
                              guide = guide_legend()) +
          scale_fill_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")),
                                       expression(italic("L. hyperborea")),
                                       expression(italic("L. ochroleuca"))),
                            guide = guide_legend()) +
          xlab("Detrital age (d)") +
          ylab(expression("Carbon assimilation (g C m"^-2*")")) +
          scale_x_continuous(expand = c(0,0), breaks = seq(0, 40, by = 20)) +
          scale_y_continuous(expand = c(0,0)) +
          coord_cartesian(ylim = c(0, 600), xlim = c(0, 40)) +
          facet_grid(~season) +
          theme(legend.position = c(0.9, 0.87),
                strip.background = element_blank(),
                strip.text = element_text(size = 15, hjust = 0),
                panel.spacing = unit(.8, "cm")) +
          mytheme

pcas

#### 4.3   Past, present and future CSP ####
#### 4.3.1 Species ####
fcs <- ggplot(fseq.85, aes(time, CSP)) +
          geom_line(aes(colour = species), size = 0.5) +
          geom_ribbon(aes(ymin = lo, ymax = hi, fill = species), alpha = 0.5) +
          scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                              labels = c(expression(italic("L. digitata")),
                                         expression(italic("L. hyperborea")),
                                         expression(italic("L. ochroleuca"))),
                              guide = guide_legend()) +
          scale_fill_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")),
                                       expression(italic("L. hyperborea")),
                                       expression(italic("L. ochroleuca"))),
                            guide = guide_legend()) +
          xlab("Year") +
          ylab(expression("Carbon sequestration (g C m"^-2*" yr"^-1*")")) +
          scale_x_continuous(expand = c(0,0)) +
          scale_y_continuous(expand = c(0,0)) +
          coord_cartesian(ylim = c(0, 250), xlim = c(1900, 2100)) +
          theme(legend.position = c(0.89, 0.9)) +
          mytheme

fcs

#### 4.3.2 Total ####
## alternative RCP scenarios
fseqt.alt <- data.frame(rbind(fseqt.26, fseqt.60),
                        scenario = c(rep("RCP2.6", 201),
                                     rep("RCP6.0", 201)),
                        row.names = NULL)

fcst2 <- ggplot(fseqt.alt, aes(time, CSP)) +
          annotate("segment", x = Inf, y = Inf, xend = -Inf, yend = Inf, lwd = 1) +
          geom_line(colour = "#f5a54a", size = 0.5) +
          geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#f5a54a", alpha = 0.5) +
          geom_vline(xintercept = 1946) +
          geom_vline(data = fseqt.alt[fseqt.alt$scenario == "RCP2.6",], aes(xintercept = 2095)) +
          ylab(expression("Carbon sequestration (g C m"^-2*" yr"^-1*")")) +
          xlab("Year") +
          scale_x_continuous(expand = c(0,0), breaks = seq(2000, 2100, by = 50)) +
          scale_y_continuous(expand = c(0,0)) +
          coord_cartesian(ylim = c(0, 400), xlim = c(2000, 2100)) +
          facet_grid(~scenario) +
          theme(legend.position = c(0.78, 0.9),
                strip.background = element_blank(),
                strip.text = element_text(size = 15, hjust = 0),
                panel.spacing = unit(1.5, "cm")) +
          mytheme

fcst2 # dimensions 4 x 6 in

## RCP8.5
fcst <- ggplot(fseqt.85, aes(time, CSP)) +
          geom_line(colour = "#f5a54a", size = 0.5) +
          geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#f5a54a", alpha = 0.5) +
          geom_vline(xintercept = 1946) +
          geom_vline(xintercept = 2090) +
          geom_rug(fseqt.alt[c(201, 402),], mapping = aes(CSP), colour = "#f5a54a") +
          geom_text(fseqt.alt[c(201, 402),],
                    mapping = aes(1925, CSP, label = scenario), 
                    size = 4.2, vjust = c(0, 1), family = "Helvetica Neue") +
          ylab(expression("Carbon sequestration (g C m"^-2*" yr"^-1*")")) +
          xlab("Year") +
          scale_x_continuous(expand = c(0,0)) +
          scale_y_continuous(expand = c(0,0)) +
          coord_cartesian(ylim = c(0, 400), xlim = c(1900, 2100)) +
          theme(legend.position = c(0.54, 0.94)) +
          mytheme

fcst


#### 4.4   Past, present and future CE ####
fce <- ggplot(fexp, aes(time, CE)) +
          annotate("segment", x = Inf, y = Inf, xend = -Inf, yend = Inf, lwd = 1) +
          geom_line(aes(colour = species), size = 0.5) +
          geom_line(data = fexpt, aes(time, CE)) +
          geom_ribbon(data = fexpt, aes(ymin = lo, ymax = hi), alpha = 0.5) +
          geom_ribbon(aes(ymin = lo, ymax = hi, fill = species), alpha = 0.5) +
          scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                              labels = c(expression(italic("L. digitata")),
                                         expression(italic("L. hyperborea")),
                                         expression(italic("L. ochroleuca"))),
                              guide = guide_legend()) +
          scale_fill_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")),
                                       expression(italic("L. hyperborea")),
                                       expression(italic("L. ochroleuca"))),
                            guide = guide_legend()) +
          facet_grid(~scenario) +
          xlab("Year") +
          ylab(expression("Carbon export (g C m"^-2*" yr"^-1*")")) +
          scale_x_continuous(expand = c(0,0), breaks = seq(1900, 2100, by = 100)) +
          scale_y_continuous(expand = c(0,0), breaks = seq(0, 700, by = 100)) +
          coord_cartesian(ylim = c(0, 700), xlim = c(1900, 2100)) +
          theme(legend.position = c(0.098, 0.87),
                strip.background = element_blank(),
                strip.text = element_text(size = 15, hjust = 0),
                panel.spacing = unit(1.5, "cm")) +
          mytheme

fce # save as 4 x 8.5 in



#### 4.5   Final combination of plots ####
require(cowplot)
pcas <- pcas + theme(legend.position = "none")
Fig.S3 <- plot_grid(pcss, pcas, labels = "auto", ncol = 1,
                    label_size = 15, label_fontfamily = "Helvetica Neue")
Fig.S3 # dimensions 8 x 8.5 in

fcst <- fcst + theme(axis.title.y = element_blank(),
                     axis.text.y = element_blank(),
                     plot.margin = unit(c(.2, .5, .2, .5),"cm"))
Fig.4 <- plot_grid(pcs, fcst, labels = "auto", label_fontfamily = "Helvetica Neue",
                   rel_widths = c(1.12, 1), hjust = 0.5, label_size = 15)
Fig.4 # dimensions 4 x 8.5 in

#### 5.   Clean up ####
detach(package:raster)
detach(package:cowplot)
detach(package:ggplot2)
detach(package:geosphere)
rm(list = ls())
graphics.off()
cat("\014")

