##################################################################
##### Project: Sequestration of carbon from Laminaria kelp   #####
#####          forests may be diminished in a warmer climate #####
##### Script purpose: Visualisation of carbon export         #####
##### Author: Luka Seamus Wright                             #####
##################################################################

#### 1.   Data preparation ####
#### 1.1  Load data ####
E <- read.csv("~/Desktop/Plymouth University/Dissertation/Publication/Data/Export.csv")
E <- E[-c(307:323),] # remove March 2017 where there is no data for L. ochroleuca
Carb <- read.csv("~/Desktop/Plymouth University/Dissertation/Publication/Data/Carbon.csv")
M <- read.csv("~/Desktop/Plymouth University/Dissertation/Publication/Data/Mass.csv")
DW <- read.csv("~/Desktop/Plymouth University/Dissertation/Publication/Data/DW.csv")

#### 1.2  Rename variables ####
b <- E$dw.export # biomass export (g plant-1 d-1)
C <- E$C.dw.export # carbon export (g C plant-1 d-1)
sp <- E$species
time <- E$time

carbon <- Carb$carbon # lamina carbon content (%)
sp2 <- Carb$species

mass <- M$mass # plant wet mass (g)
time2 <- M$season
sp3 <- M$species

dw <- DW$d.w # dry-wet mass ratio
water <- (1-dw)*100 # lamina tissue water content (%)
month <- DW$month
time3 <- DW$season
sp4 <- DW$species

#### 2.   Data analysis ####
#### 2.1  Build carbon export model ####
m1 <- lm(C ~ sp)

#### 2.2  Test model fit ####
# assumption of homogeneity
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m1, pch = sp) # heterogenous
par(mfrow = c(1,1), mar = c(2,2,2,1))
boxplot(resid(m1) ~ sp) # residual variance varies between species

# assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1)) # extremely right-skewed
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# try fitting alternative distribution

#### 2.3  Fit gamma distribution to data ####
require(fitdistrplus)
gamma <- fitdist(C+0.001, "gamma")
norm <- fitdist(C+0.001, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm),
         legendtext = c("Gamma", "Normal"),
         fitlty = 1)
cdfcomp(list(gamma, norm),
        legendtext = c("Gamma", "Normal"),
        fitlty = 1) # judging visually, gamma fits much better
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# judging statistically, the gamma distribution fits much better

#### 2.4  Determine random components ####
m2 <- glm(C+0.001 ~ sp, family = Gamma(link = "log"))

#### 2.5  Test model fit ####
# assumption of homogeneity
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m2, pch = sp) # homogenous
par(mfrow = c(1,1), mar = c(2,2,2,1))
boxplot(resid(m2) ~ sp) # no variance between species

# assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m2))
qqnorm(resid(m2))
qqline(resid(m2)) # normality is improved
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m2 is chosen as the optimal model

#### 2.6  Interpret model ####
require(car)
Anova(m2, type = 2) # Type II hypothesis test
# Response: C + 0.001
#    LR Chisq Df Pr(>Chisq)
# sp   9.9898  2   0.006772 **

summary(m2)
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  -2.9437     0.1312 -22.429  < 2e-16 ***
# sph           0.3817     0.1835   2.081  0.03830 *
# spo           0.5979     0.1895   3.155  0.00176 **

sp <- factor(sp, levels = c("h", "o", "d"))
m2 <- glm(C+0.001 ~ sp, family = Gamma(link = "log"))
summary(m2)
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  -2.5620     0.1282 -19.989   <2e-16 ***
# spo           0.2161     0.1874   1.154   0.2496
# spd          -0.3817     0.1835  -2.081   0.0383 *

sp <- factor(sp, levels = c("d", "h", "o"))

#### 2.7 Build mass model ####
m3 <- lm(mass ~ sp3)

#### 2.8  Test model fit ####
# assumption of homogeneity
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m3, pch = sp3) # homogenous
par(mfrow = c(1,1), mar = c(2,2,2,1))
boxplot(resid(m3) ~ sp3) # homogenous

# assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m3))
qqnorm(resid(m3))
qqline(resid(m3)) # normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m3 is chosen as the optimal model

#### 2.9  Interpret model ####
Anova(m3, type = 2) # Type II hypothesis test
# Response: mass
#             Sum Sq  Df F value    Pr(>F)
# sp3        2093367   2  23.572 2.982e-10 ***
# Residuals 13676209 308

summary(m3)
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)   315.20      20.37  15.473  < 2e-16 ***
# sp3h          196.78      28.74   6.846 4.11e-11 ***
# sp3o          112.39      29.62   3.794 0.000178 ***

sp3 <- factor(sp3, levels = c("h", "o", "d"))
m3 <- lm(mass ~ sp3)
summary(m3)
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)   511.97      20.28  25.249  < 2e-16 ***
# sp3o          -84.39      29.56  -2.855   0.0046 **
# sp3d         -196.78      28.74  -6.846 4.11e-11 ***

sp3 <- factor(sp3, levels = c("d", "h", "o"))

#### 2.10 Build water content model ####
m4 <- lm(water ~ sp4)

#### 2.11  Test model fit ####
# assumption of homogeneity
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m4, pch = sp4)
par(mfrow = c(1,1), mar = c(2,2,2,1))
boxplot(resid(m4) ~ sp4) # heterogenous

# assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m4))
qqnorm(resid(m4))
qqline(resid(m4)) # slightly left-skewed
# nothing can be done -> accepted as fairly normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# model variance by species

#### 2.12 Build water content model ####
require(nlme)
m5 <- gls(water ~ sp4,
          weights = varIdent(form = ~1|sp4))

#### 2.13  Test model fit ####
# assumption of homogeneity
plot(m5)
boxplot(resid(m5, type = "normalized") ~ sp4) # homogenous

# assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m5))
qqnorm(resid(m5))
qqline(resid(m5)) # normality is the same
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m5 is chosen as the optimal model

#### 2.14  Interpret model ####
Anova(m5, type = 2) # Type II hypothesis test
# Response: water
#     Df  Chisq Pr(>Chisq)    
# sp4  2 210.03  < 2.2e-16 ***

summary(m5)
#                Value Std.Error  t-value p-value
# (Intercept) 85.11495 0.1427667 596.1823  0.0000
# sp4h        -0.71172 0.2897195  -2.4566  0.0146
# sp4o         2.20469 0.1807153  12.1998  0.0000

sp4 <- factor(sp4, levels = c("h", "o", "d"))
m5 <- gls(water ~ sp4,
          weights = varIdent(form = ~1|sp4))
summary(m5)
#                Value Std.Error  t-value p-value
# (Intercept) 84.40324 0.2521013 334.7989  0.0000
# sp4o         2.91640 0.2753738  10.5907  0.0000
# sp4d         0.71172 0.2897195   2.4566  0.0146

sp4 <- factor(sp4, levels = c("d", "h", "o"))

#### 3.   Data visualisation ####
#### 3.1  Descriptive statistics ####
require(psych)
b.stat <- describeBy(b, list(sp, time), mat = T, digits = 10)
b.stat$group2 <- as.integer(b.stat$group2) # turn time back into integer

b2.stat <- describeBy(b, sp, mat = T, digits = 10) 
# annual average

C.stat <- describeBy(C, list(sp, time), mat = T, digits = 10)
C.stat$group2 <- as.integer(C.stat$group2) # turn time back into integer

C2.stat <- describeBy(C, sp, mat = T, digits = 10)
# annual average

Carb.stat <- describeBy(carbon, sp2, mat = T, digits = 10)

mass.stat <- describeBy(mass, list(sp3, time2), mat = T, digits = 10)
mass.stat2 <- describeBy(mass, sp3, mat = T, digits = 10)

dw.stat <- describeBy(dw, list(sp4, month), mat = T, digits = 10)
dw.stat2 <- describeBy(dw, list(sp4, time3), mat = T, digits = 10)

water.stat <- describeBy(water, sp4, mat = T, digits = 10)

#### 3.2  Customised theme ####
require(ggplot2)
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15),
                 axis.title.x = element_blank(),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title = element_blank(),
                 text = element_text(family = "Helvetica Neue"))

#### 3.3  Plots ####
require(ggalt)

bp <- ggplot(data = b.stat) +
        geom_xspline(aes(group2, mean, colour = group1),
                     size = 0.5, spline_shape = -0.3) +
        geom_pointrange(aes(group2, mean, ymin = mean - se,
                            ymax = mean + se, colour = group1), size = 0.5) +
        geom_rug(data = b2.stat, aes(0.8, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm"), outside = TRUE) +
        scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")),
                                       expression(italic("L. hyperborea")),
                                       expression(italic("L. ochroleuca"))),
                            guide = guide_legend()) +
        ylab(expression("Biomass export (g plant"^-1*" d"^-1*")")) +
        scale_x_continuous(breaks = 1:12, 
                           labels = c("Mar","Apr","May","Jun","Jul","Aug",
                                      "Sep","Oct","Nov","Dec","Jan","Feb")) +
        coord_cartesian(ylim = c(0.09, 1.91), clip = "off") +
        theme(legend.position = c(0.89, 0.9)) +
        mytheme

bp # dimsenisons: 4 x 7 in

Cp <- ggplot() +
        geom_xspline(data = C.stat, aes(group2, mean, colour = group1),
                     size = 0.5, spline_shape = -0.3) +
        geom_point(aes(time, C, colour = sp),
                   alpha = 0.4, shape = 16, size = 3, position = position_dodge(width = 0.5)) +
        geom_pointrange(data = C.stat, aes(group2, mean, ymin = mean - se*qnorm(0.975),
                            ymax = mean + se*qnorm(0.975), colour = group1), size = 0.5,
                        position = position_dodge(width = 0.5)) +
        geom_rug(data = C2.stat, aes(0.8, mean, colour = group1), sides = "l",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")),
                                       expression(italic("L. hyperborea")),
                                       expression(italic("L. ochroleuca"))),
                            guide = guide_legend()) +
        ylab(expression("Carbon export (g C plant"^-1*" d"^-1*")")) +
        scale_x_continuous(breaks = 1:12, 
                           labels = c("Mar","Apr","May","Jun","Jul","Aug",
                                      "Sep","Oct","Nov","Dec","Jan","Feb")) +
        scale_y_continuous(expand = c(0,0)) +
        coord_cartesian(ylim = c(0, 0.8)) +
        theme(legend.position = c(0.9, 0.9)) +
        mytheme

Cp # dimsenisons: 4 x 8.5 in

#### 4.    Annual means and variances ####
#### 4.1   Species ####
#### 4.1.1 Biomass export ####
ba <- b * c(rep(31, 28), rep(30, 28), rep(31, 24), rep(30, 27), 
            rep(31, 22), rep(31, 24), rep(30, 22), rep(31, 25),
            rep(30, 24), rep(31, 27), rep(31, 28), rep(28, 27)) 
ba # biomass export (g plant-1 mo-1)
ba.stat <- describeBy(ba, list(sp, time), mat = T, digits = 10)
ba.sum <- aggregate(mean ~ group1, ba.stat, sum) # sum means across months
ba.sum$var <- aggregate(se^2 ~ group1, ba.stat, sum)[,2] # sum variance of the mean across months
ba.sum$se <- sqrt(ba.sum$var) # calculate standard error of the sum
ba.sum$ci <- ba.sum$se * qnorm(0.975) # calculate 95% confidence interval
ba.sum

#### 4.1.2 Carbon export ####
Ca <- C * c(rep(31, 28), rep(30, 28), rep(31, 24), rep(30, 27), 
            rep(31, 22), rep(31, 24), rep(30, 22), rep(31, 25),
            rep(30, 24), rep(31, 27), rep(31, 28), rep(28, 27)) 
Ca # biomass export (g C plant-1 mo-1)
Ca.stat <- describeBy(Ca, list(sp, time), mat = T, digits = 10)
Ca.sum <- aggregate(mean ~ group1, Ca.stat, sum) # sum means across months
aggregate(n ~ group1, Ca.stat, sum)
Ca.sum$var <- aggregate(se^2 ~ group1, Ca.stat, sum)[,2] # sum variance of the mean across months
Ca.sum$se <- sqrt(Ca.sum$var) # calculate standard error of the sum
Ca.sum$ci <- Ca.sum$se * qnorm(0.975) # calculate 95% confidence interval
Ca.sum

#### 4.2   Seasons ####
#### 4.2.1 Biomass export ####
# Spring
bsp.sum <- aggregate(mean ~ group1, ba.stat[1:9,], sum) # sum means across spring months
bsp.sum$var <- aggregate(se^2 ~ group1, ba.stat[1:9,], sum)[,2] # sum variance of the mean across spring months
bsp.sum$se <- sqrt(bsp.sum$var) # calculate standard error of the sum
bsp.sum$ci <- bsp.sum$se * qnorm(0.975) # calculate 95% confidence interval

# Summer
bsu.sum <- aggregate(mean ~ group1, ba.stat[10:18,], sum) # sum means across summer months
bsu.sum$var <- aggregate(se^2 ~ group1, ba.stat[10:18,], sum)[,2] # sum variance of the mean across summer months
bsu.sum$se <- sqrt(bsu.sum$var) # calculate standard error of the sum
bsu.sum$ci <- bsu.sum$se * qnorm(0.975) # calculate 95% confidence interval

# Autumn
bau.sum <- aggregate(mean ~ group1, ba.stat[19:27,], sum) # sum means across autumn months
bau.sum$var <- aggregate(se^2 ~ group1, ba.stat[19:27,], sum)[,2] # sum variance of the mean across autumn months
bau.sum$se <- sqrt(bau.sum$var) # calculate standard error of the sum
bau.sum$ci <- bau.sum$se * qnorm(0.975) # calculate 95% confidence interval

# Winter
bwi.sum <- aggregate(mean ~ group1, ba.stat[28:36,], sum) # sum means across winter months
bwi.sum$var <- aggregate(se^2 ~ group1, ba.stat[28:36,], sum)[,2] # sum variance of the mean across winter months
bwi.sum$se <- sqrt(bwi.sum$var) # calculate standard error of the sum
bwi.sum$ci <- bwi.sum$se * qnorm(0.975) # calculate 95% confidence interval

# Combined
bs.sum <- data.frame(season = c(rep("Spring",3), rep("Summer",3), rep("Autumn",3), rep("Winter",3)),
                     rbind(bsp.sum, bsu.sum, bau.sum, bwi.sum))


#### 4.2.2 Carbon export ####
# Spring
Csp.sum <- aggregate(mean ~ group1, Ca.stat[1:9,], sum) # sum means across spring months
aggregate(n ~ group1, Ca.stat[1:9,], sum)
Csp.sum$var <- aggregate(se^2 ~ group1, Ca.stat[1:9,], sum)[,2] # sum variance of the mean across spring months
Csp.sum$se <- sqrt(Csp.sum$var) # calculate standard error of the sum
Csp.sum$ci <- Csp.sum$se * qnorm(0.975) # calculate 95% confidence interval

# Summer
Csu.sum <- aggregate(mean ~ group1, Ca.stat[10:18,], sum) # sum means across summer months
aggregate(n ~ group1, Ca.stat[10:18,], sum)
Csu.sum$var <- aggregate(se^2 ~ group1, Ca.stat[10:18,], sum)[,2] # sum variance of the mean across summer months
Csu.sum$se <- sqrt(Csu.sum$var) # calculate standard error of the sum
Csu.sum$ci <- Csu.sum$se * qnorm(0.975) # calculate 95% confidence interval

# Autumn
Cau.sum <- aggregate(mean ~ group1, Ca.stat[19:27,], sum) # sum means across autumn months
aggregate(n ~ group1, Ca.stat[19:27,], sum)
Cau.sum$var <- aggregate(se^2 ~ group1, Ca.stat[19:27,], sum)[,2] # sum variance of the mean across autumn months
Cau.sum$se <- sqrt(Cau.sum$var) # calculate standard error of the sum
Cau.sum$ci <- Cau.sum$se * qnorm(0.975) # calculate 95% confidence interval

# Winter
Cwi.sum <- aggregate(mean ~ group1, Ca.stat[28:36,], sum) # sum means across winter months
aggregate(n ~ group1, Ca.stat[28:36,], sum)
Cwi.sum$var <- aggregate(se^2 ~ group1, Ca.stat[28:36,], sum)[,2] # sum variance of the mean across winter months
Cwi.sum$se <- sqrt(Cwi.sum$var) # calculate standard error of the sum
Cwi.sum$ci <- Cwi.sum$se * qnorm(0.975) # calculate 95% confidence interval

# Combined
Cs.sum <- data.frame(season = c(rep("Spring",3), rep("Summer",3), rep("Autumn",3), rep("Winter",3)),
                     rbind(Csp.sum, Csu.sum, Cau.sum, Cwi.sum))


#### 5.    Standard error of the product of export and density ####
#### 5.1   Species ####
#### 5.1.1 Biomass export ####
D <- read.csv("~/Desktop/Plymouth University/Dissertation/Publication/Data/Density.csv")
d.stat <- with(D, describeBy(density, species, mat = T, digits = 10)) # annual means and standard errors

d.stat$mean20[1] <- d.stat$mean[1] * 4/20
d.stat$mean20[2:3] <- d.stat$mean[2:3] * 16/20 
d.stat$se20[1] <- d.stat$se[1] * 4/20
d.stat$se20[2:3] <- d.stat$se[2:3] * 16/20 
# updated annual means and standard errors, assuming a 20-m kelp forest band 
# with 16 m of mixed L. ochroleuca and L. hyperborea and 4 m of L. digitata

ba.sum$se.prod <- sqrt(ba.sum$se^2 * d.stat$se20^2 + ba.sum$se^2 * d.stat$mean20^2 +
                       d.stat$se20^2 * ba.sum$mean^2) # combined standard error
ba.sum$ci.prod <- ba.sum$se.prod * qnorm(0.975) # combined 95% confidence interval
ba.sum

#### 5.1.2 Carbon export ####
Ca.sum$se.prod <- sqrt(Ca.sum$se^2 * d.stat$se20^2 + Ca.sum$se^2 * d.stat$mean20^2 +
                       d.stat$se20^2 * Ca.sum$mean^2) # combined standard error
Ca.sum$ci.prod <- Ca.sum$se.prod * qnorm(0.975) # combined 95% confidence interval
Ca.sum


#### 5.2   Season ####
#### 5.2.1 Biomass export ####
ds.stat <- with(D, describeBy(density, list(species, season), mat = T, digits = 10)) # annual means and standard errors

ds.stat$mean20[c(2:3,5:6,8:9,11:12)] <- ds.stat$mean[c(2:3,5:6,8:9,11:12)] * 16/20 
ds.stat$mean20[c(1,4,7,10)] <- ds.stat$mean[c(1,4,7,10)] * 4/20
ds.stat$se20[c(2:3,5:6,8:9,11:12)] <- ds.stat$se[c(2:3,5:6,8:9,11:12)] * 16/20 
ds.stat$se20[c(1,4,7,10)] <- ds.stat$se[c(1,4,7,10)] * 4/20
# updated annual means and standard errors, assuming a 20-m kelp forest band 
# with 16 m of mixed L. ochroleuca and L. hyperborea and 4 m of L. digitata
ds.stat <- rbind(ds.stat[4:9,], ds.stat[1:3,], ds.stat[10:12,]) # reorder to match export dataframe

bs.sum$se.prod <- sqrt(bs.sum$se^2 * ds.stat$se20^2 + bs.sum$se^2 * ds.stat$mean20^2 +
                       ds.stat$se20^2 * bs.sum$mean^2) # combined standard error
bs.sum$ci.prod <- bs.sum$se.prod * qnorm(0.975) # combined 95% confidence interval
bs.sum

#### 5.2.2 Carbon export ####
Cs.sum$se.prod <- sqrt(Cs.sum$se^2 * ds.stat$se20^2 + Cs.sum$se^2 * ds.stat$mean20^2 +
                       ds.stat$se20^2 * Cs.sum$mean^2) # combined standard error
Cs.sum$ci.prod <- Cs.sum$se.prod * qnorm(0.975) # combined 95% confidence interval
Cs.sum

#### 6.   Master data frame ####
consta <- data.frame(period = c(rep("Year",3),rep("Spring",3),rep("Summer",3),
                                rep("Autumn",3),rep("Winter",3)),
                     species = rep(c("d","h","o"),5),
                     biomass = c(ba.sum$mean, bs.sum$mean),
                     carbon = c(Ca.sum$mean, Cs.sum$mean),
                     density = c(d.stat$mean20, ds.stat$mean20),
                     bCI = c(ba.sum$ci.prod, bs.sum$ci.prod),
                     cCI = c(Ca.sum$ci.prod, Cs.sum$ci.prod))

write.csv(consta, "~/Desktop/Plymouth University/Dissertation/Publication/Data/Constants.csv", row.names = F)

#### 7.   Clean up ####
detach(package:fitdistrplus)
detach(package:car)
detach(package:nlme)
detach(package:psych)
detach(package:ggalt)
detach(package:ggplot2)
rm(list = ls())
graphics.off()
cat("\014")
