##################################################################
##### Project: Sequestration of carbon from Laminaria kelp   #####
#####          forests may be diminished in a warmer climate #####
##### Script purpose: Analysis of carbon assimilation        #####
##### Author: Luka Seamus Wright                             #####
##################################################################

#### 1.   Data preparation ####
#### 1.1  Load data ####
PP <- read.csv("~/Desktop/Plymouth University/Dissertation/Publication/Data/Assimilation.csv")

#### 1.2  Remove unusable data ####
PP <- PP[-c(19:72),]
# at age 13 days strong irradiance upon retrieval caused permanent damage
# to the photosynthetic apparatus
rownames(PP) <-  NULL

#### 1.3  Reorder levels of species factor ####
PP <- within(PP,{
      species <- factor(species, levels = c("o","h","d"))
})

#### 1.4  Rename  and convert variables ####
nppO2 <- PP$NPP # net primary production (μmol O2 gFW-1 h-1)
rO2 <- PP$R# respiration (μmol O2 gFW-1 h-1)
gppO2 <- PP$GPP # gross primary production (μmol O2 gFW-1 h-1)
d.w <- PP$d.w # dry mass to wet mass ratio

# convert μmol O2 gFW-1 to g C gDW-1 assuming a photosynthetic and respiratory quotient of 1
npp <- nppO2 / d.w * 1e-6 * 12.0107 # net primary production (gC gDW-1 h-1)
gpp <- gppO2 / d.w * 1e-6 * 12.0107 # gross primary production (gC gDW-1 h-1)
r <- rO2 / d.w * 1e-6 * 12.0107 # respiration (gC gDW-1 h-1)

age <- PP$age # tissue age in days
sp <- PP$species
bag <- PP$bag

#### 2.   Data analysis NPP ####
#### 2.1  Determine random components ####
require(lme4)
m1 <- lm(npp ~ sp * age) # fixed effects model
m2 <- lmer(npp ~ sp * age + (age|bag), REML = F) # mixed effects model
anova(m2, m1) # models are not different
# continue with m1

#### 2.2  Determine fixed components ####
m3 <- update(m1, .~. - sp : age) # remove interaction
anova(m1, m3) # m1 is the better model

#### 2.3  Test model fit ####
# assumption of homogeneity
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m1, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m1) ~ age) # residual variance only varies slightly with age
boxplot(resid(m1) ~ sp) # and barely with species
# overall quite homogenous

# assumption of normality
hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1))
# overall quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m1 is chosen as the optimal model

#### 2.4  Interpret model ####
# Laminaria ochroleuca
require(car)
Anova(m1, type = 3) # Type III Sums of Squares test
# Response: npp
#                 Sum Sq  Df F value    Pr(>F)    
# (Intercept) 1.9876e-06   1 19.0347 2.733e-05 ***
# sp          8.2580e-07   2  3.9544 0.0217197 *  
# age         3.8397e-06   1 36.7720 1.585e-08 ***
# sp:age      1.9362e-06   2  9.2712 0.0001803 ***
# Residuals   1.2530e-05 120  

# functions, slopes and pairwise contrasts of intercepts and slopes
summary(m1)
# L. och. intercept vs. L. hyp. intercept, t = 2.344, p = 0.02 *
# L. och. intercept vs. L. dig. intercept, t = 2.518, p = 0.01 *
# L. och. slope vs. L. hyp. slope, t = 2.992, p = 0.003 **
# L. och. slope vs. L. dig. slope, t = 4.178, p < 0.001 ***

coef(m1)
# L. och. function: f(x) = -2.883424e-05x + 5.514990e-04
-5.514990e-04/-2.883424e-05 # after 19.12653 days L. och. respiration outweighs
# photosynthesis at 50.4 μmol photons m-2 s-1

# Laminaria hyperborea
sp <- factor(sp, levels = c("h", "d", "o"))
m1 <- lm(npp ~ sp * age)
Anova(m1, type = 3)
# Response: npp
#                 Sum Sq  Df F value    Pr(>F)    
# (Intercept) 6.1544e-06   1 58.9389 4.833e-12 ***
# sp          8.2580e-07   2  3.9544 0.0217197 *  
# age         3.5060e-07   1  3.3573 0.0693880 .  
# sp:age      1.9362e-06   2  9.2712 0.0001803 ***
# Residuals   1.2530e-05 120    

summary(m1)
# L. hyp. intercept vs. L. dig. intercept, t = 0.174, p = 0.86
# L. hyp. slope vs. L. dig. slope, t = 1.186, p = 0.24

coef(m1)
# L. hyp. function: f(x) = -8.712518e-06x + 9.704492e-04
-9.704492e-04/-8.712518e-06 # after 111.3856 days L. hyp. respiration outweighs
# photosynthesis at 50.4 μmol photons m-2 s-1

# Laminaria digitata
sp <- factor(sp, levels = c("d", "h", "o"))
m1 <- lm(npp ~ sp * age)
Anova(m1, type = 3)
# Response: npp
#                 Sum Sq  Df F value    Pr(>F)    
# (Intercept) 6.5563e-06   1 62.7883 1.321e-12 ***
# sp          8.2580e-07   2  3.9544 0.0217197 *  
# age         2.5000e-09   1  0.0242 0.8765794    
# sp:age      1.9362e-06   2  9.2712 0.0001803 ***
# Residuals   1.2530e-05 120  

coef(m1)
# L. dig. function: f(x) = -7.400604e-07x + 1.001638e-03
-1.001638e-03/-7.400604e-07 # after 1353.454 days L. dig. respiration outweighs
# photosynthesis at 50.4 μmol photons m-2 s-1

sp <- factor(sp, levels = c("o", "h", "d"))


#### 3.   Data analysis R ####
#### 3.1  Determine random components ####
m4 <- lm(r ~ sp * age) # fixed effects model
m5 <- lmer(r ~ sp * age + (age|bag), REML = F) # mixed effects model
anova(m5, m4) # m4 fits better

#### 3.2  Determine fixed components ####
m6 <- update(m4, .~. - sp : age) # remove interaction
anova(m6, m4) # m4 is the better model

#### 3.3  Test model fit ####
# assumption of homogeneity
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m4, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m4) ~ age) # residual variance only varies slightly with age
boxplot(resid(m4) ~ sp) # but dramatically with species
# homogeneity could be improved

# assumption of normality
hist(resid(m4))
qqnorm(resid(m4))
qqline(resid(m4))
# overall quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# optimise model with weighted least squares

#### 3.4  Determine random components ####
require(nlme)
m7 <- gls(r ~ sp * age, weights = varIdent(form = ~1|sp)) # fixed effects model
m8 <- lme(r ~ sp * age, random = ~age|bag, 
          weights = varIdent(form = ~1|sp)) # mixed effects model
anova(m8, m7) # m7 fits better

#### 3.5  Determine fixed components ####
m7 <- gls(r ~ sp * age, weights = varIdent(form = ~1|sp),
          method = "ML")
m9 <- update(m7, .~. - sp : age) # remove interaction
anova(m9, m7) # m7 is the better model
m7 <- gls(r ~ sp * age, weights = varIdent(form = ~1|sp),
          method = "REML")

#### 3.6  Test model fit ####
# assumption of homogeneity
plot(m7, col = sp)
plot(resid(m7, type = "normalized") ~ sp)
# residual spread no longer varies between species
# homogeneity is improved

# assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m7))
qqnorm(resid(m7))
qqline(resid(m7))
# overall quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m7 is chosen as the optimal model

#### 3.7  Interpret model ####
# Laminaria ochroleuca
Anova(m7, type = 3) # Type III Sums of Squares test
# Response: r
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 32.4989  1.193e-08 ***
# sp           2 42.6757  5.409e-10 ***
# age          1  1.5904     0.2073    
# sp:age       2 25.5263  2.864e-06 ***

summary(m7)
# L. och. intercept vs. L. hyp. intercept, t = 0.569637, p = 0.57
# L. och. intercept vs. L. dig. intercept, t = -3.385381, p = 0.001 **
# L. och. slope vs. L. hyp. slope, t = -0.335553, p = 0.74
# L. och. slope vs. L. dig. slope, t = 2.707743, p = 0.008 **

coef(m7)
# L. och. function: f(x) = -5.932082e-06x + 7.128592e-04
-7.128592e-04/-5.932082e-06 # after 120.1702 days L. och. respiration is zero

# Laminaria hyperborea
sp <- factor(sp, levels = c("h", "d", "o"))
m7 <- gls(r ~ sp * age, weights = varIdent(form = ~1|sp),
          method = "REML")
Anova(m7, type = 3)
# Response: r
#             Df    Chisq Pr(>Chisq)    
# (Intercept)  1 108.5349  < 2.2e-16 ***
# sp           2  42.6757  5.409e-10 ***
# age          1   7.3247   0.006801 ** 
# sp:age       2  25.5263  2.864e-06 ***

summary(m7)
# L. hyp. intercept vs. L. dig. intercept, t = -6.055639, p < 0.001 ***
# L. hyp. slope vs. L. dig. slope, t = 4.642899, p < 0.001 ***

coef(m7)
# L. hyp. function: f(x) = -7.781989e-06x + 7.963444e-04
-7.963444e-04/-7.781989e-06 # after 102.3317 days L. hyp. respiration is zero

# Laminaria digitata
sp <- factor(sp, levels = c("d", "h", "o"))
m7 <- gls(r ~ sp * age, weights = varIdent(form = ~1|sp),
          method = "REML")
Anova(m7, type = 3)
# Response: r
#             Df  Chisq Pr(>Chisq)    
# (Intercept)  1 37.992  7.102e-10 ***
# sp           2 42.676  5.409e-10 ***
# age          1 21.692  3.201e-06 ***
# sp:age       2 25.526  2.864e-06 ***

coef(m7)
# L. dig. function: f(x) = 7.537290e-06x + 2.651775e-04
# L. dig. respiration never reaches zero

sp <- factor(sp, levels = c("o", "h", "d"))


#### 4.   Data analysis GPP ####
#### 4.1  Determine random components ####
m10 <- lm(gpp ~ sp * age) # fixed effects model
m11 <- lmer(gpp ~ sp * age + (age|bag), REML = F) # mixed effects model
anova(m11, m10) # m10 fits better

#### 4.2  Determine fixed components ####
m12 <- update(m10, .~. - sp : age) # remove interaction
anova(m12, m10) # m10 is the better model

#### 4.3  Test model fit ####
# assumption of homogeneity
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m10, col = sp) 
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m10) ~ age) # residual variance barely varies with age
boxplot(resid(m10) ~ sp) # but varies slightly with species
# overall quite homogenous

# assumption of normality
hist(resid(m10))
qqnorm(resid(m10))
qqline(resid(m10)) # normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m10 is chosen as the optimal model

#### 4.4  Interpret model ####
# Laminaria ochroleuca
Anova(m10, type = 3) # Type III Sums of Squares test
# Response: gpp
#                 Sum Sq  Df F value    Pr(>F)    
# (Intercept) 1.0447e-05   1 71.0057 9.090e-14 ***
# sp          1.0944e-06   2  3.7194   0.02709 *  
# age         5.5821e-06   1 37.9415 1.004e-08 ***
# sp:age      4.0085e-06   2 13.6229 4.659e-06 ***
# Residuals   1.7655e-05 120 

summary(m10)
# L. och. intercept vs. L. hyp. intercept, t = 2.368, p = 0.02 *
# L. och. intercept vs. L. dig. intercept, t = 0.012, p = 0.99
# L. och. slope vs. L. hyp. slope, t = 2.289, p = 0.02 *
# L. och. slope vs. L. dig. slope, t = 5.207, p < 0.001 ***

coef(m10)
# L. och. function: f(x) = -3.476632e-05x + 1.264358e-03
-1.264358e-03/-3.476632e-05 # at 36.36732 days the photosynthetic apparatus of
# L. och. shuts down at 50.4 μmol photons m-2 s-1

# Laminaria hyperborea
sp <- factor(sp, levels = c("h", "d", "o"))
m10 <- lm(gpp ~ sp * age)
Anova(m10, type = 3)
# Response: gpp
#                 Sum Sq  Df  F value    Pr(>F)    
# (Intercept) 2.0399e-05   1 138.6516 < 2.2e-16 ***
# sp          1.0944e-06   2   3.7194  0.027088 *  
# age         1.2565e-06   1   8.5404  0.004152 ** 
# sp:age      4.0085e-06   2  13.6229 4.659e-06 ***
# Residuals   1.7655e-05 120   

summary(m10)
# L. hyp. intercept vs. L. dig. intercept, t = -2.356, p = 0.02 *
# L. hyp. slope vs. L. dig. slope, t = 2.918, p = 0.004 **

coef(m10)
# L. hyp. function: f(x) = -1.649451e-05x + 1.766794e-03
-1.766794e-03/-1.649451e-05 # at 107.1141 days the photosynthetic apparatus of
# L. hyp. shuts down at 50.4 μmol photons m-2 s-1

# Laminaria digitata
sp <- factor(sp, levels = c("d", "h", "o"))
m10 <- lm(gpp ~ sp * age)
Anova(m10, type = 3)
# Response: gpp
#                 Sum Sq  Df F value    Pr(>F)    
# (Intercept) 1.0487e-05   1 71.2821 8.326e-14 ***
# sp          1.0944e-06   2  3.7194   0.02709 *  
# age         2.1340e-07   1  1.4503   0.23085    
# sp:age      4.0085e-06   2 13.6229 4.659e-06 ***
# Residuals   1.7655e-05 120 

coef(m10)
# L. dig. function: f(x) = 6.797230e-06x + 1.266816e-03
# due to the positive slope, the photosynthetic apparatus of L. dig.
# apparently never shuts down at 50.4 μmol photons m-2 s-1

#### 5.   Data visualisation ####
#### 5.1  Descriptive statistics ####
require(psych)
npp.stat <- describeBy(npp, list(sp, age), mat = T)
npp.stat$group2 <- as.integer(npp.stat$group2) # turn age factor back into integer

r.stat <- describeBy(r, list(sp, age), mat = T)
r.stat$group2 <- as.integer(r.stat$group2) # turn age factor back into integer

gpp.stat <- describeBy(gpp, list(sp, age), mat = T)
gpp.stat$group2 <- as.integer(gpp.stat$group2) # turn age factor back into integer

#### 5.2  Model predictions ####
new <- data.frame(age = rep(0:32, 3),
                  sp = c(rep("d", 33), rep("h", 33), rep("o", 33)))


fit <- data.frame(predict(m1, interval = "confidence",
                          newdata = new))
new$npp.fit <- fit$fit
new$npp.lo <- fit$lwr
new$npp.hi <- fit$upr

new$r.fit <- predict(m7, newdata = new)
modmat <-  model.matrix(formula(m7)[-2], new)
int <- diag(modmat %*% vcov(m7) %*% t(modmat))
new$r.lo <- with(new, r.fit - qnorm(0.975)*sqrt(int))
new$r.hi <- with(new, r.fit + qnorm(0.975)*sqrt(int))

fit <- data.frame(predict(m10, interval = "confidence",
                          newdata = new))
new$gpp.fit <- fit$fit
new$gpp.lo <- fit$lwr
new$gpp.hi <- fit$upr

#### 5.3  Customised theme ####
require(ggplot2)
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .3, .2, .2),"cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title = element_blank(),
                 text = element_text(family = "Helvetica Neue"))

#### 5.4  Plots ####
# NB number is g C are too small -> convert to mg C for plotting
np <- ggplot() +
        geom_hline(yintercept = 0, lwd = .5, colour = "#000000") +
        geom_line(data = new, aes(age, npp.fit*1000, colour = sp, lty = sp), size = 0.5) +
        geom_ribbon(data = new, aes(age, ymin = npp.lo*1000, ymax = npp.hi*1000, fill = sp),
                    alpha = 0.5) +
        geom_pointrange(data = npp.stat, aes(group2, mean*1000, ymin = mean*1000 - se*1000,
                                             ymax = mean*1000 + se*1000, colour = group1),
                        size = 0.5, position = position_dodge(width = 1.5)) +
        scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")*"        y = –0.001x + 1"),
                                       expression(italic("L. hyperborea")*"  y = –0.009x + 0.97"),
                                       expression(italic("L. ochroleuca")*"  y = –0.029x + 0.55")),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#333b08","#627d0e","#f1c700"),
                          labels = c(expression(italic("L. digitata")*"        y = –0.001x + 1"),
                                     expression(italic("L. hyperborea")*"  y = –0.009x + 0.97"),
                                     expression(italic("L. ochroleuca")*"  y = –0.029x + 0.55")),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(5, 5, 1),
                              guide = "none") +
        ylab(expression("Carbon assimilation (mg C g"^-1*" h"^-1*")")) +
        xlab("Detrital age (d)") +
        scale_x_continuous(breaks = seq(0, 35, by = 5)) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(-0.5, 2), xlim = c(0.5, 33.36)) +
        theme(legend.position = c(0.41, 0.905)) +
        mytheme

np # dimsenisons: 4 x 5 in

rp <- ggplot() +
        geom_line(data = new, aes(age, r.fit*1000, colour = sp, lty = sp), size = 0.5) +
        geom_ribbon(data = new, aes(age, ymin = r.lo*1000, ymax = r.hi*1000, fill = sp), 
                    alpha = .5) +
        geom_pointrange(data = r.stat, aes(group2, mean*1000, ymin = mean*1000 - se*1000,
                                           ymax = mean*1000 + se*1000, colour = group1),
                        size = 0.5, position = position_dodge(width = 1.5)) +
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
        scale_linetype_manual(values = c(1, 1, 5),
                              guide = "none") +
        ylab(expression("Carbon emission (mg C g"^-1*" h"^-1*")")) +
        xlab("Detrital age (d)") +
        scale_x_continuous(breaks = seq(0, 35, by = 5)) +
        scale_y_continuous(breaks = seq(0, 1.2, by = 0.3), expand = c(0, 0)) +
        coord_cartesian(ylim = c(0, 1.2), xlim = c(0, 33.34)) +
        theme(legend.position = c(0.89, 0.905)) +
        mytheme

rp # dimsenisons: 4 x 5 in

gp <- ggplot() +
        geom_line(data = new, aes(age, gpp.fit*1000, colour = sp, lty = sp), size = 0.5) +
        geom_ribbon(data = new, aes(age, ymin = gpp.lo*1000, ymax = gpp.hi*1000, fill = sp),
                    alpha = 0.5) +
        geom_pointrange(data = gpp.stat, aes(group2, mean*1000, ymin = mean*1000 - se*1000,
                                             ymax = mean*1000 + se*1000, colour = group1),
                        size = 0.5, position = position_dodge(width = 1.5)) +
        scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")*"        y = 0.007x + 1.27"),
                                       expression(italic("L. hyperborea")*"  y = –0.016x + 1.77"),
                                       expression(italic("L. ochroleuca")*"  y = –0.035x + 1.26")),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#333b08","#627d0e","#f1c700"),
                          labels = c(expression(italic("L. digitata")*"        y = 0.007x + 1.27"),
                                     expression(italic("L. hyperborea")*"  y = –0.016x + 1.77"),
                                     expression(italic("L. ochroleuca")*"  y = –0.035x + 1.26")),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(5, 1, 1),
                              guide = "none") +
        ylab(expression("Carbon assimilation (mg C g"^-1*" h"^-1*")")) +
        xlab("Detrital age (d)") +
        scale_x_continuous(breaks = seq(0, 35, by = 5)) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(0, 2.5), xlim = c(0, 33.34)) +
        theme(legend.position = c(0.66, 0.905)) +
        mytheme

gp # dimsenisons: 4 x 5 in

#### 5.5  Combined ####
require(cowplot)
# plot pp needs to be loaded from the Decomposition.R file
final <- plot_grid(np, pp, rel_widths = c(1.4, 1), 
                   labels = "auto", hjust = 0.5)
final # dimensions 4 x 8.5 in

#### 6.   Clean up ####
detach(package:fitdistrplus)
detach(package:car)
detach(package:lme4)
detach(package:nlme)
detach(package:psych)
detach(package:cowplot)
detach(package:ggplot2)
rm(list = ls())
graphics.off()
cat("\014")
