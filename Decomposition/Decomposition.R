##################################################################
##### Project: Sequestration of carbon from Laminaria kelp   #####
#####          forests may be diminished in a warmer climate #####
##### Script purpose: Analysis of decomposition, phenols and #####
#####                 elemental stoichiometry                #####
##### Author: Luka Seamus Wright                             #####
##################################################################

#### 1.   Data preparation ####
#### 1.1  Load data ####
Deco <- read.csv("~/PATH/Decomposition.csv")
Biochem <- read.csv("~/PATH/Biochemical.csv")
Grazing <- read.csv("~/PATH/Grazing.csv")

#### 1.2  Reorder levels of species factor ####
Deco <- within(Deco,{
  species <- factor(species, levels = c("o","h","d"))
})

Biochem <- within(Biochem,{
  species <- factor(species, levels = c("o","h","d"))
})


Grazing <- within(Grazing,{
  species <- factor(species, levels = c("o","h","d"))
})

#### 1.3  Rename variables ####
loss <- Biochem$perc.loss # biomass loss (% d-1)
phen <- Biochem$phenols # phenols (%)
C <- Biochem$C # carbon (%)
N <- Biochem$N # nitrogen (%)
CN <- Biochem$CN # carbon:nitrogen
sp <- Biochem$species
age <- Biochem$age
bag <- Biochem$bag

ex <- Grazing$excavation # excavated tissue (%)
per <- Grazing$perforation # perforated tissue (%)
sp2 <- Grazing$species
bag2 <- Grazing$bag

loss2 <- Deco$perc.loss # biomass loss (% d-1)
mesh <- Deco$mesh # mesh aperture (mm)
site <- Deco$site
subs <- Deco$substratum
sp3 <- Deco$species

#### 2.   Decomposition data analysis ####
#### 2.1  Determine fixed components ####
m1 <- lm(loss2 ~ sp3 + subs + mesh + site)
drop1(m1, test = "F")
m1 <- lm(loss2 ~ sp3 + subs + mesh)
drop1(m1, test = "F")
m1 <- lm(loss2 ~ sp3 * subs)
drop1(m1, test = "F")
m1 <- lm(loss2 ~ sp3 + subs)

#### 2.2  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m1, col = sp3)
par(mfrow = c(1,1), mar = c(2,2,2,1))
boxplot(resid(m1) ~ sp3)
# overall very homogenous

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1))
# overall perfectly normal

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m1 is chosen as the optimal model based on its good homogeneity
# and perfect normality

#### 2.3  Interpret model ####
require(car)
Anova(m1, type = 2) # Type II Sums of Squares test
# Response: loss2
#            Sum Sq  Df F value   Pr(>F)   
# sp3        19.418   2  7.1261 0.001173 **
# subs        7.191   1  5.2780 0.023259 * 
# Residuals 170.307 125 

summary(m1)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.7207     0.1948   8.833 7.76e-15 ***
# sp3h         -0.9357     0.2517  -3.717 0.000303 ***
# sp3d         -0.6119     0.2517  -2.431 0.016487 *  
# subsSand     -0.4885     0.2126  -2.297 0.023259 * 

sp3 <- factor(sp3, levels = c("h", "d", "o"))
m1 <- lm(loss2 ~ sp3 + subs)
summary(m1)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.7850     0.1948   4.030 9.63e-05 ***
# sp3d          0.3238     0.2517   1.286 0.200759    
# sp3o          0.9357     0.2517   3.717 0.000303 ***
# subsSand     -0.4885     0.2126  -2.297 0.023259 *  

sp3 <- factor(sp3, levels = c("o", "h", "d"))

#### 3.   Phenol data analysis ####
#### 3.1  Determine random components ####
require(lme4)
m2 <- lm(phen ~ sp * age) # fixed effects model
m3 <- lmer(phen ~ sp * age + (age|bag), REML = F) # mixed effects model
# singular fit -> simplify mixed model
m3 <- lmer(phen ~ sp * age + (1|bag), REML = F) # mixed effects model
# still too complex
anova(m3, m2) # continue with m2

#### 3.2  Determine fixed components ####
drop1(m2, test = "F")
# keep interaction term

#### 3.3  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m2, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m2) ~ age) # residual variance is quite constant across ages
plot(resid(m2) ~ sp) # but is strongly variable between species
# overall heterogenous

hist(resid(m2))
qqnorm(resid(m2))
qqline(resid(m2))
# residuals deviate from normality at distribution edges

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# mostly homogeneity needs improving
# -> model heterogeneity with generalised least squares

#### 3.4  Determine random components ####
require(nlme)
m4 <- gls(phen ~ sp * age,
          weights = varIdent(form = ~1|sp),
          method = "REML") # fixed effects model
m5 <- lme(phen ~ sp * age,
          random = ~age|bag,
          weights = varIdent(form = ~1|sp),
          method = "REML") # mixed effects model
anova(m5, m4) # m4 is better
m4 <- gls(phen ~ sp * age,
          weights = varIdent(form = ~1|sp),
          method = "ML")

#### 3.5  Determine fixed components ####
drop1(m4, test = "Chisq")
# keep interaction term

m4 <- gls(phen ~ sp * age,
          weights = varIdent(form = ~1|sp),
          method = "REML")

#### 3.6  Test model fit ####
plot(m4, col = sp)
boxplot(resid(m4, type = "normalized") ~ sp)
# homogeneity improved

par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m4))
qqnorm(resid(m4))
qqline(resid(m4))
# normality is the same as m2

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m4 is chosen as the optimal model on the basis of good homogeneity
# and a balanced, although not perfectly normal, data distribution

#### 3.7  Interpret model ####
Anova(m4, type = 3) # Type III Sums of Squares test
# Response: phen
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1  2.0304  0.1541827    
# sp           2  4.8453  0.0886882 .  
# age          1  1.3975  0.2371457    
# sp:age       2 15.2352  0.0004917 ***

summary(m4)
# L. och. intercept vs. L. hyp. intercept, t = 2.068532, p = 0.04 *
# L. och. intercept vs. L. dig. intercept, t = 0.022617, p = 0.98
# L. och. slope vs. L. hyp. slope, t = 3.264467, p = 0.002 **
# L. och. slope vs. L. dig. slope, t = -0.855525, p = 0.4

coef(m4)
# L. och. function: 0.003021621x + 0.089658453

sp <- factor(sp, levels = c("h", "d", "o"))
m4 <- gls(phen ~ sp * age,
          weights = varIdent(form = ~1|sp),
          method = "REML")
Anova(m4, type = 3)
# Response: phen
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1  7.6066  0.0058156 ** 
# sp           2  4.8453  0.0886882 .  
# age          1 15.5665  7.966e-05 ***
# sp:age       2 15.2352  0.0004917 ***
  
summary(m4)
# L. hyp. intercept vs. L. dig. intercept, t = -2.200702, p = 0.03 *
# L. hyp. slope vs. L. dig. slope, t = -3.817406, p < 0.001 ***

coef(m4)
# L. hyp. function: 0.02641994x + 0.45463995

sp <- factor(sp, levels = c("d", "h", "o"))
m4 <- gls(phen ~ sp * age,
          weights = varIdent(form = ~1|sp),
          method = "REML")
Anova(m4, type = 3)
# Response: phen
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 72.3064  < 2.2e-16 ***
# sp           2  4.8453  0.0886882 .  
# age          1  3.4076  0.0648957 .  
# sp:age       2 15.2352  0.0004917 ***

coef(m4)
# L. dig. function: f(x) = 0.0008033959x + 0.0911020732

sp <- factor(sp, levels = c("o", "h", "d"))


#### 4.   Decomposition~phenol data analysis ####
#### 4.1  Determine random components ####
m6 <- lm(loss ~ phen * sp) # fixed effects model
m7 <- lmer(loss ~ phen * sp + (phen|bag), REML = F) # mixed effects model
# singular fit -> simplify mixed effects model
m7 <- lmer(loss ~ phen * sp + (1|bag), REML = F) # mixed effects model
# still too complex
anova(m7, m6) # continue with m6

#### 4.2  Determine fixed components ####
drop1(m6, test = "F") # retain interaction

#### 4.3  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m6, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m6) ~ age) # residual spread increases with age
plot(resid(m6) ~ sp) # and varies between species
# homogeneity could be improved

hist(resid(m6))
qqnorm(resid(m6))
qqline(resid(m6))
# almost perfect normality

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# try generalised least squares

#### 4.4  Determine random components ####
m8 <- gls(loss ~ phen * sp, 
          weights = varIdent(form = ~1|sp)) # fixed effects model
m9 <- lme(loss ~ phen * sp,
          random = ~1|bag, # simplify because no convergence
          weights = varIdent(form = ~1|sp)) # mixed effects model
anova(m9, m8) # continue with m8

#### 4.5  Determine fixed components ####
m8 <- gls(loss ~ phen * sp, 
          weights = varIdent(form = ~1|sp),
          method = "ML")
drop1(m8, test = "Chisq") # retain interaction
m8 <- gls(loss ~ phen * sp, 
          weights = varIdent(form = ~1|sp),
          method = "REML")

#### 4.6  Test model fit ####
plot(m8, col = sp)
boxplot(resid(m8, type = "normalized") ~ sp) # no large variance between species
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m8))
qqnorm(resid(m8))
qqline(resid(m8))
# almost perfect normality

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m8 is chosen as the optimal model

#### 4.7  Interpret model ####
Anova(m8, type = 3) # Type III Sums of Squares test
# Response: loss
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 37.0270  1.165e-09 ***
# phen         1 10.5494   0.001162 ** 
# sp           2  0.0171   0.991486    
# phen:sp      2  8.1594   0.016913 * 

summary(m8)
#                  Value Std.Error   t-value p-value
# (Intercept)   3.293228  0.541206  6.084979  0.0000
# phen         -9.230244  2.841839 -3.247982  0.0017
# sph           0.104560  0.805774  0.129764  0.8971
# spd           0.025830  1.365323  0.018919  0.9850
# phen:sph      6.791463  2.891203  2.349010  0.0215
# phen:spd    -11.754796 11.610703 -1.012410  0.3146

coef(m8)
# L. och. function: f(x) = -9.2302443x + 3.2932281


sp <- factor(sp, levels = c("h", "d", "o"))
m8 <- gls(loss ~ phen * sp, 
          weights = varIdent(form = ~1|sp),
          method = "REML")
Anova(m8, type = 3) 
# Response: loss
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 32.3962  1.257e-08 ***
# phen         1 21.0162  4.554e-06 ***
# sp           2  0.0171    0.99149    
# phen:sp      2  8.1594    0.01691 * 
  
summary(m8)
#                  Value Std.Error   t-value p-value
# (Intercept)   3.397788  0.596965  5.691768  0.0000
# phen         -2.438781  0.531980 -4.584348  0.0000
# spd          -0.078730  1.388369 -0.056707  0.9549
# spo          -0.104560  0.805774 -0.129764  0.8971
# phen:spd    -18.546259 11.270110 -1.645615  0.1040
# phen:spo     -6.791463  2.891203 -2.349010  0.0215

coef(m8)
# L. hyp. function: f(x) = -2.43878127x + 3.39778831

sp <- factor(sp, levels = c("d", "h", "o"))
m8 <- gls(loss ~ phen * sp, 
          weights = varIdent(form = ~1|sp),
          method = "REML")
Anova(m8, type = 3) 
# Response: loss
#             Df  Chisq Pr(>Chisq)   
# (Intercept)  1 7.0113    0.00810 **
# phen         1 3.4748    0.06231 . 
# sp           2 0.0171    0.99149   
# phen:sp      2 8.1594    0.01691 * 

coef(m8)
# L. dig. function: f(x) = -20.98503999x + 3.31905836

sp <- factor(sp, levels = c("o", "h", "d"))

# This model shows that decomposition speed correlates with phenolic
# content intraspecifically. Now it would be interesting to see if 
# decomposition correlates with phenolic content interspecifically.

#### 4.8  Determine random components ####
m10 <- lm(loss ~ phen) # fixed effects model
m11 <- lmer(loss ~ phen + (phen|bag), REML = F) # mixed effects model
m11 <- lmer(loss ~ phen + (1|bag), REML = F)
anova(m11, m10) # continue with m10

#### 4.9  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m10, col = sp)
par(mfrow = c(1,1), mar = c(2,2,2,1))
plot(resid(m10) ~ phen) 
# residual variance decreases with increasing phenolic content
# homogeneity could be improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m10))
qqnorm(resid(m10))
qqline(resid(m10))
# perfect normality

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity can be improved with generalised least squares

#### 4.10  Determine random components ####
m12 <- gls(loss ~ phen, 
           weights = varPower(form = ~phen)) # fixed effects model
m13 <- lme(loss ~ phen,
           random = ~1|bag,
           weights = varPower(form = ~phen)) # mixed effects model
anova(m13, m12) # continue with m12

#### 4.11  Test model fit ####
plot(m12, col = sp)
plot(resid(m12, type = "normalized") ~ phen) # decrease in variance improved
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m12))
qqnorm(resid(m12))
qqline(resid(m12))
# perfect normality

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m12 is chosen as the optimal model

#### 4.12  Interpret model ####
Anova(m12, type = 2) # Type II Sums of Squares test
# Response: loss
#      Df  Chisq Pr(>Chisq)    
# phen  1 15.632  7.695e-05 ***

coef(m12) 
# y = -1.048487x + 1.662985


#### 5.   Elemental data analysis ####
#### 5.1  Determine random components ####
m14 <- lm(C ~ sp * age) # fixed effects model
m15 <- lmer(C ~ sp * age + (age|bag), REML = F) # mixed effects model
# singular fit
m15 <- lmer(C ~ sp * age + (1|bag), REML = F) # mixed effects model
anova(m15, m14) # models are not different
# continue with m14

#### 5.2  Determine fixed components ####
drop1(m14, test = "F") # retain interaction term

#### 5.3  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m14, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m14) ~ age) # residual variance increases with age
boxplot(resid(m14) ~ sp) # and varies between species
# overall heterogeous

hist(resid(m14))
qqnorm(resid(m14))
qqline(resid(m14))
# deviation from normality at distribution edges

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity needs improving
# -> model heterogeneity with generalised least squares

#### 5.4  Determine random components ####
m16 <- gls(C ~ sp * age,
           weights = varExp(form = ~age|sp),
           method = "REML") # fixed effects model
m17 <- lme(C ~ sp * age,
           random = ~age|bag,
           weights = varExp(form = ~age|sp),
           method = "REML") # mixed effects model
anova(m17, m16) # m17 is a better fit

m17 <- lme(C ~ sp * age,
           random = ~age|bag,
           weights = varExp(form = ~age|sp),
           method = "ML")

#### 5.5  Determine fixed components ####
drop1(m17, test = "Chisq") # retain interaction term

m17 <- lme(C ~ sp * age,
           random = ~age|bag,
           weights = varExp(form = ~age|sp))

#### 5.6  Test model fit ####
plot(m17, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m17, type = "normalized") ~ age)
boxplot(resid(m17, type = "normalized") ~ sp)
# homogeneity is perfect

hist(resid(m17))
qqnorm(resid(m17))
qqline(resid(m17))
# normality similar to m14

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m17 is chosen as the optimal model because of good homogeneity
# and a balanced data distribution albeit deviation at the distribution edges

#### 5.7  Interpret model ####
Anova(m17, type = 3) # Type III Sums of Squares test
# Response: C
#                Chisq Df Pr(>Chisq)    
# (Intercept) 124.0719  1  < 2.2e-16 ***
# sp            4.4934  2   0.105747    
# age           2.5079  1   0.113279    
# sp:age        9.9751  2   0.006822 ** 

summary(m17)
#                 Value Std.Error DF   t-value p-value
# (Intercept) 26.400310 2.3701309 69 11.138756  0.0000
# sph          2.055624 2.7290956  6  0.753225  0.4798
# spd         -1.633484 2.6113174  6 -0.625540  0.5546
# age         -0.186600 0.1178308 69 -1.583625  0.1179
# sph:age      0.195964 0.1277298 69  1.534207  0.1296
# spd:age      0.322677 0.1224947 69  2.634210  0.0104

# y = -0.186600x + 26.400310

coef(m17)
#    (Intercept)      sph       spd        age   sph:age   spd:age
# B1    27.25034 2.055624 -1.633484 -0.1865998 0.1959639 0.3226769
# B2    26.37077 2.055624 -1.633484 -0.1865998 0.1959639 0.3226769
# B3    25.57983 2.055624 -1.633484 -0.1865998 0.1959639 0.3226769
# B4    27.29081 2.055624 -1.633484 -0.1865998 0.1959639 0.3226769
# B5    26.95732 2.055624 -1.633484 -0.1865998 0.1959639 0.3226769
# B6    24.95280 2.055624 -1.633484 -0.1865998 0.1959639 0.3226769
# B7    27.11138 2.055624 -1.633484 -0.1865998 0.1959639 0.3226769
# B8    25.76476 2.055624 -1.633484 -0.1865998 0.1959639 0.3226769
# B9    26.32479 2.055624 -1.633484 -0.1865998 0.1959639 0.3226769

sp <- factor(sp, levels = c("h", "d", "o"))
m17 <- lme(C ~ sp * age,
           random = ~age|bag,
           weights = varExp(form = ~age|sp))
Anova(m17, type = 3)
# Response: C
#                Chisq Df Pr(>Chisq)    
# (Intercept) 442.3741  1  < 2.2e-16 ***
# sp            4.4934  2   0.105747    
# age           0.0361  1   0.849365    
# sp:age        9.9751  2   0.006822 ** 

summary(m17)
#                 Value Std.Error DF   t-value p-value
# (Intercept) 28.455934 1.3529384 69 21.032690  0.0000
# spd         -3.689108 1.7412352  6 -2.118673  0.0784
# spo         -2.055624 2.7290956  6 -0.753225  0.4798
# age          0.009364 0.0493031 69  0.189929  0.8499
# spd:age      0.126713 0.0595957 69  2.126210  0.0371
# spo:age     -0.195964 0.1277298 69 -1.534207  0.1296

# y = 0.009364x + 28.455934

coef(m17)
#    (Intercept)       spd       spo         age  spd:age    spo:age
# B1    29.30596 -3.689108 -2.055624 0.009364081 0.126713 -0.1959639
# B2    28.42639 -3.689108 -2.055624 0.009364082 0.126713 -0.1959639
# B3    27.63545 -3.689108 -2.055624 0.009364083 0.126713 -0.1959639
# B4    29.34644 -3.689108 -2.055624 0.009364080 0.126713 -0.1959639
# B5    29.01294 -3.689108 -2.055624 0.009364083 0.126713 -0.1959639
# B6    27.00842 -3.689108 -2.055624 0.009364083 0.126713 -0.1959639
# B7    29.16701 -3.689108 -2.055624 0.009364081 0.126713 -0.1959639
# B8    27.82038 -3.689108 -2.055624 0.009364083 0.126713 -0.1959639
# B9    28.38041 -3.689108 -2.055624 0.009364082 0.126713 -0.1959639

sp <- factor(sp, levels = c("d", "h", "o"))
m17 <- lme(C ~ sp * age,
           random = ~age|bag,
           weights = varExp(form = ~age|sp))
Anova(m17, type = 3)
# Response: C
#                Chisq Df Pr(>Chisq)    
# (Intercept) 510.5428  1  < 2.2e-16 ***
# sp            4.4934  2   0.105747    
# age          16.5204  1  4.813e-05 ***
# sp:age        9.9751  2   0.006822 ** 

summary(m17)
# y = 0.136077x + 24.766826

sp <- factor(sp, levels = c("o", "h", "d"))


#### 5.8  Determine random components ####
m18 <- lm(N ~ sp * age) # fixed effects model
m19 <- lmer(N ~ sp * age + (age|bag), REML = F) # mixed effects model
# singular fit
m19 <- lmer(N ~ sp * age + (1|bag), REML = F) # mixed effects model
anova(m19, m18) # models are not different
# continue with m18

#### 5.9  Determine fixed components ####
drop1(m18, test = "F") # retain interaction term

#### 5.10  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m18, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m18) ~ age) # residual variance varies between ages
boxplot(resid(m18) ~ sp) # and between species
# overall heterogeous

hist(resid(m18))
qqnorm(resid(m18))
qqline(resid(m18))
# good normality

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity needs improving
# -> model heterogeneity with generalised least squares

#### 5.11  Determine random components ####
m20 <- gls(N ~ sp * age,
           weights = varIdent(form = ~1|sp),
           method = "REML") # fixed effects model
m21 <- lme(N ~ sp * age,
           random = ~age|bag,
           weights = varIdent(form = ~1|sp),
           method = "REML") # mixed effects model
anova(m21, m20) # m20 is a better fit

m20 <- gls(N ~ sp * age,
           weights = varIdent(form = ~1|sp),
           method = "ML")

#### 5.12  Determine fixed components ####
drop1(m20, test = "Chisq") # retain interaction term

m20 <- gls(N ~ sp * age,
           weights = varIdent(form = ~1|sp))

#### 5.13  Test model fit ####
plot(m20, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m20, type = "normalized") ~ age)
boxplot(resid(m20, type = "normalized") ~ sp)
# homogeneity is better

hist(resid(m20))
qqnorm(resid(m20))
qqline(resid(m20))
# normality identical to m18

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m20 is chosen as the optimal model because of good homogeneity and normality

#### 5.14  Interpret model ####
Anova(m20, type = 3) # Type III Sums of Squares test
# Response: N
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 43.6712  3.884e-11 ***
# sp           2  1.4938    0.47384    
# age          1  0.5895    0.44263    
# sp:age       2  6.1595    0.04597 * 

summary(m20)
#                  Value  Std.Error   t-value p-value
# (Intercept)  1.5014053 0.22719577  6.608421   0.000
# sph         -0.0140179 0.25904450 -0.054114   0.957
# spd         -0.2167617 0.26108258 -0.830242   0.409
# age         -0.0070859 0.00922920 -0.767767   0.445
# sph:age      0.0178465 0.01052296  1.695960   0.094
# spd:age      0.0261450 0.01060575  2.465174   0.016

coef(m20)
# y = -0.007085874x + 1.501405262

sp <- factor(sp, levels = c("h", "d", "o"))
m20 <- gls(N ~ sp * age,
           weights = varIdent(form = ~1|sp))
Anova(m20, type = 3)
# Response: N
#             Df    Chisq Pr(>Chisq)    
# (Intercept)  1 142.8582    < 2e-16 ***
# sp           2   1.4938    0.47384    
# age          1   4.5311    0.03328 *  
# sp:age       2   6.1595    0.04597 *

summary(m20)
#                  Value  Std.Error   t-value p-value
# (Intercept)  1.4873873 0.12444329 11.952331  0.0000
# spd         -0.2027437 0.17897576 -1.132800  0.2609
# spo          0.0140179 0.25904450  0.054114  0.9570
# age          0.0107607 0.00505516  2.128646  0.0366
# spd:age      0.0082985 0.00727039  1.141411  0.2573
# spo:age     -0.0178465 0.01052296 -1.695960  0.0940

coef(m20)
# y = 0.010760652 + 1.487387337

sp <- factor(sp, levels = c("d", "h", "o"))
m20 <- gls(N ~ sp * age,
           weights = varIdent(form = ~1|sp))
Anova(m20, type = 3)
# Response: N
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 99.7395  < 2.2e-16 ***
# sp           2  1.4938  0.4738387    
# age          1 13.3040  0.0002648 ***
# sp:age       2  6.1595  0.0459706 * 

coef(m20)
# y = 0.019059162 + 1.284643598

sp <- factor(sp, levels = c("o", "h", "d"))

#### 5.13  Determine random components ####
m22 <- lm(CN ~ sp * age) # fixed effects model
m23 <- lmer(CN ~ sp * age + (age|bag), REML = F) # mixed effects model
anova(m23, m22) # models are not different
# continue with m22

#### 5.14  Determine fixed components ####
drop1(m22, test = "F") # remove interaction
m22 <- lm(CN ~ sp + age)

#### 5.15 Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m22, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m22) ~ age) # residual variance varies with age
boxplot(resid(m22) ~ sp) # and somewhat with species
# overall heterogenous

hist(resid(m22))
qqnorm(resid(m22))
qqline(resid(m22))
# residuals are right-skewed

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity and normality should be improved
# -> try fitting a gamma distribution to the data distribution

#### 5.16  Fit gamma distribution ####
require(fitdistrplus)
gamma <- fitdist(CN, "gamma") 
norm <- fitdist(CN, "norm")

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"),
        fitlty = 1) # judging visually, gamma fits slightly better
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# judging statistically, the gamma distribution fits slightly better

#### 5.17  Determine random components ####
m24 <- glm(CN ~ sp * age, family = Gamma(link = "log")) # fixed effects model
m25 <- glmer(CN ~ sp * age + (age|bag), family = Gamma(link = "log")) # mixed effects model
# too complex -> failed to converge
anova(m25, m24) # continue with m24

#### 5.18  Determine fixed components ####
drop1(m24, test = "Chisq") # remove interaction
m24 <- glm(CN ~ sp + age, family = Gamma(link = "log"))

#### 5.19  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m24, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m24) ~ age)
boxplot(resid(m24) ~ sp)
# homogeneity improved

hist(resid(m24))
qqnorm(resid(m24))
qqline(resid(m24))
# normality improved

#### 5.20  Interpret model ####
Anova(m24, type = 2) # Type II Sums of Squares test
# Response: CN
#     LR Chisq Df Pr(>Chisq)  
# sp    0.5357  2    0.76504  
# age   4.6207  1    0.03159 *

coef(m24)
# y = exp(-0.005495675x + 2.966531527)

sp <- factor(sp, levels = c("h", "d", "o"))
m24 <- glm(CN ~ sp + age, family = Gamma(link = "log"))
coef(m24)
# y = exp(-0.005495675x + 2.944914770)

sp <- factor(sp, levels = c("d", "h", "o"))
m24 <- glm(CN ~ sp + age, family = Gamma(link = "log"))
coef(m24)
# y = exp(-0.005495675x + 2.931236558)

sp <- factor(sp, levels = c("o", "h", "d"))

#### 6.   Decomposition~carbon data analysis ####
#### 6.1  Determine random components ####
m26 <- lm(loss ~ C * sp) # fixed effects model
m27 <- lmer(loss ~ C * sp + (C|bag), REML = F) # mixed effects model
# singular fit -> simplify mixed effects model
m27 <- lmer(loss ~ C * sp + (1|bag), REML = F) # mixed effects model
# still too complex
anova(m27, m26) # continue with m26

#### 6.2  Determine fixed components ####
drop1(m26, test = "F")
# remove interaction and species factor
m26 <- lm(loss ~ C)

#### 6.3  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m26, col = sp)
par(mfrow = c(1,1), mar = c(2,2,2,1))
plot(resid(m26) ~ C)
# homogeneity could be improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m26))
qqnorm(resid(m26))
qqline(resid(m26))
# good normality

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity can be improved with generalised least squares

#### 6.4  Determine random components ####
m28 <- gls(loss ~ C * sp,
           weights = varExp(form = ~C)) # fixed effects model
m29 <- lme(loss ~ C * sp,
           random = ~C|bag,
           weights = varExp(form = ~C)) # mixed effects model
anova(m29, m28) # continue with m28

#### 6.5  Determine fixed components ####
m28 <- gls(loss ~ C * sp,
           weights = varExp(form = ~C),
           method = "ML")
drop1(m28, test = "Chisq")
# remove interaction and species factor
m28 <- gls(loss ~ C,
           weights = varExp(form = ~C))

#### 6.6  Test model fit ####
plot(m28, col = sp)
plot(resid(m28, type = "normalized") ~ C)
# homogeneity is not significantly improved

#### 6.7  Interpret model ####
Anova(m26, type = 2) # Type II Sums of Squares test
# Response: loss
#            Sum Sq Df F value  Pr(>F)  
# C          11.515  1  5.7504 0.01884 *
# Residuals 158.198 79

coef(m26) # y = -0.08091949x + 3.32610403

#### 7.   Decomposition~CN data analysis ####
#### 7.1  Determine random components ####
m30 <- lm(loss ~ CN * sp) # fixed effects model
m31 <- lmer(loss ~ CN * sp + (CN|bag), REML = F) # mixed effects model
# singular fit -> simplify mixed effects model
m31 <- lmer(loss ~ CN * sp + (1|bag), REML = F) # mixed effects model
# still too complex
anova(m31, m30) # continue with m30

#### 7.2  Determine fixed components ####
drop1(m30, test = "F")
# remove interaction and species factor
m30 <- lm(loss ~ CN)

#### 7.3  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m30, col = sp)
par(mfrow = c(1,1), mar = c(2,2,2,1))
plot(resid(m30) ~ CN)
# good homogeneity

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m30))
qqnorm(resid(m30))
qqline(resid(m30))
# good normality

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 7.4 Interpret model ####
Anova(m30, type = 2) # Type II Sums of Squares test
# Response: loss
#            Sum Sq Df F value Pr(>F)
# CN          0.178  1  0.0831 0.7739
# Residuals 169.535 79

coef(m30) # y = 0.01541559x + 0.94607894

#### 8.   Decomposition~nitrogen data analysis ####
#### 8.1  Determine random components ####
m32 <- lm(loss ~ N * sp) # fixed effects model
m33 <- lmer(loss ~ N * sp + (N|bag), REML = F) # mixed effects model
# singular fit -> simplify mixed effects model
m33 <- lmer(loss ~ N * sp + (1|bag), REML = F) # mixed effects model
# still too complex
anova(m33, m32) # continue with m32

#### 8.2  Determine fixed components ####
drop1(m32, test = "F")
# remove interaction and species factor
m32 <- lm(loss ~ N)

#### 8.3  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m32, col = sp)
par(mfrow = c(1,1), mar = c(2,2,2,1))
plot(resid(m32) ~ N)
# good homogeneity

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m32))
qqnorm(resid(m32))
qqline(resid(m32))
# good normality

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 8.4  Interpret model ####
Anova(m32, type = 2) # Type II Sums of Squares test
# Response: loss
#            Sum Sq Df F value  Pr(>F)  
# N           7.577  1  3.6918 0.05829 .
# Residuals 162.136 79 

coef(m32) # y = -0.8947009x + 2.6374751

#### 9.   Grazing data analysis ####
#### 9.1  Determine random components ####
m34 <- lm(ex ~ sp2) # fixed effects model
m35 <- lmer(ex ~ sp2 + (1|bag2), REML = F) # mixed effects model
# singular fit -> model is too complex
anova(m35, m34) # models are not different
# continue with m34

#### 9.2  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m34, col = sp2)
par(mfrow = c(1,1), mar = c(2,2,2,1))
boxplot(resid(m34) ~ sp2)
# heterogenous: residual variance differs between species

par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m34))
qqnorm(resid(m34))
qqline(resid(m34))
# residuals deviate strongly at lower distribution edge

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity and normality need improving
# -> fit gamma distribution to data distribution

#### 9.3  Determine random components ####
m36 <- glm(ex+1 ~ sp2, family = Gamma(link = "log"))
m37 <- glmer(ex+1 ~ sp2 + (1|bag2), family = Gamma(link = "log"))
anova(m37, m36) # models are not different
# continue with m36

#### 9.4  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m36, col = sp2)
par(mfrow = c(1,1), mar = c(2,2,2,1))
boxplot(resid(m36) ~ sp2)
# homogeneity improved

par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m36))
qqnorm(resid(m36))
qqline(resid(m36))
# perfectly normal

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m32 is chosen as the optimal model based on its perfect normality
# and its good homogeneity

#### 9.5  Interpret model ####
Anova(m36, type = 2) # Type II Sums of Squares test
# Response: ex + 1
#     LR Chisq Df Pr(>Chisq)    
# sp2   104.94  2  < 2.2e-16 ***

summary(m36)
# L. och. vs. L. hyp., t = -7.392, p < 0.001 ***
# L. och. vs. L. dig., t = -9.867, p < 0.001 ***

sp2 <- factor(sp2, levels = c("h", "d", "o"))
m36 <- glm(ex+1 ~ sp2, family = Gamma(link = "log"))
summary(m36)
# L. hyp. vs. L. dig., t = -2.475, p = 0.02 *

sp2 <- factor(sp2, levels = c("o", "h", "d"))


#### 9.6  Determine random components ####
m38 <- lm(per ~ sp2) # fixed effects model
m39 <- lmer(per ~ sp2 + (1|bag2), REML = F) # mixed effects model
# singular fit -> model is too complex
anova(m39, m38) # models are not different
# continue with m38

#### 9.7  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m38, col = sp2)
par(mfrow = c(1,1), mar = c(2,2,2,1))
boxplot(resid(m38) ~ sp2)
# heterogenous: residual variance differs between species

par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m38))
qqnorm(resid(m38))
qqline(resid(m38))
# residuals are right-skewed

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# homogeneity and normality need improving
# -> fit gamma distribution to data distribution

#### 9.8  Determine random components ####
m40 <- glm(per+1 ~ sp2, family = Gamma(link = "log"))
m41 <- glmer(per+1 ~ sp2 + (1|bag2), family = Gamma(link = "log"))
anova(m41, m40) # models are not different
# continue with m40

#### 9.9  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m40, col = sp2)
par(mfrow = c(1,1), mar = c(2,2,2,1))
boxplot(resid(m40) ~ sp2)
# homogeneity improved

par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m40))
qqnorm(resid(m40))
qqline(resid(m40))
# normality improved

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m36 is chosen as the optimal model

#### 9.10  Interpret model ####
Anova(m40, type = 2) # Type II Sums of Squares test
# Response: per + 1
#     LR Chisq Df Pr(>Chisq)   
# sp2   12.903  2   0.001578 **

summary(m40)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.59234    0.09238   6.412 4.62e-08 ***
# sp2h        -0.45541    0.13065  -3.486  0.00102 ** 
# sp2d        -0.30739    0.13065  -2.353  0.02253 * 

sp2 <- factor(sp2, levels = c("h", "d", "o"))
m40 <- glm(per+1 ~ sp2, family = Gamma(link = "log"))
summary(m40)
#             Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  0.13694    0.09238   1.482  0.14442   
# sp2d         0.14802    0.13065   1.133  0.26255   
# sp2o         0.45541    0.13065   3.486  0.00102 **

sp2 <- factor(sp2, levels = c("o", "h", "d"))


#### 10.   Data visualisation ####
#### 10.1  Descriptive statistics ####
require(psych)
deco.stat <- describeBy(loss2, list(sp3, subs), mat = T)
deco.stat2 <- describeBy(loss2, sp3, mat = T)
ex.stat <- describeBy(ex, sp2, mat = T)
per.stat <- describeBy(per, sp2, mat = T)

phen.stat <- describeBy(phen, list(sp, age), mat = T)
phen.stat$group2 <- as.integer(phen.stat$group2) # turn age back into integer
phen.stat2 <- describeBy(phen, sp, mat = T)

loss.stat <- describeBy(loss, sp, mat = T)
CN.stat <- describeBy(CN, sp, mat = T)
CNl.stat <- cbind(CN.stat[,c(2,5,15)], loss.stat[,c(5, 15)])
rownames(CNl.stat) <- NULL
colnames(CNl.stat) <- c("species","CN","CNse","D","Dse")

pl.stat <- cbind(phen.stat2[,c(2,5,15)], loss.stat[,c(5, 15)])
rownames(pl.stat) <- NULL
colnames(pl.stat) <- c("species","P","Pse","D","Dse")

#### 10.2  Model predictions ####
new <- data.frame(age = rep(seq(13, 32, by = 0.5), 3),
                  sp = c(rep("d", 39), rep("h", 39), rep("o", 39)))

sp <- factor(sp, levels = c("d", "h", "o"))

# Phenolic content
m4 <- gls(phen ~ sp * age,
          weights = varIdent(form = ~1|sp))
new$phen.fit <- predict(m4, newdata = new)
modmat <-  model.matrix(formula(m4)[-2], new)
int <- diag(modmat %*% vcov(m4) %*% t(modmat))
new$phen.lo <- with(new, phen.fit - qnorm(0.975)*sqrt(int))
new$phen.hi <- with(new, phen.fit + qnorm(0.975)*sqrt(int))


# Intraspecific decomposition~phenol
m8 <- gls(loss ~ phen * sp, 
          weights = varIdent(form = ~1|sp))
Biochem$intrap.fit <- predict(m8)
modmat <-  model.matrix(formula(m8)[-2])
int <- diag(modmat %*% vcov(m8) %*% t(modmat))
Biochem$intrap.lo <- with(Biochem, intrap.fit - qnorm(0.975)*sqrt(int))
Biochem$intrap.hi <- with(Biochem, intrap.fit + qnorm(0.975)*sqrt(int))

# Interspecific decomposition~phenol
m12 <- gls(loss ~ phen, 
           weights = varPower(form = ~phen))
Biochem$interp.fit <- predict(m12)
modmat <-  model.matrix(formula(m12)[-2])
int <- diag(modmat %*% vcov(m12) %*% t(modmat))
Biochem$interp.lo <- with(Biochem, interp.fit - qnorm(0.975)*sqrt(int))
Biochem$interp.hi <- with(Biochem, interp.fit + qnorm(0.975)*sqrt(int))

# Decomposition~CN
m30 <- lm(loss ~ CN)
Biochem$CN.fit <- predict(m30)
Biochem$CN.lo <- with(Biochem, CN.fit - qnorm(0.975)*predict(m30, se.fit = T)$se.fit)
Biochem$CN.hi <- with(Biochem, CN.fit + qnorm(0.975)*predict(m30, se.fit = T)$se.fit)

#### 10.3  Customised theme ####
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
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title = element_blank(),
                 text = element_text(family = "Helvetica Neue"))

#### 10.4  Plots ####
# Decomposition
deco.stat$group1 <- factor(deco.stat$group1, levels = c("o", "d", "h"))
deco.stat2$group1 <- factor(deco.stat2$group1, levels = c("o", "d", "h"))

require(ggnewscale)
dp <- ggplot(deco.stat, aes(group1, mean)) +
        geom_col(aes(fill = group1, alpha = group2), width = 0.72,
                 position = position_dodge(width = 0.8)) +
        geom_errorbar(aes(colour = group2, ymin = mean - se, ymax = mean + se), width = 0.1,
                      position = position_dodge(width = 0.8)) +
        scale_colour_manual(values = c("#000000", "#000000"),
                            guide = "none") +
        new_scale_colour() +
        geom_text(aes(group1, 0.01, label = n, colour = group2), hjust = 0,
                  position = position_dodge(width = 0.8), size = 4.2) +
        scale_colour_manual(values = c("#ffffff", "#000000"),
                            guide = "none") +
        new_scale_colour() +
        geom_rug(data = deco.stat2, aes(0.8, mean, colour = group1), sides = "r",
                 length = unit(.25, "cm")) +
        scale_colour_manual(values = c("#f1c700","#333b08","#627d0e"),
                            guide = "none") +
        geom_text(aes(group1, rep(c(1.8148850+0.03,0.7855972+0.03,1.0138941+0.03),2),
                      label = rep(c("b","a","a"),2)),
                  hjust = 0, vjust = 0.35, position = position_dodge(width = 0.8),
                  size = 4.2, check_overlap = T) +
        scale_fill_manual(values = c("#f1c700","#333b08","#627d0e"),
                          guide = "none") +
        scale_alpha_manual(values = c(1, 0.5)) +
        ylab(expression("Decomposition (% d"^-1*")")) +
        scale_y_continuous(expand = c(0,0), position = "right") +
        coord_flip(ylim = c(0, 2.5)) +
        theme(legend.position = c(0.85, 0.923),
              axis.title.y = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank()) +
        mytheme

dp

require(ggridges)
Deco$species <- factor(Deco$species, levels = c("d", "h", "o"))
Deco$experiment <- c(rep("a2019b",27), rep("b2016b",16), rep("a2019c",27),
                     rep("b2016c",16), rep("a2019a",27), rep("b2016a",16))
deco.stat$group1 <- factor(deco.stat$group1, levels = c("d", "h", "o"))
deco.stat2$group1 <- factor(deco.stat2$group1, levels = c("d", "h", "o"))
deco.stat$experiment <- c("a2019a","a2019c","a2019b","b2016a","b2016c","b2016b")

dp.alt <- ggplot() +
        geom_vline(aes(xintercept = 0)) +
        geom_density_ridges(data = Deco,
                            aes(x = perc.loss, y = experiment,
                                fill = species, colour = species),
                            scale = 0.7, bandwidth = 0.4, calc_ecdf = TRUE, alpha = 0.3, point_size = 3,
                            point_shape = 16, point_alpha = 0.5, quantiles = c(0.025, 0.5, 0.975),
                            quantile_lines = TRUE, jittered_points = TRUE, position = "raincloud") +
        geom_pointrange(data = deco.stat, aes(x = mean, xmin = mean - qnorm(0.975)*se,
                                              xmax = mean + qnorm(0.975)*se, y = experiment)) +
        geom_rug(data = deco.stat2, aes(mean, 1, colour = group1), outside = T,
                 sides = "t", length = unit(.25, "cm")) +
        scale_fill_manual(values = c("#333b08","#627d0e","#f1c700"),
                          labels = c(expression(italic("L. digitata")),
                                     expression(italic("L. hyperborea")),
                                     expression(italic("L. ochroleuca"))),
                          guide = guide_legend()) +
        scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")),
                                       expression(italic("L. hyperborea")),
                                       expression(italic("L. ochroleuca"))),
                            guide = guide_legend()) +
        # annotate("text", x = -2.92, y = c(1, 4), label = c("2019 experiment", "2016 experiment"),
        #          hjust = 0, size = 4.2, angle = 90, family = "Helvetica Neue") +
        annotate("segment", y = c(1, 4), yend = c(3, 6), x = -2.8, xend = -2.8) +
        xlab(expression("Decomposition (% d"^-1*")")) +
        scale_y_discrete(labels = c("", "2019 experiment", "", "", "2016 experiment", "")) +
        scale_x_continuous(breaks = seq(-3, 6, by = 3), expand = c(0, 0), position = "top") +
        coord_cartesian(xlim = c(-3, 6), clip = "off") +
        mytheme +
        theme(legend.position = "none",
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_text(angle = 90, hjust = 0.5),
              strip.text = element_blank(),
              strip.background = element_blank(),)

dp.alt

# Excavation
ex.stat$group1 <- factor(ex.stat$group1, levels = c("o", "d", "h"))
exp <- ggplot(ex.stat, aes(group1, mean)) +
        geom_col(aes(fill = group1), width = 0.9) +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
        geom_text(aes(group1, mean + se + 1, label = c("c","a","b")),
                  hjust = 0, vjust = 0.35, size = 4.2) +
        scale_fill_manual(values = c("#f1c700","#333b08","#627d0e"),
                          guide = "none") +
        scale_alpha_manual(values = c(1, 0.5)) +
        ylab("Excavated tissue (%)") +
        scale_y_continuous(breaks = seq(0, 90, by = 30),
                           expand = c(0,0), position = "right") +
        coord_flip(ylim = c(0, 90)) +
        theme(axis.title.y = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank()) +
        mytheme

exp

Grazing$species <- factor(Grazing$species, levels = c("o", "d", "h"))
ex.stat$group1 <- factor(ex.stat$group1, levels = c("o", "d", "h"))
exp.alt <- ggplot() +
        geom_density_ridges(data = Grazing,
                            aes(x = excavation, y = species,
                                fill = species, colour = species),
                            scale = 0.7, bandwidth = 2, calc_ecdf = TRUE, alpha = 0.3, point_size = 3,
                            point_shape = 16, point_alpha = 0.5, quantiles = c(0.025, 0.5, 0.975),
                            quantile_lines = TRUE, jittered_points = TRUE, position = "raincloud") +
        geom_pointrange(data = ex.stat, aes(x = mean, xmin = mean - qnorm(0.975)*se,
                                            xmax = mean + qnorm(0.975)*se, y = group1)) +
        scale_fill_manual(values = c("#f1c700", "#333b08","#627d0e"),
                          labels = c(expression(italic("L. ochroleuca")),
                                     expression(italic("L. digitata")),
                                     expression(italic("L. hyperborea"))),
                          guide = guide_legend()) +
        scale_colour_manual(values = c("#f1c700", "#333b08","#627d0e"),
                            labels = c(expression(italic("L. ochroleuca")),
                                       expression(italic("L. digitata")),
                                       expression(italic("L. hyperborea"))),
                            guide = guide_legend()) +
        xlab("Excavated tissue (%)") +
        scale_x_continuous(breaks = seq(0, 100, by = 20), expand = c(0, 0), position = "top") +
        coord_cartesian(xlim = c(0, 100), clip = "off") +
        mytheme +
        theme(legend.position = "none",
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              strip.text = element_blank(),
              strip.background = element_blank(),)

exp.alt

# Perforation
per.stat$group1 <- factor(per.stat$group1, levels = c("o", "d", "h"))
perp <- ggplot(per.stat, aes(group1, mean)) +
        geom_col(aes(fill = group1), width = 0.9) +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
        geom_text(aes(group1, mean + se + 0.02, label = c("b","a","a")),
                  hjust = 0, vjust = 0.35, size = 4.2) +
        scale_fill_manual(values = c("#f1c700","#333b08","#627d0e"),
                          guide = "none") +
        scale_alpha_manual(values = c(1, 0.5)) +
        ylab("Perforated tissue (%)") +
        scale_y_continuous(breaks = seq(0, 1.2, by = 0.4),
                           expand = c(0,0), position = "right") +
        coord_flip(ylim = c(0, 1.2)) +
        theme(axis.title.y = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank()) +
        mytheme

perp

per.stat$group1 <- factor(per.stat$group1, levels = c("o", "d", "h"))
perp.alt <- ggplot() +
        geom_density_ridges(data = Grazing,
                            aes(x = perforation, y = species,
                                fill = species, colour = species),
                            scale = 0.7, bandwidth = 0.06, calc_ecdf = TRUE, alpha = 0.3, point_size = 3,
                            point_shape = 16, point_alpha = 0.5, quantiles = c(0.025, 0.5, 0.975),
                            quantile_lines = TRUE, jittered_points = TRUE, position = "raincloud") +
        geom_pointrange(data = per.stat, aes(x = mean, xmin = mean - qnorm(0.975)*se,
                                            xmax = mean + qnorm(0.975)*se, y = group1)) +
        scale_fill_manual(values = c("#f1c700", "#333b08","#627d0e"),
                          labels = c(expression(italic("L. ochroleuca")),
                                     expression(italic("L. digitata")),
                                     expression(italic("L. hyperborea"))),
                          guide = guide_legend()) +
        scale_colour_manual(values = c("#f1c700", "#333b08","#627d0e"),
                            labels = c(expression(italic("L. ochroleuca")),
                                       expression(italic("L. digitata")),
                                       expression(italic("L. hyperborea"))),
                            guide = guide_legend()) +
        xlab("Perforated tissue (%)") +
        scale_x_continuous(breaks = seq(0, 4, by = 1), expand = c(0, 0), position = "top") +
        coord_cartesian(xlim = c(0, 4), clip = "off") +
        mytheme +
        theme(legend.position = "none",
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              strip.text = element_blank(),
              strip.background = element_blank(),)

perp.alt

# Phenolic content across detrital ages
pp <- ggplot()+
        geom_line(data = new, aes(age, phen.fit, colour = sp, lty = sp), size = 0.5) +
        geom_ribbon(data = new, aes(age, ymin = phen.lo, ymax = phen.hi, fill = sp),
                    alpha = .5) +
        geom_pointrange(data = phen.stat, aes(group2, mean, ymin = mean - se,
                                              ymax = mean + se, colour = group1),
                        size = 0.5, position = position_dodge(width = 2.16)) +
        scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")*"        y = 0.001x + 0.09"),
                                       expression(italic("L. hyperborea")*"  y = 0.026x + 0.45"),
                                       expression(italic("L. ochroleuca")*"  y = 0.003x + 0.09")),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#333b08","#627d0e","#f1c700"),
                          labels = c(expression(italic("L. digitata")*"        y = 0.001x + 0.09"),
                                     expression(italic("L. hyperborea")*"  y = 0.026x + 0.45"),
                                     expression(italic("L. ochroleuca")*"  y = 0.003x + 0.09")),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(5, 1, 5),
                              guide = "none") +
        ylab("Phenolic content (%)") +
        xlab("Detrital age (d)") +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0)) +
        coord_cartesian(ylim = c(0, 2), xlim = c(10, 35)) +
        theme(legend.position = c(0.63, 0.905)) +
        mytheme

pp

Biochem$species <- factor(Biochem$species, levels = c("d", "h", "o"))
# Biochem$spage <- with(Biochem, paste(species, age, sep = ""))
pp.alt <- ggplot()+
        geom_line(data = new, aes(age, phen.fit, colour = sp, lty = sp), size = 0.5) +
        geom_ribbon(data = new, aes(age, ymin = phen.lo, ymax = phen.hi, fill = sp),
                    alpha = .5) +
        # geom_violin(data = Biochem, aes(age, phenols, group = spage, colour = species),
        #             fill = NA, width = 10, position = position_dodge(width = 2.16)) +
        geom_point(data = Biochem, aes(age, phenols, colour = species), 
                   alpha = 0.4, shape = 16, size = 3,
                   position = position_dodge(width = 2.16)) +
        geom_pointrange(data = phen.stat, aes(group2, mean, ymin = mean - se*qnorm(0.975),
                                              ymax = mean + se*qnorm(0.975), colour = group1),
                        size = 0.5, position = position_dodge(width = 2.16)) +
        scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")*"        y = 0.001x + 0.09"),
                                       expression(italic("L. hyperborea")*"  y = 0.026x + 0.45"),
                                       expression(italic("L. ochroleuca")*"  y = 0.003x + 0.09")),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#333b08","#627d0e","#f1c700"),
                          labels = c(expression(italic("L. digitata")*"        y = 0.001x + 0.09"),
                                     expression(italic("L. hyperborea")*"  y = 0.026x + 0.45"),
                                     expression(italic("L. ochroleuca")*"  y = 0.003x + 0.09")),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(5, 1, 5),
                              guide = "none") +
        ylab("Phenolic content (%)") +
        xlab("Detrital age (d)") +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0)) +
        coord_cartesian(ylim = c(0, 2), xlim = c(10, 35)) +
        theme(legend.position = c(0.63, 0.905)) +
        mytheme + theme(plot.margin = unit(c(.2, .3, .2, .2),"cm"))

pp.alt

pp2 <- ggplot() +
        geom_hline(yintercept = 0) +
        geom_line(data = Biochem, aes(phenols, intrap.fit, colour = sp, lty = sp), size = 0.5) +
        geom_ribbon(data = Biochem, aes(phenols, ymin = intrap.lo, ymax = intrap.hi, fill = sp), 
                    alpha = .5) +
        geom_point(data = Biochem, aes(phenols, loss, colour = sp),
                   size = 3, shape = 16, alpha = 0.4) +
        scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")*"        y = –20.99x + 3.32"),
                                       expression(italic("L. hyperborea")*"  y = –2.44x + 3.4"),
                                       expression(italic("L. ochroleuca")*"  y = –9.23x + 3.29")),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#333b08","#627d0e","#f1c700"),
                          labels = c(expression(italic("L. digitata")*"        y = –20.99x + 3.32"),
                                     expression(italic("L. hyperborea")*"  y = –2.44x + 3.4"),
                                     expression(italic("L. ochroleuca")*"  y = –9.23x + 3.29")),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(5, 1, 1),
                              guide = "none") +
        ylab(expression("Decomposition (% d"^-1*")")) +
        xlab("Phenolic content (%)") +
        scale_y_continuous(expand = c(0,0), breaks = seq(-3, 6, by = 3)) +
        scale_x_continuous(expand = c(0,0)) +
        coord_cartesian(xlim = c(0, 2), ylim = c(-3, 6), clip = "off") +
        theme(legend.position = c(0.67, 0.905)) +
        mytheme

pp2 # dimensions 4 x 5 in

CNp <- ggplot() +
        geom_hline(yintercept = 0) +
        geom_point(data = Biochem, aes(CN, perc.loss, colour = species), 
                   alpha = 0.4, shape = 16, size = 3) +
        # stat_ellipse(data = Biochem, aes(CN, perc.loss, colour = species), level = 0.075, type = "norm") +
        geom_line(data = Biochem, aes(CN, CN.fit), lty = 5) +
        geom_ribbon(data = Biochem, aes(CN, ymin = CN.lo, ymax = CN.hi), 
                    fill = NA, colour = "#000000") +
        geom_pointrange(data = CNl.stat, 
                        aes(CN, D, ymin = D - Dse*qnorm(0.975), ymax = D + Dse*qnorm(0.975), colour = species)) +
        geom_errorbarh(data = CNl.stat, 
                       aes(y = D, xmin = CN - CNse*qnorm(0.975), xmax = CN + CNse*qnorm(0.975), colour = species), 
                       height = 0) +
        geom_text(data = Biochem, aes(10.6, 5.8), label = "y = 0.02x + 0.95", 
                  family = "Helvetica Neue", check_overlap = T, hjust = 0, size = 4.2) +
        scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")),
                                       expression(italic("L. hyperborea")),
                                       expression(italic("L. ochroleuca"))),
        guide = guide_legend()) +
        ylab(expression("Decomposition (% d"^-1*")")) +
        xlab("Carbon-nitrogen ratio") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0), breaks = seq(10, 35, by = 5)) +
        coord_cartesian(ylim = c(-2, 6), xlim = c(10, 35)) +
        theme(legend.position = "none") +
        mytheme
CNp

pp3 <- ggplot() +
        geom_hline(yintercept = 0) +
        geom_point(data = Biochem, aes(phenols, perc.loss, colour = species), 
                   alpha = 0.4, shape = 16, size = 3) +
        geom_line(data = Biochem, aes(phenols, interp.fit)) +
        geom_ribbon(data = Biochem, aes(phenols, ymin = interp.lo, ymax = interp.hi), 
                    fill = NA, colour = "#000000") +
        geom_pointrange(data = pl.stat, 
                        aes(P, D, ymin = D - Dse*qnorm(0.975), ymax = D + Dse*qnorm(0.975), colour = species)) +
        geom_errorbarh(data = pl.stat, 
                       aes(y = D, xmin = P - Pse*qnorm(0.975), xmax = P + Pse*qnorm(0.975), colour = species), 
                       height = 0) +
        geom_text(data = Biochem, aes(0.05, 5.8), label = "y = –1.05x + 1.66", 
                  family = "Helvetica Neue", check_overlap = T, hjust = 0, size = 4.2) +
        scale_colour_manual(values = c("#333b08","#627d0e","#f1c700"),
                            labels = c(expression(italic("L. digitata")),
                                       expression(italic("L. hyperborea")),
                                       expression(italic("L. ochroleuca"))),
                            guide = guide_legend()) +
        ylab(expression("Decomposition (% d"^-1*")")) +
        xlab("Phenolic content (%)") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1.8, by = 0.6)) +
        coord_cartesian(ylim = c(-2, 6), xlim = c(0, 1.8)) +
        theme(legend.position = "none") +
        mytheme
pp3

#### 10.5  Combined plots ####
require(cowplot)
legend <- get_legend(dp.alt + theme(legend.position = "top",
                                    legend.justification = "left"))
plots <- align_plots(dp.alt, exp.alt + theme(plot.margin = unit(c(.2, .5, .2, .5),"cm")), 
                     perp.alt + theme(plot.margin = unit(c(.2, .5, .2, .5),"cm")), 
                     align = "h", axis = "t")
plots2 <- align_plots(dp.alt, CNp, align = "v", axis = "l")

right_col <- plot_grid(plots[[2]], plots[[3]], 
                       labels = c("b", "c"), label_size = 15, label_fontfamily = "Helvetica Neue",
                       ncol = 1, hjust = 0.5)
left_col <- plot_grid(plots2[[1]], legend, labels = c("a", ""), 
                      label_size = 15, ncol = 1, rel_heights = c(1, 0.1),
                      label_fontfamily = "Helvetica Neue")
top_row <- plot_grid(left_col, right_col, nrow = 1, rel_widths = c(1, 0.92))
  
bottom_row <- plot_grid(plots2[[2]], pp3 + theme(axis.text.y = element_blank(),
                                                 axis.title.y = element_blank(),
                                                 plot.margin = unit(c(.2, .5, .2, .5),"cm")),
                        labels = c("d", "e"), label_size = 15, label_fontfamily = "Helvetica Neue",
                        nrow = 1, rel_widths = c(1, 0.92), hjust = c(-0.5, 0.5))

Fig.3 <- plot_grid(top_row, bottom_row, ncol = 1)
Fig.3 # dimensions 8 x 8.5 in

#### 11.   Clean up ####
detach(package:fitdistrplus)
detach(package:car)
detach(package:lme4)
detach(package:nlme)
detach(package:psych)
detach(package:ggnewscale)
detach(package:cowplot)
detach(package:ggplot2)
rm(list = ls())
graphics.off()
cat("\014")
