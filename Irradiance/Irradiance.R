##################################################################
##### Project: Sequestration of carbon from Laminaria kelp   #####
#####          forests may be diminished in a warmer climate #####
##### Script purpose: Analysis of PAR data from station L4   #####
##### Author: Luka Seamus Wright                             #####
##################################################################

#### 1.   Data preparation ####
#### 1.1  Load data ####
# Data source: https://www.westernchannelobservatory.org.uk/l4_ctdf/index.php
L4 <- read.csv("~/Desktop/PATH/L4.csv")

#### 1.2  Exploratory data visualisation ####
with(L4, plot(PAR ~ Depth)) # clearly there is one unlikely outlier above 2000
which(L4$PAR > 2000) # the outlier is in row 17940
L4 <- L4[-c(17940),] # delete row with outlier
with(L4, plot(PAR ~ Depth)) # looks better

#### 1.3  Rename variables ####
PAR <- L4$PAR # irradiance (μmol photons m-2 s-1)
depth <- L4$Depth # depth (m)
year <- L4$Year
season <- L4$Season

#### 2.   Data analysis ####
#### 2.1  Build annual model ####
require(lme4)
m1 <- lm(log(PAR) ~ depth)
m2 <- lmer(log(PAR) ~ depth + (1|year), REML = F)
anova(m2, m1) # continue with m2
m2 <- lmer(log(PAR) ~ depth + (1|year), REML = T)

#### 2.2  Test model fit ####
# assumption of homogeneity
plot(m2) # some outlying trends where irradiance bottoms out
# but overall quite homogenous

# assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m2))
qqnorm(resid(m2))
qqline(resid(m2)) # somewhat left-skewed, but acceptable
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 2.3  Interpret model ####
require(car)
Anova(m2, type = 2) # Type II Sums of Squares test
# Response: log(PAR)
#        Chisq Df Pr(>Chisq)    
# depth 190169  1  < 2.2e-16 ***

summary(m2)
# function: y = -0.1520991x + 5.0476945

#### 2.4  Build seasonal model ####
m3 <- lm(log(PAR) ~ depth * season)
m4 <- lmer(log(PAR) ~ depth * season + (1|year), REML = F)
anova(m4, m3) # continue with m4

m5 <- lmer(log(PAR) ~ depth + season + (1|year), REML = F)
anova(m5, m4) # continue with m4

m4 <- lmer(log(PAR) ~ depth * season + (1|year), REML = T)

#### 2.5  Test model fit ####
# assumption of homogeneity
plot(m4) # some outlying trends where irradiance bottoms out
# but overall quite homogenous

# assumption of normality
par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m4))
qqnorm(resid(m4))
qqline(resid(m4)) # somewhat left-skewed, but acceptable
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 2.6  Interpret model ####
Anova(m4, type = 3) # Type III Sums of Squares test
# Response: log(PAR)
#                 Chisq Df Pr(>Chisq)    
# (Intercept)    821.43  1  < 2.2e-16 ***
# depth        61751.65  1  < 2.2e-16 ***
# season        9920.59  3  < 2.2e-16 ***
# depth:season   676.14  3  < 2.2e-16 ***

summary(m4)
# Autumn function: y = -0.1447491x + 4.6215877

season <- factor(season, levels = c("Spring", "Summer", "Autumn", "Winter"))
m4 <- lmer(log(PAR) ~ depth * season + (1|year), REML = T)
Anova(m4, type = 3)
# Response: log(PAR)
#                 Chisq Df Pr(>Chisq)    
# (Intercept)   1255.51  1  < 2.2e-16 ***
# depth        93229.54  1  < 2.2e-16 ***
# season        9920.59  3  < 2.2e-16 ***
# depth:season   676.14  3  < 2.2e-16 ***

summary(m4)
# Spring function: y = -0.1636840x + 5.7091730

season <- factor(season, levels = c("Summer", "Autumn", "Winter", "Spring"))
m4 <- lmer(log(PAR) ~ depth * season + (1|year), REML = T)
Anova(m4, type = 3)
# Response: log(PAR)
#                 Chisq Df Pr(>Chisq)    
# (Intercept)   1304.75  1  < 2.2e-16 ***
# depth        74103.59  1  < 2.2e-16 ***
# season        9920.59  3  < 2.2e-16 ***
# depth:season   676.14  3  < 2.2e-16 ***

summary(m4)
# Summer function: y = -0.1496887x + 5.8214088

season <- factor(season, levels = c("Winter", "Autumn", "Summer", "Spring"))
m4 <- lmer(log(PAR) ~ depth * season + (1|year), REML = T)
Anova(m4, type = 3)
# Response: log(PAR)
#                 Chisq Df Pr(>Chisq)    
# (Intercept)    475.19  1  < 2.2e-16 ***
# depth        50870.17  1  < 2.2e-16 ***
# season        9920.59  3  < 2.2e-16 ***
# depth:season   676.14  3  < 2.2e-16 ***

summary(m4)
# Winter function: y = -0.1484400x + 3.5215956

#### 4.   Data visualisation ####
#### 4.1  Model predictions ####
L4$annu.fit <- predict(m2, re.form = NA)
L4$seas.fit <- predict(m4, re.form = NA)

#### 4.2  Customised theme ####
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

#### 4.3  Plot ####
ap <- ggplot(L4, aes(Depth, PAR)) +
        geom_point(colour = "#dbdddf") +
        geom_line(aes(y = exp(annu.fit)), size = 0.5) +
        ylab(expression("Irradiance ("*mu*"mol photons m"^-2*" s"^-1*")")) +
        xlab("Depth (m)") +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(0, 1500), xlim = c(0, 60)) +
        mytheme

ap # dimsenisons: 4 x 7 in

L4$Season <- factor(L4$Season, levels = c("Spring", "Summer", "Autumn", "Winter"))
sp <- ggplot(L4, aes(Depth, PAR, colour = Season)) +
        geom_point(colour = "#dbdddf") +
        geom_line(aes(y = exp(seas.fit)), size = 0.5) +
        scale_colour_manual(values = c("#c7b300", "#50590d", "#6b4d8d", "#003875"),
                            labels = c("Spring", "Summer", "Autumn", "Winter"),
                            guide = guide_legend()) +
        ylab(expression("Irradiance ("*mu*"mol photons m"^-2*" s"^-1*")")) +
        xlab("Depth (m)") +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim = c(0, 1500), xlim = c(0, 60)) +
        theme(legend.position = c(0.92, 0.87)) +
        mytheme

sp # dimsenisons: 4 x 7 in

#### 5.   Clean up ####
detach(package:car)
detach(package:lme4)
detach(package:ggplot2)
rm(list = ls())
graphics.off()
cat("\014")
