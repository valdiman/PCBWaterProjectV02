## Water PCB concentrations data analysis per site
# Chesapeake Bay & Delaware Canal

# Install packages
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("robustbase")
install.packages("dplyr")
install.packages("tibble")
install.packages("Matrix")
install.packages("lme4")
install.packages("MuMIn")
install.packages("lmerTest")
install.packages("Matrix")
install.packages("zoo")
install.packages("dataRetrieval")

# Load libraries
library(ggplot2)
library(scales) # function trans_breaks
library(stringr) # str_detect
library(robustbase) # function colMedians
library(dplyr) # performs %>%
library(tibble) # adds a column
library(lme4) # performs lme
library(MuMIn) # gets Rs from lme
library(lmerTest) # gets the p-value from lme
library(zoo) # yields seasons
library(dataRetrieval) # read data from USGS

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("WaterDataCongenerAroclor08052022.csv")

# Select Chesapeake Bay & Delaware Canal data ---------------------------------------------------
che.0 <- wdc[str_detect(wdc$SiteName, 'Chesapeake Bay'),]

# Data preparation --------------------------------------------------------
# Calculate total PCB
tpcb.che <- rowSums(che.0[, c(14:117)], na.rm = T)
# Change date format
che.0$SampleDate <- as.Date(che.0$SampleDate, format = "%m/%d/%y")
# Calculate sampling time
time.day <- data.frame(as.Date(che.0$SampleDate) - min(as.Date(che.0$SampleDate)))
# Create individual code for each site sampled
site.numb <- che.0$LocationID %>% as.factor() %>% as.numeric
# Include season
yq.s <- as.yearqtr(as.yearmon(che.0$SampleDate, "%m/%d/%Y") + 1/12)
season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                   labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
# Create data frame
che.tpcb <- cbind(factor(che.0$LocationID), che.0$SampleDate,
                  che.0$Latitude, che.0$Longitude, as.matrix(tpcb.che),
                  data.frame(time.day), site.numb, season.s)
# Add column names
colnames(che.tpcb) <- c("LocationID", "date", "Latitude", "Longitude",
                        "tPCB", "time", "site.code", "season")

# Get coordinates per site to plot in Google Earth
che.location <- che.tpcb[c('LocationID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
che.location <- aggregate(tPCB ~ LocationID + Latitude + Longitude,
                          data = che.location, mean)

# (2) Calculate total log PCB
# Remove metadata
che.log <- subset(che.0, select = -c(SampleID:AroclorCongener))
# Remove Aroclor data
che.log <- subset(che.log, select = -c(A1016:A1260))
# Log 10 individual PCBs 
che.log <- log10(che.log)
# Replace -inf to NA
che.log <- do.call(data.frame,
                   lapply(che.log,
                          function(x) replace(x, is.infinite(x), NA)))
# Sum individual log 10 PCBs
che.log.tpcb <- rowSums(che.log, na.rm = T)
# Generate data.frame for analysis and plots
che.log.tpcb <- cbind(factor(che.0$LocationID), che.0$SampleDate,
                      as.matrix(che.log.tpcb), data.frame(time.day),
                      site.numb, season.s)
colnames(che.log.tpcb) <- c("LocationID", "date", "logtPCB", "time",
                            "site.code", "season")

# General plots -------------------------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(che.tpcb$tPCB)
hist(log10(che.tpcb$tPCB))
# (1.2) log.tPCB
hist(che.log.tpcb$logtPCB)
hist(log10(che.log.tpcb$logtPCB))

# (2) Time trend plots
# (2.1) tPCB
ggplot(che.tpcb, aes(y = tPCB,
                     x = format(date,'%Y'))) +
  geom_point() +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black"))

# (2.2) log.tPCB
ggplot(che.log.tpcb, aes(y = logtPCB,
                         x = format(date,'%Y'))) +
  geom_point() +
  xlab("") +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black"))

# (3) Seasonality
# (3.1) tPCB
ggplot(che.tpcb, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Concentration " *Sigma*"PCB 2012 - 2018 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (3.2) log.tPCB
ggplot(che.log.tpcb, aes(x = season, y = logtPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Concentration " *Sigma*"PCB 2012 - 2018 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (4) Sites
# From ~north to ~south
sites <- c("SusquehannaRiverMDPABorder", "LittleElkCreekTelegraphRd",
           "SusquehannaRiverConowingoBridge", "BigElkCreekHWY40",
           "NorthEastRiverNEIslesDr", "ScottRunBiddlePoint",
           "NorthEastRiverRoachsShore", "PennysShoalParkIsland",
           "ChesapeakeAndDelawareCanalGoosePt",
           "CranberryRunCranberryRd", "ChurchCreekPulaskiHwy",
           "GreatBohemiaCreekOldTelegraphRd",
           "GreatBohemiaCreekWoodstockFarmLn", "BohemiaRiverFreeSchoolPt",
           "TurkeyPoint", "LakeRolandSunsetRock", "HenIslandCreekWilsonPt",
           "BackRiverHWY40", "BayAtBowleyPt", "BullneckCreekMerrittPoint",
           "HawkCovePleasureIsland", "BayAtArcadia",
           "MagothyRiverGibsonIsland", "CorsicaRiverYellowBankStream",
           "SevernRiverMouth", "SouthRiverMouth", "RhodeRiverCadleCreekConflux",
           "RhodeWestConflux", "PatuxentRiverEagleHarbor", "StLeonardCreekParranRd",
           "PatuxentRiverBarrettIsland", "PatuxentRiverMouth", "ChesapeakeBayMouth")

ggplot(che.tpcb, aes(x = factor(site, levels = sites), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Water Concentration " *Sigma*"PCB 2012 - 2018 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0, color = "black") +
  theme(legend.position = "none")

# (4.1) tPCB
ggplot(che.tpcb, aes(x = factor(site), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Water Concentration " *Sigma*"PCB 2012 - 2018 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# (4.2) log.tPCB
ggplot(che.log.tpcb, aes(x = factor(site), y = logtPCB)) + 
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Water Concentration " *Sigma*"PCB 2012 - 2018 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# Regressions -------------------------------------------------------------
# All data
# (1) Perform linear regression (lr)
# (1.1) tPCB vs. time
lr.che.tpcb.t <- lm(log10(tPCB) ~ time, data = che.tpcb)
# See results
summary(lr.che.tpcb.t)
# Look at residuals
res <- resid(lr.che.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.che.log.tpcb.t <- lm(logtPCB ~ time, data = che.log.tpcb)
# See results
summary(lr.che.log.tpcb.t)
# Look at residuals
res <- resid(lr.che.log.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB vs. season
lr.che.tpcb.s <- lm(log10(tPCB) ~ season, data = che.tpcb)
# See results
summary(lr.che.tpcb.s)
# Look at residuals
res <- resid(lr.che.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.che.log.tpcb.s <- lm(logtPCB ~ season, data = che.log.tpcb)
# See results
summary(lr.che.log.tpcb.s)
# Look at residuals
res <- resid(lr.che.log.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) MLR
# (2.1) tPCB vs. time + season (che.tpcb)
mlr.che.tpcb <- lm(log10(tPCB) ~ time + season, data = che.tpcb)
# See results
summary(mlr.che.tpcb)
# Look at residuals
res <- resid(mlr.che.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2.2) log.tPCB vs. time + season + flow + temp (che.log.tpcb)
mlr.che.log.tpcb <- lm(logtPCB ~ time + season,
                       data = che.log.tpcb)
# See results
summary(mlr.che.log.tpcb)
# Look at residuals
res <- resid(mlr.che.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (3) Perform Linear Mixed-Effects Model (LMEM)
# (3.1) tPCB vs. time + season + site (che.tpcb)
tpcb <- che.tpcb$tPCB
time <- che.tpcb$time
site <- che.tpcb$site.code
season <- che.tpcb$season

lmem.che.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.che.tpcb)
# Look at residuals
res.che.tpcb <- resid(lmem.che.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.che.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.che.tpcb)
# Shapiro test
shapiro.test(res.che.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.che.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.che.tpcb))[1, 'R2c']

# Extract coefficient values
time.coeff <- summary(lmem.che.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lmem.che.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -log(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# (3.2) log.tPCB vs. time + season + site (che.log.tpcb)
log.tpcb <- che.log.tpcb$logtPCB
time <- che.log.tpcb$time
site <- che.log.tpcb$site.code
season <- che.log.tpcb$season

lmem.che.log.tpcb <- lmer(log.tpcb ~ 1 + time + season + (1|site),
                          REML = FALSE,
                          control = lmerControl(check.nobs.vs.nlev = "ignore",
                                                check.nobs.vs.rankZ = "ignore",
                                                check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.che.log.tpcb)
# Look at residuals
res.che.log.tpcb <- resid(lmem.che.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.che.log.tpcb, main = expression(paste("Normal Q-Q Plot", " (", Sigma,
                                                 "log"[10]*"PCB)")))
# Add a straight diagonal line to the plot
qqline(res.che.log.tpcb)
# Shapiro test
shapiro.test(res.che.log.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res.che.log.tpcb, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.che.log.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.che.log.tpcb))[1, 'R2c']

# Predictions -------------------------------------------------------------
# Modeling plots
# (1) Get predicted values tpcb
fit.values.che.tpcb <- as.data.frame(fitted(lmem.che.tpcb))
# Add column name
colnames(fit.values.che.tpcb) <- c("predicted")
# Add predicted values to data.frame
che.tpcb$predicted <- 10^(fit.values.che.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(che.tpcb, aes(x = tPCB, y = predicted)) +
  geom_point() +
  scale_x_log10(limits = c(1, 1e6)) +
  scale_y_log10(limits = c(1, 1e6)) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 0.8) + # 1:1 line
  geom_abline(intercept = 0.5, slope = 1, col = "blue", size = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = -0.5, slope = 1, col = "blue", size = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 30, y = 1000000,
           label = 'Chesapeake Bay', colour = 'black', size = 4,
           fontface = 2)

# Plot residuals vs. predictions
plot(che.tpcb$predicted, res.che.tpcb)
abline(0, 0)

# (2) Get predicted values log.tpcb
fit.values.che.log.tpcb <- as.data.frame(fitted(lmem.che.log.tpcb))
# Add column name
colnames(fit.values.che.log.tpcb) <- c("predicted")
# Add predicted values to data.frame
che.log.tpcb$predicted <- fit.values.che.log.tpcb$predicted

# Plot prediction vs. observations, 1:1 line
ggplot(che.log.tpcb, aes(x = logtPCB, y = predicted)) +
  geom_point() +
  scale_x_continuous(limits = c(1, 75)) +
  scale_y_continuous(limits = c(1, 75)) +
  xlab(expression(bold("Observed concentration " *Sigma*"log"[10]*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"log"[10]*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3) +
  geom_abline(intercept = 2, slope = 1, col = "blue", size = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = -2, slope = 1, col = "blue", size = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotate('text', x = 25, y = 74,
           label = expression("Chesapeake Bay (R"^2*"= 0.66)"), colour = 'black', size = 4,
           fontface = 2)

# Plot residuals vs. predictions
plot(che.log.tpcb$predicted, res.che.log.tpcb,
     ylim = c(-25, 25),
     xlab = "Preditions",
     ylab = "Residual")
abline(0, 0)

# Selected sites -------------------------------------------------------------
# Remove predictions
che.tpcb.1 <- che.tpcb %>% select(-predicted)

# Near USGS station
che.tpcb.1 <- che.tpcb.1[che.tpcb$site %in% c('SusquehannaRiverConowingoBridge',
                                              'SusquehannaRiverMDPABorder',
                                              'PennysShoalParkIsland'), ]

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station
sitecheN1 <- "01578310" # SusquehannaRiverConowingoBridge, SusquehannaRiverMDPABorder, PennysShoalParkIsland
sitecheN2 <- "01482800" # Site Name: Delaware River at Reedy Island Jetty, DE just for water temp
# Codes to retrieve data
paramflow <- "00060" # discharge, ft3/s
paramtemp <- "00010" # water temperature, C
# Retrieve USGS data
flow <- readNWISdv(sitecheN1, paramflow,
                   min(che.tpcb.1$date), max(che.tpcb.1$date))
temp <- readNWISdv(sitecheN2, paramtemp,
                   min(che.tpcb.1$date), max(che.tpcb.1$date))
# Add USGS data to che.tpcb, matching dates
che.tpcb.1$flow <- flow$X_00060_00003[match(che.tpcb.1$date, flow$Date)]
che.tpcb.1$temp <- temp$X_00010_00003[match(che.tpcb.1$date, temp$Date)]
# Remove samples with temp = NA
che.tpcb.1 <- na.omit(che.tpcb.1)

# Add USGS data to che.log.tpcb, matching dates
che.log.tpcb.1$flow <- flow$X_.Primary.Stream.Flow._00060_00003[match(che.log.tpcb.1$date,
                                                                      flow$Date)]
che.log.tpcb.1$temp <- temp$X_00010_00003[match(che.log.tpcb.1$date, temp$Date)]
# Remove samples with temp = NA
che.log.tpcb.2 <- na.omit(che.log.tpcb.2)

# Plot
ggplot(che.tpcb.1, aes(y = tPCB,
                       x = format(date,'%Y%m'))) +
  geom_point() +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black"))

ggplot(che.tpcb.1, aes(x = factor(site), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 2012 - 2018 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0)

# Regression with selected sites
# (1.1) tPCB vs. flow (che.tpcb.1)
lr.che.tpcb.1.f <- lm(log10(tPCB) ~ flow, data = che.tpcb.1)
# See results
summary(lr.che.tpcb.1.f)
# Look at residuals
res <- resid(lr.che.tpcb.1.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. flow (che.log.tpcb.2)
lr.che.log.tpcb.1.f <- lm(logtPCB ~ flow, data = che.log.tpcb.1)
# See results
summary(lr.che.log.tpcb.1.f)
# Look at residuals
res <- resid(lr.che.log.tpcb.1.f) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB vs. water temperature (che.tpcb.2)
lr.che.tpcb.1.te <- lm(log10(tPCB) ~ temp, data = che.tpcb.1)
# See results
summary(lr.che.tpcb.1.te)
# Look at residuals
res <- resid(lr.che.tpcb.1.te) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. temperature (che.log.tpcb.1)
lr.che.log.tpcb.1.te <- lm(logtPCB ~ temp, data = che.log.tpcb.1)
# See results
summary(lr.che.log.tpcb.1.te)
# Look at residuals
res <- resid(lr.che.log.tpcb.1.te) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) MLR with selected sites
# (2.1) tPCB vs. time + season (che.tpcb)
mlr.che.tpcb.1 <- lm(log10(tPCB) ~ time + season + flow + temp, data = che.tpcb.1)
# See results
summary(mlr.che.tpcb.1)
# Look at residuals
res <- resid(mlr.che.tpcb.1) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (3) Perform Linear Mixed-Effects Model (LMEM) with selected sites
# (3.1) tPCB vs. time + season + flow + temp + site (che.tpcb.1)
tpcb <- che.tpcb.1$tPCB
time <- che.tpcb.1$time
site <- che.tpcb.1$site.code
season <- che.tpcb.1$season
flow <- che.tpcb.1$flow
tem <- che.tpcb.1$temp

lmem.che.tpcb.1 <- lmer(log10(tpcb) ~ 1 + time + season + flow + tem + (1|site),
                        REML = FALSE,
                        control = lmerControl(check.nobs.vs.nlev = "ignore",
                                              check.nobs.vs.rankZ = "ignore",
                                              check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.che.tpcb.1)
# Look at residuals
res.che.tpcb.1 <- resid(lmem.che.tpcb.1) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.che.tpcb.1, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.che.tpcb.1)
# Shapiro test
shapiro.test(res.che.tpcb.1)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.che.tpcb.1))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.che.tpcb.1))[1, 'R2c']

# Modeling plots
# (1) Get predicted values tpcb
fit.values.che.tpcb.1 <- as.data.frame(fitted(lmem.che.tpcb.1))
# Add column name
colnames(fit.values.che.tpcb.1) <- c("predicted")
# Add predicted values to data.frame
che.tpcb.1$predicted <- 10^(fit.values.che.tpcb.1$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(che.tpcb.1, aes(x = tPCB, y = predicted)) +
  geom_point() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 9),
        axis.title.x = element_text(face = "bold", size = 9)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "bl") +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3)

ggplot(che.tpcb.1, aes(x = tPCB, y = predicted)) +
  geom_point() +
  scale_x_log10(limits = c(1, 1e6)) +
  scale_y_log10(limits = c(1, 1e6)) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 100, y = 1000000,
           label = 'Chesapeake Bay', colour = 'black', size = 4,
           fontface = 2)
