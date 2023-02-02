## Pre water PCB concentrations data analysis

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
install.packages("zoo")

# Load libraries
{library(ggplot2)
library(scales) # function trans_breaks
#library(gridExtra)
#library(tidyverse)
library(stringr) # str_detect
library(robustbase) # function colMedians
library(dplyr) # performs %>%
library(tibble) # adds a column
library(lme4) # performs lme
library(MuMIn) # gets Rs from lme
library(lmerTest) # gets the p-value from lme
library(zoo) # yields seasons
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor08052022.csv")

# (I) Entire data ---------------------------------------------------------
# Prepare data -----------------------------------------------------------
# Remove samples (rows) with total PCBs  = 0
wdc.2 <- wdc[!(rowSums(wdc[, c(14:117)], na.rm = TRUE)==0),]
# Calculate total PCB
wdc.tpcb <- rowSums(wdc.2[, c(14:117)], na.rm = T)
# Calculate total log PCB
# Remove metadata
wdc.3 <- subset(wdc.2, select = -c(SampleID:AroclorCongener))
# Remove Aroclor data
wdc.3 <- subset(wdc.3, select = -c(A1016:A1260))
# Log10 individual PCBs 
wdc.log <- log10(wdc.3)
# Replace -inf to NA
wdc.log <- do.call(data.frame,
                   lapply(wdc.log,
                          function(x) replace(x, is.infinite(x), NA)))

# Sum individual log 10 PCBs
wdc.log.tpcb <- rowSums(wdc.log, na.rm = T)
# Change date format
wdc.2$SampleDate <- as.Date(wdc.2$SampleDate, format = "%m/%d/%y")
# Calculate sampling time
time.day <- data.frame(as.Date(wdc.2$SampleDate) - min(as.Date(wdc.2$SampleDate)))
# Create individual code for each site sampled
site.numb <- wdc.2$SiteID %>% as.factor() %>% as.numeric
# Include season
yq.s <- as.yearqtr(as.yearmon(wdc.2$SampleDate, "%m/%d/%Y") + 1/12)
season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                   labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
# Create data frame
wdc.tpcb <- cbind(factor(wdc.2$SiteID), wdc.2$SampleDate,
                  wdc.2$Latitude, wdc.2$Longitude, as.matrix(wdc.tpcb),
                  as.matrix(wdc.log.tpcb), data.frame(time.day),
                  site.numb, season.s)
# Add column names
colnames(wdc.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                        "tPCB", "logtPCB", "time", "site.code", "season")

# Number of sites
num.site <- aggregate(wdc.tpcb$SiteID, by = list(wdc.tpcb$SiteID), FUN = length)
# Add names to the columns
colnames(num.site) <- c("Site Name", "#")
# Number of site replicates 
site.replic <- aggregate(wdc.tpcb$SiteID, by = list(wdc.tpcb$SiteID), FUN = length)
# Add names to the columns
colnames(site.replic) <- c("SiteID", "#")
median(site.replic$`#`)
# Number of locations replicates
location.replic <- aggregate(wdc$LocationName, by = list(wdc$LocationName),
                             FUN = length)
# Add names to the columns
colnames(location.replic) <- c("LocationID", "#")
mean(location.replic$`#`)

# General plots -------------------------------------------------------

# Reviewed until here!
# (1) Histograms
# (1.1) tpcb
hist(tpcb$tPCB)
hist(log10(tpcb$tPCB))
# (1.2) log.tPCB
hist(log.tpcb$logtPCB)
hist(log10(log.tpcb$logtPCB))

# (2) One box plot
# (2.1) tPCB
ggplot(tpcb, aes(x = "", y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(aspect.ratio = 14/2) +
  xlab(expression(bold(Sigma*"PCB (n = 5300)")))+
  ylab(expression(bold("Water Concentration 1990 - 2020 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10,
                                    angle = 45, hjust = 1.8,
                                    vjust = 2)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l")

# (2.2) log.tPCB
ggplot(log.tpcb, aes(x = "", y = logtPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(aspect.ratio = 14/2) +
  xlab(expression(bold(Sigma*"PCB (n = 6098)")))+
  ylab(expression(bold("Water Concentration 1990 - 2020 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10,
                                    angle = 45, hjust = 1.8,
                                    vjust = 2)) +
  theme(axis.ticks = element_line(size = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l")

# (3) Box plots per site
# (3.1) tPCB
ggplot(tpcb, aes(x = factor(Site), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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
ggplot(log.tpcb, aes(x = factor(Site), y = logtPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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

# (4) Time series
# (4.1) tPCB
ggplot(tpcb, aes(x = format(date,'%Y'), y = tPCB)) +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  geom_hline(yintercept = 0.64*1000, color = "#9999CC",
             size = 0.8) + # U.S. EPA Water Quality Criterion for Human Health from fish consumption, associated with an incremental cancer risk of 10−5
  geom_hline(yintercept = 0.064*1000, color = "#CC6666",
             size = 0.8) # associated with an incremental cancer risk of 10−6.

# (4.2) log.tPCB
ggplot(log.tpcb, aes(x = format(date,'%Y'), y = logtPCB)) +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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

# (5) Seasonality
# (5.1) tPCB
ggplot(tpcb, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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

# (5.2) log.tPCB
ggplot(log.tpcb, aes(x = season, y = logtPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1990 - 2020 (pg/L)"))) +
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

# General regressions -----------------------------------------------------
# Remove tPCB = 0
tpcb.1 <- tpcb[!(tpcb$tPCB == 0), ]
log.tpcb.1 <- log.tpcb[!(log.tpcb$logtPCB == 0), ]

# (1) Perform linear regression (lr)
# (1.1) tPCB vs. time
lr.tpcb <- lm(log10(tPCB) ~ time, data = tpcb.1)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.tpcb <- lm(logtPCB ~ time, data = log.tpcb.1)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB + 1 vs. season
lr.tpcb <- lm(log10(tPCB) ~ time + season, data = tpcb.1)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.tpcb <- lm(logtPCB ~ season, data = log.tpcb.1)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) MLR
# (2.1) tPCB vs. time + season (tpcb.1)
mlr.tpcb <- lm(log10(tPCB) ~ time + season, data = tpcb.1)
# See results
summary(mlr.tpcb)
# Look at residuals
res <- resid(mlr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2.2) log.tPCB vs. time + season + flow + temp (log.tpcb.1)
mlr.log.tpcb <- lm(logtPCB ~ time + season,
                       data = log.tpcb.1)
# See results
summary(mlr.log.tpcb)
# Look at residuals
res <- resid(mlr.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (3) Perform Linear Mixed-Effects Model (LMEM)
# (3.1) tPCB vs. time + season + site (tpcb.1)
# Create individual code for each site sampled
site.numb <- tpcb.1$LocationID %>% as.factor() %>% as.numeric
tpcb.1$site.code <- site.numb
tpcb <- tpcb.1$tPCB
time <- tpcb.1$time
site <- tpcb.1$site.code
season <- tpcb.1$season

lmem.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.tpcb)
# Look at residuals
res.tpcb <- resid(lmem.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.tpcb, main = expression(paste("Normal Q-Q Plot",
                                         " (log"[10], Sigma, "PCB)")))
# Add a straight diagonal line to the plot
qqline(res.tpcb, col = "blue", lwd = 2)
# Shapiro test
shapiro.test(res.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.tpcb))[1, 'R2c']

# Extract coefficient values
time.coeff <- summary(lmem.che.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lmem.che.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -log(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# (3.2) log.tPCB vs. time + season + site (che.log.tpcb)
site.numb.log <- log.tpcb.1$LocationID %>% as.factor() %>% as.numeric
log.tpcb.1$site.code <- site.numb.log
log.tpcb <- log.tpcb.1$logtPCB
time <- log.tpcb.1$time
site <- log.tpcb.1$site.code
season <- log.tpcb.1$season

lmem.log.tpcb <- lmer(log.tpcb ~ 1 + time + season + (1|site),
                          REML = FALSE,
                          control = lmerControl(check.nobs.vs.nlev = "ignore",
                                                check.nobs.vs.rankZ = "ignore",
                                                check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.log.tpcb)
# Look at residuals
res.log.tpcb <- resid(lmem.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.log.tpcb, main = expression(paste("Normal Q-Q Plot", " (", Sigma,
                                         "log"[10]*"PCB)")))
# Add a straight diagonal line to the plot
qqline(res.log.tpcb, col = "blue", lwd = 2)
# Shapiro test
shapiro.test(res.log.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res.log.tpcb, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.log.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.log.tpcb))[1, 'R2c']

# Predictions -------------------------------------------------------------
# Modeling plots
# (1) Get predicted values tpcb
fit.values.tpcb <- as.data.frame(fitted(lmem.tpcb))
# Add column name
colnames(fit.values.tpcb) <- c("predicted")
# Add predicted values to data.frame
tpcb.1$predicted <- 10^(fit.values.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(tpcb.1, aes(x = tPCB, y = predicted)) +
  geom_point() +
  scale_x_log10(limits = c(0.5, 1e7)) +
  scale_y_log10(limits = c(0.5, 1e7)) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 0.8) + # 1:1 line
  geom_abline(intercept = 0.5, slope = 1, col = "blue", size = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = -0.5, slope = 1, col = "blue", size = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 25, y = 7000000,
           label = "All samples (R2 = 0.89)", colour = 'black', size = 4,
           fontface = 2)

# Plot residuals vs. predictions
plot(tpcb.1$predicted, res.tpcb,
     xlab = "Preditions",
     ylab = "Residual")
abline(0, 0)

# (2) Get predicted values log.tpcb
fit.values.log.tpcb <- as.data.frame(fitted(lmem.log.tpcb))
# Add column name
colnames(fit.values.log.tpcb) <- c("predicted")
# Add predicted values to data.frame
log.tpcb.1$predicted <- fit.values.log.tpcb$predicted

# Plot prediction vs. observations, 1:1 line
ggplot(log.tpcb.1, aes(x = logtPCB, y = predicted)) +
  geom_point() +
  scale_x_continuous(limits = c(-150, 350)) +
  scale_y_continuous(limits = c(-150, 350)) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotate('text', x = -50, y = 350,
           label = "All samples (R2 = 0.85)", colour = 'black', size = 4,
           fontface = 2)
  
# Plot residuals vs. predictions
plot(log.tpcb.1$predicted, res.log.tpcb,
     xlab = "Preditions",
     ylab = "Residual")
abline(0, 0)

