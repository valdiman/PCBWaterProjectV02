## Water PCB concentrations data analysis per site
# Kalamazoo River
# Aroclors 1242, 1254 and 1260

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

# Select Kalamazoo data ---------------------------------------------------
kal <- wdc[str_detect(wdc$SiteName, 'Kalamazoo River'),]
# Superfund site from Morrow Dam (Kalamazoo River) to Lake Michigan
# and 30 miles of Portage Creek (south), Cork St and Portage Creek Cork St sites
# Dredging occurred at Plainwell Dam site.

# Data preparation --------------------------------------------------------
# Calculate total PCB
tpcb.kal <- rowSums(kal[, c(14:117)], na.rm = T)
# Change date format
kal$SampleDate <- as.Date(kal$SampleDate, format = "%m/%d/%y")
# Calculate sampling time
time.day <- data.frame(as.Date(kal$SampleDate) - min(as.Date(kal$SampleDate)))
# Create individual code for each site sampled
site.numb <- kal$LocationID %>% as.factor() %>% as.numeric
# Include season
yq.s <- as.yearqtr(as.yearmon(kal$SampleDate, "%m/%d/%Y") + 1/12)
season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                   labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
# Create data frame
kal.tpcb <- cbind(factor(kal$LocationID), kal$SampleDate,
                kal$Latitude, kal$Longitude, as.matrix(tpcb.kal),
                data.frame(time.day), site.numb, season.s)
# Add column manes
colnames(kal.tpcb) <- c("LocationID", "date", "Latitude", "Longitude",
                        "tPCB", "time", "site.code", "season")

# Get coordinates per site
kal.location <- kal.tpcb[c('LocationID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
kal.location <- aggregate(tPCB ~ LocationID + Latitude + Longitude,
                            data = kal.location, mean)

# (2) Calculate total log PCB
# Remove metadata
kal.log <- subset(kal, select = -c(SampleID:AroclorCongener))
# Remove Aroclor data
kal.log <- subset(kal.log, select = -c(A1016:A1260))
# Log 10 individual PCBs 
kal.log <- log10(kal.log)
# Replace -inf to NA
kal.log <- do.call(data.frame,
                       lapply(kal.log,
                              function(x) replace(x, is.infinite(x), NA)))
# Sum individual log 10 PCBs
kal.log.tpcb <- rowSums(kal.log, na.rm = T)
# Generate data.frame for analysis and plots
kal.log.tpcb <- cbind(factor(kal$LocationID), kal$SampleDate,
                    as.matrix(kal.log.tpcb), data.frame(time.day),
                    site.numb, season.s)
colnames(kal.log.tpcb) <- c("LocationID", "date", "logtPCB", "time",
                          "site.code", "season")

# General plots -------------------------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(kal.tpcb$tPCB)
hist(log10(kal.tpcb$tPCB))
# (1.2) log.tPCB
hist(kal.log.tpcb$logtPCB)

# (2) Time trend plots
# (2.1) tPCB
ggplot(kal.tpcb, aes(y = tPCB,
                     x = format(date,'%Y'))) +
  geom_point(shape = 1, col = "#66ccff") +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Conncetration",
                            paste(Sigma*"PCB 1993 - 2010 (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9))

# (2.2) log.tPCB
ggplot(kal.log.tpcb, aes(y = logtPCB,
                         x = format(date,'%Y'))) +
  geom_point() +
  xlab("") +
  theme_bw() +
  theme(aspect.ratio = 5/15)

# (3) Seasonality
# (3.1) tPCB
ggplot(kal.tpcb, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1993 - 2010 (pg/L)"))) +
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
ggplot(kal.log.tpcb, aes(x = season, y = logtPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1993 - 2010 (pg/L)"))) +
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
# (4.1) tPCB
ggplot(kal.tpcb, aes(x = factor(LocationID), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1993 - 2010 (pg/L)"))) +
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

# Remove site -------------------------------------------------------------
# Remove site PlainwellDam. Dredging = WCPCB-KAL023
kal.tpcb.1 <- subset(kal.tpcb, LocationID != c("WCPCB-KAL023"))
kal.log.tpcb.1 <- subset(kal.log.tpcb, LocationID != c("WCPCB-KAL023"))

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Kalamazoo River
siteKalN1 <- "04108660" # KALAMAZOO RIVER AT NEW RICHMOND, MI
siteKalN2 <- "04106000" # KALAMAZOO RIVER AT COMSTOCK, MI
# Codes to retrieve data
paramflow <- "00060" # discharge, ft3/s
# paramtemp <- "00010" # water temperature, C Not available

# Flow (ft3/s)
flow.1 <- readNWISdv(siteKalN1, paramflow,
                     min(kal.tpcb.1$date), max(kal.tpcb.1$date))
flow.2 <- readNWISdv(siteKalN2, paramflow,
                     min(kal.tpcb.1$date), max(kal.tpcb.1$date))

kal.tpcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(kal.tpcb.1$date,
                                                flow.1$Date)]
kal.tpcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(kal.tpcb.1$date,
                                                      flow.2$Date)]
# Create flow, flow.3
kal.tpcb.1$flow.3 <- ifelse(is.na(kal.tpcb.1$flow.1) == TRUE,
                            kal.tpcb.1$flow.2/0.46, kal.tpcb.1$flow.1)
# Remove samples with flow.1 = NA
kal.tpcb.2 <- na.omit(kal.tpcb.1)

# 2nd plots -----------------------------------------------------------
# (1.1) tPCB
hist(kal.tpcb.1$tPCB)
hist(log10(kal.tpcb.1$tPCB))
# (1.2) log.tPCB
hist(kal.log.tpcb.1$logtPCB)

# (2) Time trend plots
# (2.1) tPCB
ggplot(kal.tpcb.1, aes(y = tPCB,
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
ggplot(kal.log.tpcb.1, aes(y = logtPCB,
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
ggplot(kal.tpcb.1, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1993 - 2010 (pg/L)"))) +
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
ggplot(kal.log.tpcb.1, aes(x = season, y = logtPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold("Water Conncetration " *Sigma*"PCB 1993 - 2010 (pg/L)"))) +
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

# (3.3) Flow plots
ggplot(kal.tpcb.2, aes(y = tPCB, x = flow.1)) +
  geom_point() +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15)

ggplot(kal.tpcb.1, aes(y = tPCB, x = flow.3)) +
  geom_point() +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15)

# (3.4) box plot
ggplot(kal.tpcb.1, aes(x = "", y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(aspect.ratio = 14/2) +
  xlab(expression(bold(Sigma*"PCB (n = 556)")))+
  ylab(expression(bold("Water Concentration 1993 - 2010 (pg/L)"))) +
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

# Regressions -------------------------------------------------------------
# (1) Perform linear regression (lr)
# (1.1) tPCB vs. time
lr.kal.tpcb.t <- lm(log10(tPCB) ~ time, data = kal.tpcb.1)
# See results
summary(lr.kal.tpcb.t)
# Look at residuals
res <- resid(lr.kal.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.kal.log.tpcb.t <- lm(logtPCB ~ time, data = kal.log.tpcb.1)
# See results
summary(lr.kal.log.tpcb.t)
# Look at residuals
res <- resid(lr.kal.log.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB vs. season
lr.kal.tpcb.s <- lm(log10(tPCB) ~ season, data = kal.tpcb.1)
# See results
summary(lr.kal.tpcb.s)
# Look at residuals
res <- resid(lr.kal.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.kal.log.tpcb.s <- lm(logtPCB ~ season, data = kal.log.tpcb.1)
# See results
summary(lr.kal.log.tpcb.s)
# Look at residuals
res <- resid(lr.kal.log.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) MLR
# (2.1) tPCB vs. time + season
mlr.kal.tpcb <- lm(log10(tPCB) ~ time + season, data = kal.tpcb.1)
# See results
summary(mlr.kal.tpcb)
# Look at residuals
res <- resid(mlr.kal.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2.2) log.tPCB vs. time + season 
mlr.kal.log.tpcb <- lm(logtPCB ~ time + season, data = kal.log.tpcb.1)
# See results
summary(mlr.kal.log.tpcb)
# Look at residuals
res <- resid(mlr.kal.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (3) Perform Linear Mixed-Effects Model (LMEM)
# (3.1) tPCB vs. time + season + flow + temp + site (kal.tpcb.2)
tpcb <- kal.tpcb.1$tPCB
time <- kal.tpcb.1$time
site <- kal.tpcb.1$site.code
season <- kal.tpcb.1$season
flow <- kal.tpcb.1$flow.3

lmem.kal.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + season + flow + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.kal.tpcb)
# Look at residuals
res.kal.tpcb <- resid(lmem.kal.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.kal.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.kal.tpcb)
# Shapiro test
shapiro.test(res.kal.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res.kal.tpcb, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.kal.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.kal.tpcb))[1, 'R2c']

# Extract coefficient values
time.coeff <- summary(lmem.kal.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lmem.kal.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -log(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# (3.2) log.tPCB vs. time + season + flow + temp + site (kal.log.tpcb.2)
log.tpcb <- kal.log.tpcb.1$logtPCB
time <- kal.log.tpcb.1$time
site <- kal.log.tpcb.1$site.code
season <- kal.log.tpcb.1$season
flow <- kal.log.tpcb.1$flow
tem <- kal.log.tpcb.1$temp

lmem.kal.log.tpcb <- lmer(log.tpcb ~ 1 + time + season + season + flow + tem + (1|site),
                          REML = FALSE,
                          control = lmerControl(check.nobs.vs.nlev = "ignore",
                                                check.nobs.vs.rankZ = "ignore",
                                                check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.kal.log.tpcb)
# Look at residuals
res.kal.log.tpcb <- resid(lmem.kal.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.kal.log.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.kal.log.tpcb)
# Shapiro test
shapiro.test(res.kal.log.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res.kal.log.tpcb, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.kal.log.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.kal.log.tpcb))[1, 'R2c']

# Plots
# Get predicted values
fit.values <- as.data.frame(fitted(lmem.kal.tpcb))
# Add column name
colnames(fit.values) <- c("predicted")
# Add predicted values to data.frame kal.tpcb.3
kal.tpcb.1$predicted <- fit.values$predicted

# Plot prediction vs. observations, 1:1 line
ggplot(kal.tpcb.1, aes(x = tPCB, y = 10^predicted)) +
  geom_point() +
  scale_x_log10(limits = c(10, 1e6)) +
  scale_y_log10(limits = c(10, 1e6)) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3) +
  geom_abline(intercept = 0.5, slope = 1, col = "blue", size = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = -0.5, slope = 1, col = "blue", size = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 10^2.5, y = 10^6,
           label = expression("Kalamazoo River (R"^2*"= 0.82)"), colour = 'black', size = 4,
           fontface = 2)


# Plot residuals vs. predictions
plot(log10(kal.tpcb.1$predicted), res.kal.tpcb,
     ylim = c(-1, 1),
     xlab = "Preditions",
     ylab = "Residual")
abline(0, 0)
