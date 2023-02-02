## Water PCB concentrations analysis.
# Data were obtained from EPA and contractors from PCB Superfund
# sites in USA

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
{library(ggplot2)
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
}
# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor08052022.csv") # WaterDataConcentration

# Data preparation --------------------------------------------------------
# (1) All data, including 0s
# Remove metadata
wdc.1 <- subset(wdc, select = -c(SampleID:AroclorCongener))
# Remove Aroclor data
wdc.1 <- subset(wdc.1, select = -c(A1016:A1260))
# (2) Only consider congener data
wdc.cong <- subset(wdc, AroclorCongener == "Congener")
# Remove metadata
wdc.cong.1 <- subset(wdc.cong, select = -c(SampleID:AroclorCongener))
# Remove Aroclor data
wdc.cong.1 <- subset(wdc.cong.1, select = -c(A1016:A1260))

# Frequency analysis ------------------------------------------------------
# Just congener data
# Create a frequency detection plot
{wdc.cong.freq <- colSums(! is.na(wdc.cong.1) & (wdc.cong.1 !=0))/nrow(wdc.cong.1)
wdc.cong.freq <- data.frame(wdc.cong.freq)
colnames(wdc.cong.freq) <- c("PCB.frequency")
congener <- row.names(wdc.cong.freq)
wdc.cong.freq <- cbind(congener, wdc.cong.freq$PCB.frequency)
colnames(wdc.cong.freq) <- c("congener", "PCB.frequency")
wdc.cong.freq <- data.frame(wdc.cong.freq)
wdc.cong.freq$congener <- as.character(wdc.cong.freq$congener)
wdc.cong.freq$congener <- gsub('\\.', '+', wdc.cong.freq$congener) # replace dot for +
wdc.cong.freq$PCB.frequency <- as.numeric(as.character(wdc.cong.freq$PCB.frequency))
wdc.cong.freq$congener <- factor(wdc.cong.freq$congener,
                            levels = rev(wdc.cong.freq$congener)) # change the order to be plotted.
}

# Summary statistic of frequency of detection
summary(wdc.cong.freq$PCB.frequency)

# Frequency detection plot
ggplot(wdc.cong.freq, aes(x = 100*PCB.frequency, y = congener)) +
  geom_bar(stat = "identity", fill = "#66ccff") +
  geom_vline(xintercept = 100*mean(wdc.cong.freq$PCB.frequency),
             color = "red") +
  ylab("") +
  theme_bw() +
  xlim(c(0,100)) +
  theme(aspect.ratio = 20/5) +
  xlab(expression(bold("Frequency detection (%)"))) +
  theme(axis.text.x = element_text(face = "bold", size = 5),
        axis.title.x = element_text(face = "bold", size = 5)) +
  theme(axis.text.y = element_text(face = "bold", size = 3))

# Total Concentration Analysis --------------------------------------------
# Summary statistic of total PCB (congeners + Aroclor) in pg/L
summary(rowSums(wdc.1, na.rm = T))

# Histogram
hist(rowSums(wdc.1, na.rm = T))
hist(log10(rowSums(wdc.1, na.rm = T)))

# Total PCBs in 1 box plot
ggplot(wdc.1, aes(x = "", y = rowSums(wdc.1, na.rm = T))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(aspect.ratio = 14/2) +
  xlab(expression(bold(Sigma*"PCB (n = 5265)")))+
  ylab(expression(bold("Water Concentration 1990 - 2020 (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10,
                                    angle = 45, hjust = 1.8,
                                    vjust = 2)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l")

# Include 64 pg/L from EPA
ggplot(wdc.1, aes(x = "", y = rowSums(wdc.1, na.rm = T))) + 
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
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotation_logticks(sides = "l") +
  geom_hline(yintercept = 0.64*1000, color = "#9999CC",
             linewidth = 0.8) + # U.S. EPA Water Quality Criterion for Human Health from fish consumption, associated with an incremental cancer risk of 10−5
  geom_hline(yintercept = 0.064*1000, color = "#CC6666",
             linewidth = 0.8) # associated with an incremental cancer risk of 10−6.
  
# Individual congeners
summary(wdc.cong.1, na.rm = T, zero = T)
# Get the max value for each congener
cong.max <-as.numeric(sub('.*:', '',
                          summary(wdc.cong.1, na.rm = T,
                                  zero = T)[6,]))

# Obtain the median for each individual congener
cong.median <- as.numeric(sub('.*:',
                              '', summary(wdc.cong.1, na.rm = T,
                                          zero = T)[3,]))

# Individual PCB boxplot
ggplot(stack(wdc.cong.1), aes(x = ind, y = values)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_boxplot(width = 0.6, outlier.colour = "#66ccff", col = "#66ccff",
               outlier.shape = 1) +
  scale_x_discrete(labels = wdc.cong.freq$congener) + # use to change the "." to "+"
  theme_bw() +
  theme(aspect.ratio = 25/135) +
  xlab(expression("")) +
  ylab(expression(bold("PCB congener concentration (pg/L)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8,
                                   color = "black"),
        axis.title.y = element_text(face = "bold", size = 8,
                                    color = "black")) +
  theme(axis.text.x = element_text(face = "bold", size = 6,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.6, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# Spatial Plots and Analysis ----------------------------------------------
# Modify x-axis
# States. Needs work
sites <- c("CA", "DE", "ID", "IN", "MA", "MD", "MI", "MO",
           "MT", "NM", "NY", "OH", "OR", "TX", "VA", "WA", "WI")

# Total PCBs
ggplot(wdc, aes(x = factor(StateSampled, levels = sites),
                y = rowSums(wdc[, c(14:117)],  na.rm = T))) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB 1990 - 2019 (pg/L)"))))) +
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


# Selected StateSampled and individual PCB congener
wdc.pcb.sp <- subset(wdc, select = c(StateSampled, PCB4.10))

# Plot
ggplot(wdc.pcb.sp, aes(x = factor(StateSampled, levels = sites),
                   y = PCB4.10)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concentration",
                            paste("PCB 4+10 1990 - 2019 (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 1, col = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  geom_hline(yintercept = 0.03, color = "#cc0000") # median 

# Regression analysis and plots---------------------------------------------
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

# Plots -------------------------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(wdc.tpcb$tPCB)
hist(log10(wdc.tpcb$tPCB))
# (1.2) log.tPCB
hist(wdc.tpcb$logtPCB)
hist(log10(wdc.tpcb$logtPCB))

# (2) Time trend plots
# (2.1) tPCB
ggplot(wdc.tpcb, aes(y = tPCB,
                     x = format(date,'%Y'))) +
  geom_point(shape = 1, cex = 1.2, col = "#66ccff") +
  theme(aspect.ratio = 5/20) +
  xlab("") +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB 1990 - 2019 (pg/L)"))))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black")) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 11))

# (2.2) log.tPCB
ggplot(wdc.tpcb, aes(y = logtPCB,
                         x = format(date,'%Y'))) +
  geom_point(shape = 1, cex = 1.2, col = "#66ccff") +
  xlab("") +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB 1990 - 2019 (pg/L)"))))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1,
                                   color = "black"))

# (3) Seasonality
# (3.1) tPCB
ggplot(wdc.tpcb, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB 1990 - 2019 (pg/L)"))))) +
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
ggplot(wdc.tpcb, aes(x = season, y = logtPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB 1990 - 2019 (pg/L)"))))) +
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
# Get variables
tpcb <- wdc.tpcb$tPCB
log.tpcb <- wdc.tpcb$logtPCB
time <- wdc.tpcb$time
site <- wdc.tpcb$site.code
season <- wdc.tpcb$season
# (1) Perform linear regression (lr)
# (1.1) tPCB vs. time
lr.wdc.tpcb.t <- lm(log10(tpcb) ~ time)
# See results
summary(lr.wdc.tpcb.t)
# Look at residuals
res <- resid(lr.wdc.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) log.tPCB vs. time
lr.wdc.log.tpcb.t <- lm(log.tpcb ~ time)
# See results
summary(lr.wdc.log.tpcb.t)
# Look at residuals
res <- resid(lr.wdc.log.tpcb.t) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.3) tPCB vs. season
lr.wdc.tpcb.s <- lm(log10(tpcb) ~ season)
# See results
summary(lr.wdc.tpcb.s)
# Look at residuals
res <- resid(lr.wdc.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.4) log.tPCB vs. season
lr.wdc.log.tpcb.s <- lm(log.tpcb ~ season)
# See results
summary(lr.wdc.log.tpcb.s)
# Look at residuals
res <- resid(lr.wdc.log.tpcb.s) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) MLR
# (2.1) tPCB vs. time + season
mlr.wdc.tpcb <- lm(log10(tpcb) ~ time + season)
# See results
summary(mlr.wdc.tpcb)
# Look at residuals
res <- resid(mlr.wdc.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2.2) log.tPCB vs. time + season
mlr.wdc.log.tpcb <- lm(log.tpcb ~ time + season)
# See results
summary(mlr.wdc.log.tpcb)
# Look at residuals
res <- resid(mlr.wdc.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (3) Perform Linear Mixed-Effects Model (lme)
# (3.1) tPCB vs. time + season + site
lmem.wdc.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE = "ignore"))

# See results
summary(lmem.wdc.tpcb)
# Look at residuals
res.wdc.tpcb <- resid(lmem.wdc.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.wdc.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.wdc.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.wdc.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.wdc.tpcb))[1, 'R2c']

# Extract coefficient values
time.coeff <- summary(lmem.wdc.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lmem.wdc.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -log(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# Modeling plots
# (1) Get predicted values tpcb
fit.values.wdc.tpcb <- as.data.frame(fitted(lmem.wdc.tpcb))
# Add column name
colnames(fit.values.wdc.tpcb) <- c("predicted")
# Add predicted values to data.frame
wdc.tpcb$predicted <- 10^(fit.values.wdc.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(wdc.tpcb, aes(x = tPCB, y = predicted)) +
  geom_point() +
  scale_x_log10(limits = c(0.1, 1e7)) +
  scale_y_log10(limits = c(0.1, 1e7)) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3) +
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl")

# Plot residuals vs. predictions
plot(wdc.tpcb$predicted, res.wdc.tpcb)
abline(0, 0)

# (3.2) log.tPCB vs. time + season + site (wdc.log.tpcb)
lmem.wdc.log.tpcb <- lmer(log.tpcb ~ 1 + time + season + season + (1|site),
                          REML = FALSE,
                          control = lmerControl(check.nobs.vs.nlev = "ignore",
                                                check.nobs.vs.rankZ = "ignore",
                                                check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.wdc.log.tpcb)
# Look at residuals
res.wdc.log.tpcb <- resid(lmem.wdc.log.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res.wdc.log.tpcb, main = "log10(C)")
# Add a straight diagonal line to the plot
qqline(res.wdc.log.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res.wdc.log.tpcb, 'pnorm')
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.wdc.log.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.wdc.log.tpcb))[1, 'R2c']

# (2) Get predicted values log.tpcb
fit.values.wdc.log.tpcb <- as.data.frame(fitted(lmem.wdc.log.tpcb))
# Add column name
colnames(fit.values.wdc.log.tpcb) <- c("predicted")
# Add predicted values to data.frame
wdc.tpcb$predictedlog <- (fit.values.wdc.log.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(wdc.tpcb, aes(x = logtPCB, y = predictedlog)) +
  geom_point() +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", size = 1.3) +
  theme_bw() +
  theme(aspect.ratio = 15/15)

