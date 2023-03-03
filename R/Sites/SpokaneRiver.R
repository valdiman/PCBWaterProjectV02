## Water PCB concentrations data analysis per site
## Spokane River

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
{
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
  library(reshape)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor08052022.csv")

# Select Spokane River data ---------------------------------------------------
spo.0 <- wdc[str_detect(wdc$LocationName, 'Spokane River'),]

# Data preparation --------------------------------------------------------
{
  # Remove samples (rows) with total PCBs  = 0
  spo.1 <- spo.0[!(rowSums(spo.0[, c(14:117)], na.rm = TRUE)==0),]
  # Calculate total PCB
  tpcb.spo <- rowSums(spo.1[, c(14:117)], na.rm = T)
  # Change date format
  spo.1$SampleDate <- as.Date(spo.1$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(spo.1$SampleDate) - min(as.Date(spo.1$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- spo.1$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(spo.1$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  spo.tpcb <- cbind(factor(spo.1$SiteID), spo.1$SampleDate,
                    spo.1$Latitude, spo.1$Longitude, as.matrix(tpcb.spo),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(spo.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
spo.location <- spo.tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
spo.location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                          data = spo.location, mean)

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(spo.tpcb$tPCB)
hist(log10(spo.tpcb$tPCB))

# (2) Time trend plots
ggplot(spo.tpcb, aes(y = tPCB,
                     x = format(date,'%Y-%m'))) +
  geom_point(shape = 21, size = 2, fill = "#66ccff") +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotate("text", x = 4.5, y = 10^4.2, label = "Spokane River",
           size = 3)

# (3) Seasonality
ggplot(spo.tpcb, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotate("text", x = 1, y = 10^4.2, label = "Spokane River",
           size = 3)

# (4) Sites
ggplot(spo.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotate("text", x = 2.5, y = 10^4.2, label = "Spokane River",
           size = 3)

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Spokane River
siteSpoN1 <- "12417650" # SPOKANE RIVER BLW BLACKWELL NR COEUR D ALENE ID
siteSpoN2 <- "12419000" # Spokane River near Post Falls, ID
siteSpoN3 <- "12422500" # Spokane River at Spokane, WA
siteSpoN4 <- "12424000" # Hangman Creek at Spokane, WA
siteSpoN5 <- "12422000"
# Codes to retrieve data
paramflow <- "00060" # discharge, ft3/s
#paramtemp <- "00010" # water temperature, C No data
# Retrieve USGS data
flow.1 <- readNWISdv(siteSpoN1, paramflow,
                     min(spo.tpcb$date), max(spo.tpcb$date))
flow.2 <- readNWISdv(siteSpoN2, paramflow,
                     min(spo.tpcb$date), max(spo.tpcb$date))
flow.3 <- readNWISdv(siteSpoN3, paramflow,
                     min(spo.tpcb$date), max(spo.tpcb$date))
flow.4 <- readNWISdv(siteSpoN4, paramflow,
                     min(spo.tpcb$date), max(spo.tpcb$date))

# Add USGS data to spo.tpcb, matching dates
spo.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(spo.tpcb$date,
                                              flow.1$Date)]
spo.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(spo.tpcb$date,
                                              flow.2$Date)]
spo.tpcb$flow.3 <- 0.03*flow.3$X_00060_00003[match(spo.tpcb$date,
                                              flow.3$Date)]
spo.tpcb$flow.4 <- 0.03*flow.4$X_00060_00003[match(spo.tpcb$date,
                                              flow.4$Date)]

# Remove site -------------------------------------------------------------
## Sample sites not located in the Spokane River
## Coeur d'Alene WWTP, Post Falls WWTP, Liberty Lake, Kaiser Aluminum
spo.tpcb.2 <- subset(spo.tpcb, SiteID != c("WCPCB-SPR002")) #Coeur d'Alene WWTP
spo.tpcb.2 <- subset(spo.tpcb.2, SiteID != c("WCPCB-SPR005")) # Inland Empire Paper
spo.tpcb.2 <- subset(spo.tpcb.2, SiteID != c("WCPCB-SPR006")) # Kaiser Aluminum
spo.tpcb.2 <- subset(spo.tpcb.2, SiteID != c("WCPCB-SPR008")) #Liberty Lake
spo.tpcb.2 <- subset(spo.tpcb.2, SiteID != c("WCPCB-SPR011")) # Post Falls WWTP
spo.tpcb.2 <- subset(spo.tpcb.2, SiteID != c("WCPCB-SPR0012")) # Regional WRF

# (1) Histograms
hist(spo.tpcb.2$tPCB)
hist(log10(spo.tpcb.2$tPCB))

# (2) Time trend plots
ggplot(spo.tpcb.2, aes(y = tPCB,
                     x = format(date,'%Y-%m'))) +
  geom_point(shape = 21, size = 2, fill = "#66ccff") +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotate("text", x = 4.5, y = 10^4.2, label = "Spokane River",
           size = 3)

# (3) Sites
ggplot(spo.tpcb.2, aes(x = factor(SiteID), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(linewidth = 0.8, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l") +
  geom_jitter(position = position_jitter(0.3), cex = 1.2,
              shape = 21, fill = "#66ccff") +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0) +
  annotate("text", x = 2.5, y = 10^4.2, label = "Spokane River",
           size = 3)

# Select sites
spo.tpcb.PF <- subset(spo.tpcb, SiteID == "WCPCB-SPR010") # flow.2 PostFalls
spo.tpcb.WRF <- subset(spo.tpcb, SiteID == "WCPCB-SPR013") # flow.3 Spokane WRF

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- spo.tpcb.2$tPCB
time <- spo.tpcb.2$time
site <- spo.tpcb.2$site.code
season <- spo.tpcb.2$season
flow <- spo.tpcb.2$flow.4
# tPCB vs. time + season + flow + temp + site
lme.spo.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + flow + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.spo.tpcb)
# Look at residuals
{
  res.spo.tpcb <- resid(lme.spo.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.spo.tpcb, main = "log10(C)")
  qqnorm(res.spo.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.spo.tpcb)
}
# Shapiro test
shapiro.test(res.spo.tpcb)
# Random effect site Std Dev
RandonEffectSiteStdDev <- as.data.frame(VarCorr(lme.spo.tpcb))[1,'sdcor']
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lme.spo.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lme.spo.tpcb))[1, 'R2c']
# Extract coefficient values
time.coeff <- summary(lme.spo.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lme.spo.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -ln(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.spo.tpcb <- as.data.frame(fitted(lme.spo.tpcb))
# Add column name
colnames(fit.lme.values.spo.tpcb) <- c("predicted")
# Add predicted values to data.frame
spo.tpcb$predicted.1 <- 10^(fit.lme.values.spo.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(spo.tpcb, aes(x = tPCB, y = predicted.1)) +
  geom_point(shape = 21, size = 3, fill = "#66ccff") +
  scale_y_log10(limits = c(10, 10^4.5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10, 10^4.5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", linewidth = 1.3) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 50, y = 10^4.3,
           label = expression(atop("Spokane River (R"^2*"= 0.65)",
                                   paste("t"[1/2]*" = 11 ± 2 (yr)"))),
           size = 3, fontface = 2)

# Plot residuals vs. predictions
{
  plot(log10(spo.tpcb$predicted), res.spo.tpcb,
       points(log10(spo.tpcb$predicted.1), res.spo.tpcb, pch = 16, 
              col = "#66ccff"),
       xlim = c(1.5, 4),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  abline(0, 0)
  abline(h = c(-1, 1), col = "grey")
  abline(v = seq(1.5, 4, 0.5), col = "grey")
  }

# Estimate a factor of 2 between observations and predictions
spo.tpcb$factor2 <- spo.tpcb$tPCB/spo.tpcb$predicted.1
factor2.tpcb <- nrow(spo.tpcb[spo.tpcb$factor2 > 0.5 & spo.tpcb$factor2 < 2,
                                ])/length(spo.tpcb[,1])*100
