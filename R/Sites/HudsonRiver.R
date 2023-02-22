## Water PCB concentrations data analysis per site
## Hudson River
## Only Linear Mixed-Effects Model (lme) and
## log10SumPCB

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
install.packages("reshape")

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

# Select Hudson River data ---------------------------------------------------
hud.0 <- wdc[str_detect(wdc$LocationName, 'Hudson River'),]
# PCBs were discharged to the river from the General Electric
# (GE) manufacturing plants in Hudson Falls and Fort Edward, NY
# Dredging from 2009 to 2015
# https://www.epa.gov/system/files/documents/2021-08/hudson_summer2021_floodplainrifs_factsheet_final.pdf

# Data preparation --------------------------------------------------------
{
  # Remove samples (rows) with total PCBs  = 0
  hud.1 <- hud.0[!(rowSums(hud.0[, c(14:117)], na.rm = TRUE)==0),]
  # Calculate total PCB
  tpcb.hud <- rowSums(hud.1[, c(14:117)], na.rm = T)
  # Change date format
  hud.1$SampleDate <- as.Date(hud.1$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(hud.1$SampleDate) - min(as.Date(hud.1$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- hud.1$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hud.1$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  hud.tpcb <- cbind(factor(hud.1$SiteID), hud.1$SampleDate,
                    hud.1$Latitude, hud.1$Longitude, as.matrix(tpcb.hud),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(hud.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
hud.location <- hud.tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
hud.location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                          data = hud.location, mean)

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(hud.tpcb$tPCB)
hist(log10(hud.tpcb$tPCB)) # Better approach

# (2) Time trend plots
ggplot(hud.tpcb, aes(y = tPCB,
                     x = format(date,'%Y-%m'))) +
  geom_point(shape = 21, size = 2, fill = "#66ccff") +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotate("text", x = 20, y = 10^5.3, label = "Hudson River",
           size = 4)

# (3) Seasonality
ggplot(hud.tpcb, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
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
  annotate("text", x = 4, y = 10^5.7, label = "Hudson River",
           size = 4)

# (4) Sites
ggplot(hud.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concetration",
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
  annotate("text", x = 14, y = 400, label = "Hudson River",
           size = 4)

# Remove site -------------------------------------------------------------
## Remove site Bakers Falls. Upstream source
## North Bakers Falls = WCPCB-HUD006 and
## South Bakers Falls = WCPCB-HUD006.
hud.tpcb.2 <- subset(hud.tpcb, SiteID != c("WCPCB-HUD006"))
hud.tpcb.2 <- subset(hud.tpcb.2, SiteID != c("WCPCB-HUD010"))

# Plots w/o Bakers Falls sites ------------------------------------------------
# (1) Histograms
hist(hud.tpcb.2$tPCB)
hist(log10(hud.tpcb.2$tPCB))

# (2) Time trend plots
ggplot(hud.tpcb.2, aes(y = tPCB,
                       x = format(date,'%Y-%m'))) +
  geom_point(shape = 21, fill = "#66ccff") +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotate("text", x = 20, y = 10^5.2, label = "Hudson River",
           size = 4)

# (3) Seasonality
ggplot(hud.tpcb.2, aes(x = season, y = tPCB)) +
  xlab("") +
  scale_x_discrete(labels = c("winter", "spring", "summer", "fall")) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
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
  annotate("text", x = 1.2, y = 800, label = "Hudson River",
           size = 4)

# (4) Sites
ggplot(hud.tpcb.2, aes(x = factor(SiteID), y = tPCB)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab(expression("")) +
  theme(aspect.ratio = 5/20) +
  ylab(expression(bold(atop("Water Concetration",
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
  annotate("text", x = 3, y = 1000, label = "Hudson River",
           size = 4)

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Hudson River
# Hudson River
sitehudN1 <- "01331095" # HUDSON RIVER AT STILLWATER NY No temp!
sitehudN2 <- "01335754" # HUDSON RIVER ABOVE LOCK 1 NEAR WATERFORD NY, no temp!
sitehudN3 <- "01328770" # HUDSON RIVER AT THOMSON NY, no temp!
sitehudN4 <- "01327750" # HUDSON RIVER AT FORT EDWARD NY, no temp!
sitehudN5 <- "01359139" # HUDSON RIVER AT ALBANY NY No flow!

# Codes to retrieve data
paramflow <- "00060" # discharge, ft3/s
paramtemp <- "00010" # water temperature, C Not available

# Flow (ft3/s)
flow.1 <- readNWISdv(sitehudN1, paramflow,
                   min(hud.tpcb.2$date), max(hud.tpcb.2$date))
flow.2 <- readNWISdv(sitehudN2, paramflow,
                     min(hud.tpcb.2$date), max(hud.tpcb.2$date))
flow.3 <- readNWISdv(sitehudN3, paramflow,
                     min(hud.tpcb.2$date), max(hud.tpcb.2$date))
flow.4 <- readNWISdv(sitehudN4, paramflow,
                     min(hud.tpcb.2$date), max(hud.tpcb.2$date))
# Water temperature in Celsius
temp <- readNWISdv(sitehudN5, paramtemp,
                   min(hud.tpcb.2$date), max(hud.tpcb.2$date))

# Add USGS data to hud.tpcb.2, matching dates
hud.tpcb.2$flow.1 <- 0.03*flow.1$X_00060_00003[match(hud.tpcb.2$date,
                                                     flow.1$Date)]
hud.tpcb.2$flow.2 <- 0.03*flow.2$X_00060_00003[match(hud.tpcb.2$date,
                                                     flow.2$Date)]
hud.tpcb.2$flow.3 <- 0.03*flow.3$X_00060_00003[match(hud.tpcb.2$date,
                                                     flow.3$Date)]
hud.tpcb.2$flow.4 <- 0.03*flow.4$X_00060_00003[match(hud.tpcb.2$date,
                                                     flow.4$Date)]
hud.tpcb.2$temp <- 273 + temp$X_00010_00003[match(hud.tpcb.2$date,
                                                  temp$Date)]
# Remove samples with temp = NA
hud.tpcb.3 <- na.omit(hud.tpcb.2)

# Time trend plots
ggplot(hud.tpcb.3, aes(y = tPCB,
                       x = format(date,'%Y-%m'))) +
  geom_point(shape = 21, fill = "#66ccff") +
  xlab("") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotate("text", x = 10, y = 10^5.2, label = "Hudson River",
           size = 4)

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- hud.tpcb.3$tPCB
time <- hud.tpcb.3$time
site <- hud.tpcb.3$site.code
season <- hud.tpcb.3$season
flow <- hud.tpcb.3$flow.3
tem <- hud.tpcb.3$temp
# tPCB vs. time + season + flow + temp + site
lmem.hud.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + flow + tem + (1|site),
                      REML = FALSE,
                      control = lmerControl(check.nobs.vs.nlev = "ignore",
                                            check.nobs.vs.rankZ = "ignore",
                                            check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.hud.tpcb)
# Look at residuals
{
  res.hud.tpcb <- resid(lmem.hud.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.hud.tpcb, main = "log10(C)")
  qqnorm(res.hud.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.hud.tpcb)
}
# Shapiro test
shapiro.test(res.hud.tpcb)
# Random effect site Std Dev
RandonEffectSiteStdDev <- as.data.frame(VarCorr(lmem.hud.tpcb))[1,'sdcor']
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.hud.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.hud.tpcb))[1, 'R2c']

# Extract coefficient values
time.coeff <- summary(lmem.hud.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lmem.hud.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -ln(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.hud.tpcb <- as.data.frame(fitted(lmem.hud.tpcb))
# Add column name
colnames(fit.lme.values.hud.tpcb) <- c("predicted")
# Add predicted values to data.frame
hud.tpcb.3$predicted <- 10^(fit.lme.values.hud.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(hud.tpcb.3, aes(x = tPCB, y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "#66ccff") +
  scale_y_log10(limits = c(100, 10^6), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(100, 10^6), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", linewidth = 1.3) +
  geom_abline(intercept = 0.3, slope = 1, col = "blue", linewidth = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = -0.3, slope = 1, col = "blue", linewidth = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 900, y = 10^5.8,
           label = expression("Hudson River (R"^2*"= 0.80)"),
           size = 4, fontface = 2)

# Plot residuals vs. predictions
{
  plot(log10(hud.tpcb.3$predicted), res.hud.tpcb,
       points(log10(hud.tpcb.3$predicted), res.hud.tpcb, pch = 16, 
              col = "#66ccff"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  abline(0, 0)
  abline(h = c(-1, 1), col = "grey")
  abline(v = seq(2, 4.5, 0.5), col = "grey")
  }

# Plot time series with lme predictions
# Create a data frame to storage data
{
  time.serie.tpcb <- as.data.frame(matrix(nrow = length(hud.tpcb.3[,1]),
                                          ncol = 3))
  # Add name to columns
  colnames(time.serie.tpcb) <- c('date', 'tPCB', 'lmetPCB')
  # Add data
  time.serie.tpcb$date <- hud.tpcb.3$date
  time.serie.tpcb$tPCB <- hud.tpcb.3$tPCB
  time.serie.tpcb$lmetPCB <- 10^(fit.lme.values.hud.tpcb)
  # Change again the names
  colnames(time.serie.tpcb[,3]) <- c("lmetPCB")
  # Change data.frame format to be plotted
  time.serie.tpcb.2 <- melt(time.serie.tpcb, id.vars = c("date"))
}
# Plot
ggplot(time.serie.tpcb.2, aes(x = date, y = value, group = variable)) +
  geom_point(aes(shape = variable, color = variable, size = variable,
                 fill = variable)) +
  scale_shape_manual(values = c(21, 3)) +
  scale_color_manual(values = c('black','#8856a7')) +
  scale_size_manual(values = c(2, 1, 1)) +
  scale_fill_manual(values = c("#1b98e0", '#8856a7')) +
  scale_x_date(labels = date_format("%Y-%m")) +
  scale_y_log10(limits = c(100, 10^5.5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm")) +
  annotate("text", x = as.Date("2017-12-01", format = "%Y-%m-%d"),
           y = 10^5.2, label = "Hudson River", size = 3.5)



