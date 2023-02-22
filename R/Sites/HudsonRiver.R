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

# Individual PCB Analysis -------------------------------------------------
# Use hud.1 (no 0s samples)
# Prepare data.frame
{
  hud.pcb <- subset(hud.1, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  hud.pcb <- subset(hud.pcb, select = -c(A1016:A1260))
  # Log10 individual PCBs 
  hud.pcb <- log10(hud.pcb)
  # Replace -inf to NA
  hud.pcb <- do.call(data.frame,
                     lapply(hud.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  hud.pcb.1 <- hud.pcb[,
                       -which(colSums(is.na(hud.pcb))/nrow(hud.pcb) > 0.7)]
  # Add site ID
  hud.pcb.1$SiteID <- hud.1$SiteID
  # Change date format
  hud.pcb.1$SampleDate <- as.Date(hud.1$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  hud.pcb.1$time <- as.Date(hud.1$SampleDate) - min(as.Date(hud.1$SampleDate))
  # Create individual code for each site sampled
  hud.pcb.1$site.numb <- hud.1$SiteID %>% as.factor() %>% as.numeric
  # Include season
  hud.pcb.1$season <- factor(format(yq.s, "%q"), levels = 1:4,
                             labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  ## Remove site Bakers Falls. Upstream source
  ## North Bakers Falls = WCPCB-HUD006 and
  ## South Bakers Falls = WCPCB-HUD006.
  hud.pcb.1 <- subset(hud.pcb.1, SiteID != c("WCPCB-HUD006"))
  hud.pcb.1 <- subset(hud.pcb.1, SiteID != c("WCPCB-HUD010"))
  # Add USGS data to hud.tpcb.2, matching dates
  hud.pcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(hud.pcb.1$SampleDate,
                                                       flow.1$Date)]
  hud.pcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(hud.pcb.1$SampleDate,
                                                       flow.2$Date)]
  hud.pcb.1$flow.3 <- 0.03*flow.3$X_00060_00003[match(hud.pcb.1$SampleDate,
                                                       flow.3$Date)]
  hud.pcb.1$flow.4 <- 0.03*flow.4$X_00060_00003[match(hud.pcb.1$SampleDate,
                                                       flow.4$Date)]
  hud.pcb.1$temp <- 273 + temp$X_00010_00003[match(hud.pcb.1$SampleDate,
                                                    temp$Date)]
  # Remove metadata
  hud.pcb.2 <- subset(hud.pcb.1, select = -c(SiteID:temp))
}
  
# Get covariates
time <- hud.pcb.1$time
flow <- hud.pcb.1$flow.1
temper <- hud.pcb.1$temp
season <- hud.pcb.1$season
site <- hud.pcb.1$site.numb

# LME for individual PCBs -------------------------------------------------
# Create matrix to store results
lme.pcb <- matrix(nrow = length(hud.pcb.2[1,]), ncol = 24)

# Perform LME
for (i in 1:length(hud.pcb.2[1,])) {
  fit <- lmer(hud.pcb.2[,i] ~ 1 + time + flow + temper + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"))
  lme.pcb[i,1] <- fixef(fit)[1] # intercept
  lme.pcb[i,2] <- summary(fit)$coef[1,"Std. Error"] # intercept error
  lme.pcb[i,3] <- summary(fit)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.pcb[i,4] <- fixef(fit)[2] # time
  lme.pcb[i,5] <- summary(fit)$coef[2,"Std. Error"] # time error
  lme.pcb[i,6] <- summary(fit)$coef[2,"Pr(>|t|)"] # time p-value
  lme.pcb[i,7] <- fixef(fit)[3] # flow
  lme.pcb[i,8] <- summary(fit)$coef[3,"Std. Error"] # flow error
  lme.pcb[i,9] <- summary(fit)$coef[3,"Pr(>|t|)"] # flow p-value
  lme.pcb[i,10] <- fixef(fit)[4] # temperature
  lme.pcb[i,11] <- summary(fit)$coef[4,"Std. Error"] # temperature error
  lme.pcb[i,12] <- summary(fit)$coef[4,"Pr(>|t|)"] # temperature p-value
  lme.pcb[i,13] <- fixef(fit)[5] # season 2
  lme.pcb[i,14] <- summary(fit)$coef[5,"Std. Error"] # season 2 error
  lme.pcb[i,15] <- summary(fit)$coef[5,"Pr(>|t|)"] # season 2 p-value
  lme.pcb[i,16] <- fixef(fit)[6] # season 3
  lme.pcb[i,17] <- summary(fit)$coef[6,"Std. Error"] # season 3 error
  lme.pcb[i,18] <- summary(fit)$coef[6,"Pr(>|t|)"] # season 3 p-value
  lme.pcb[i,19] <- -log(2)/lme.pcb[i,4]/365 # t0.5
  lme.pcb[i,20] <- abs(-log(2)/lme.pcb[i,4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
  lme.pcb[i,21] <- as.data.frame(VarCorr(fit))[1,'sdcor']
  lme.pcb[i,22] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
  lme.pcb[i,23] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
  lme.pcb[i,24] <- shapiro.test(resid(fit))$p.value
}

# Just 3 significant figures
lme.pcb <- formatC(signif(lme.pcb, digits = 3))
# Add congener names
congeners <- colnames(hud.pcb.2)
lme.pcb <- cbind(congeners, lme.pcb)
# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "flow", "flow.error", "flow.pv", "temperature",
                       "temperature.error", "temperature.pv", "season2",
                       "season2.error", "season2, pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality")

# Export results
write.csv(lme.pcb, file = "Output/Data/csv/LmeHudPCB.csv")

# Get predicted values for selected PCBs
# tPCB vs. time + season + flow + temp
# lme
lme.hud.pcbi <- lmer(hud.pcb.2$PCB17 ~ 1 + time + flow + temper + season +
                       (1|site), REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.hud.pcbi)
# Look at residuals
{
  res.lme.hud.pcbi <- resid(lme.hud.pcbi) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.lme.hud.pcbi, main = expression(paste("Normal Q-Q Plot PCB 17")))
  # Add a straight diagonal line to the plot
  qqline(res.lme.hud.pcbi)
}
# Shapiro test
shapiro.test(res.lme.hud.pcbi)

# (1) Get predicted values pcbi
date.pcbi <- format(hud.pcb.1$SampleDate, "%Y-%m-%d")
obs <- hud.pcb.2$PCB17
hud.pcbi <- cbind(date.pcbi, obs)
# Remove NA value from observations
hud.pcbi <- na.omit(hud.pcbi)
fit.lme.values.pcbi <- as.data.frame(fitted(lme.hud.pcbi))
hud.pcbi <- cbind(hud.pcbi, fit.lme.values.pcbi)
colnames(hud.pcbi) <- c("date", "obs", 'lme')
hud.pcbi$date <- as.Date(hud.pcbi$date)
hud.pcbi$obs <- as.numeric(hud.pcbi$obs)

# Plot residuals vs. predictions
# lme
{
  plot(hud.pcbi$lme, res.lme.hud.pcbi,
       points(hud.pcbi$lme, res.lme.hud.pcbi, pch = 16, 
              col = "#66ccff"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted concentration PCB 17 (pg/L)")),
       ylab = "Residual (lme)")
  abline(0, 0)
  abline(h = seq(-2, 2, 1), col = "grey")
  abline(v = seq(0.5, 3, 0.5), col = "grey")
}

# Modeling plots
# Change data.frame format to be plotted
hud.pcbi.2 <- melt(hud.pcbi, id.vars = c("date"))
# Plot
ggplot(hud.pcbi.2, aes(x = date, y = 10^(value), group = variable)) +
  geom_point(aes(shape = variable, color = variable, size = variable,
                 fill = variable)) +
  scale_shape_manual(values = c(21, 3)) +
  scale_color_manual(values = c('black','#8856a7')) +
  scale_size_manual(values = c(2, 1)) +
  scale_fill_manual(values = c("#1b98e0", '#8856a7')) +
  scale_x_date(labels = date_format("%Y-%m")) +
  scale_y_log10(limits = c(1, 10000)) +
  xlab("") +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concetration",
                            paste("PCB 17 (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 8)) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm")) +
  annotate("text", x = as.Date("2017-01-01", format = "%Y-%m-%d"),
           y = 10000, label = "Hudson River", size = 3.5)


