## Water PCB concentrations data analysis per site
# Portland Harbor

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

# Select Portland Harbor data ---------------------------------------------------
por.0 <- wdc[str_detect(wdc$LocationName, 'Portland Harbor'),]

# Data preparation --------------------------------------------------------
{
  # Remove samples (rows) with total PCBs  = 0
  por.1 <- por.0[!(rowSums(por.0[, c(14:117)], na.rm = TRUE)==0),]
  # Calculate total PCB
  tpcb.por <- rowSums(por.1[, c(14:117)], na.rm = T)
  # Change date format
  por.1$SampleDate <- as.Date(por.1$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(por.1$SampleDate) - min(as.Date(por.1$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- por.1$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(por.1$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  por.tpcb <- cbind(factor(por.1$SiteID), por.1$SampleDate,
                    por.1$Latitude, por.1$Longitude, as.matrix(tpcb.por),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(por.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
por.location <- por.tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
por.location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                          data = por.location, mean)

# General plots -------------------------------------------------------------------
# (1) Histograms
# (1.1) tPCB
hist(por.tpcb$tPCB)
hist(log10(por.tpcb$tPCB))

# (2) Time trend plots
ggplot(por.tpcb, aes(y = tPCB,
                     x = format(date,'%Y'))) +
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
  annotate("text", x = 5.8, y = 10^3.7, label = "Portland Harbor",
           size = 3)

# (3) Seasonality
ggplot(por.tpcb, aes(x = season, y = tPCB)) +
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
  annotate("text", x = 1, y = 10^3.7, label = "Portland Harbor",
           size = 3)

# (4) Sites
ggplot(por.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 3, y = 10^3.7, label = "Portland Harbor",
           size = 3)

# Remove site -------------------------------------------------------------
# ?

# Include USGS flow data --------------------------------------------------
{
  # Include flow data from USGS station Portland Harbor
  sitePorN1 <- "14211720" # WILLAMETTE RIVER AT PORTLAND, OR
  sitePorN2 <- "14211820" # COLUMBIA SLOUGH AT PORTLAND, OR No!
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  # Flow (ft3/s)
  flow.1 <- readNWISdv(sitePorN1, paramflow,
                       min(por.tpcb$date), max(por.tpcb$date))
  flow.2 <- readNWISdv(sitePorN2, paramflow,
                       min(por.tpcb$date), max(por.tpcb$date))
  temp.1 <- readNWISdv(sitePorN1, paramtemp,
                       min(por.tpcb$date), max(por.tpcb$date))
  # Add USGS data to fox.tpcb, matching dates
  por.tpcb$flow.1 <- flow.1$X_00060_00003[match(por.tpcb$date, flow.1$Date)]
  por.tpcb$flow.2 <- flow.2$X_00060_00003[match(por.tpcb$date, flow.2$Date)]
  por.tpcb$temp.1 <- temp.1$X_00010_00003[match(por.tpcb$date, temp.1$Date)]
  # Remove samples with temp = NA
  por.tpcb.2 <- na.omit(por.tpcb)
}

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- por.tpcb.2$tPCB
time <- por.tpcb.2$time
site <- por.tpcb.2$site.code
season <- por.tpcb.2$season
flow <- por.tpcb.2$flow
tem <- por.tpcb.2$temp
# tPCB vs. time + season + flow + temp + site
lme.por.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + flow + tem + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.por.tpcb)
# Look at residuals
{
  res.por.tpcb <- resid(lme.por.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.por.tpcb, main = "log10(C)")
  qqnorm(res.por.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.por.tpcb)
}
# Shapiro test
shapiro.test(res.por.tpcb)
# Random effect site Std Dev
RandonEffectSiteStdDev <- as.data.frame(VarCorr(lme.por.tpcb))[1,'sdcor']
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lme.por.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lme.por.tpcb))[1, 'R2c']
# Extract coefficient values
time.coeff <- summary(lme.por.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lme.por.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -ln(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.por.tpcb <- as.data.frame(fitted(lme.por.tpcb))
# Add column name
colnames(fit.lme.values.por.tpcb) <- c("predicted")
# Add predicted values to data.frame
por.tpcb.2$predicted <- 10^(fit.lme.values.por.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(por.tpcb.2, aes(x = tPCB, y = predicted)) +
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
           label = expression(atop("Portland Harbor (R"^2*"= X)",
                                   paste("t"[1/2]*" = X ± Y (yr)"))),
           size = 3, fontface = 2)

# Plot residuals vs. predictions
{
  plot(log10(por.tpcb.2$predicted), res.por.tpcb,
       points(log10(por.tpcb.2$predicted), res.por.tpcb, pch = 16, 
              col = "#66ccff"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  abline(0, 0)
  abline(h = c(-1, 1), col = "grey")
  abline(v = seq(2, 3.5, 0.5), col = "grey")
  }

# Estimate a factor of 2 between observations and predictions
por.tpcb.2$factor2 <- por.tpcb.2$tPCB/por.tpcb.2$predicted
factor2.tpcb <- nrow(por.tpcb.2[por.tpcb.2$factor2 > 0.5 & por.tpcb.2$factor2 < 2,
                                ])/length(por.tpcb.2[,1])*100

# Plot time series with lme predictions
# Create a data frame to storage data
{
  time.serie.tpcb <- as.data.frame(matrix(nrow = length(por.tpcb.2[,1]),
                                          ncol = 3))
  # Add name to columns
  colnames(time.serie.tpcb) <- c('date', 'tPCB', 'lmetPCB')
  # Add data
  time.serie.tpcb$date <- por.tpcb.2$date
  time.serie.tpcb$tPCB <- por.tpcb.2$tPCB
  time.serie.tpcb$lmetPCB <- 10^(fit.lme.values.por.tpcb)
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
  scale_size_manual(values = c(2, 1)) +
  scale_fill_manual(values = c("#1b98e0", '#8856a7')) +
  scale_x_date(labels = date_format("%Y-%m")) +
  scale_y_log10(limits = c(10, 10^6), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("") +
  theme_bw() +
  theme(aspect.ratio = 5/15) +
  ylab(expression(bold(atop("Water Concentration",
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
  annotate("text", x = as.Date("2018-06-01", format = "%Y-%m-%d"),
           y = 10^3.8, label = "Portland Harbor", size = 3.5)

# Individual PCB Analysis -------------------------------------------------
# Use fox.1 (no 0s samples)
# Prepare data.frame
{
  por.pcb <- subset(por.1, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  por.pcb <- subset(por.pcb, select = -c(A1016:A1260))
  # Log10 individual PCBs 
  por.pcb <- log10(por.pcb)
  # Replace -inf to NA
  por.pcb <- do.call(data.frame,
                     lapply(por.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  por.pcb.1 <- por.pcb[,
                       -which(colSums(is.na(por.pcb))/nrow(por.pcb) > 0.7)]
  # Add site ID
  por.pcb.1$SiteID <- por.1$SiteID
  # Change date format
  por.pcb.1$SampleDate <- as.Date(por.1$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  por.pcb.1$time <- as.Date(por.1$SampleDate) - min(as.Date(por.1$SampleDate))
  # Create individual code for each site sampled
  por.pcb.1$site.numb <- por.1$SiteID %>% as.factor() %>% as.numeric
  # Include season
  por.pcb.1$season <- factor(format(yq.s, "%q"), levels = 1:4,
                             labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Include flow data from USGS station Portland Harbor
  sitePorN1 <- "14211720" # WILLAMETTE RIVER AT PORTLAND, OR
  sitePorN2 <- "14211820" # COLUMBIA SLOUGH AT PORTLAND, OR No!
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C Not available
  # Flow (ft3/s)
  flow.1 <- readNWISdv(sitePorN1, paramflow,
                       min(por.pcb.1$date), max(por.pcb.1$date))
  flow.2 <- readNWISdv(sitePorN2, paramflow,
                       min(por.pcb.1$date), max(por.pcb.1$date))
  temp.1 <- readNWISdv(sitePorN1, paramtemp,
                       min(por.pcb.1$date), max(por.pcb.1$date))
  # Add USGS data to fox.tpcb, matching dates, conversion to m3/s
  por.pcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(por.pcb.1$date, flow.1$Date)]
  por.pcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(por.pcb.1$date, flow.2$Date)]
  por.pcb.1$temp.1 <- 273.15 + temp.1$X_00010_00003[match(por.pcb.1$date, temp.1$Date)]
  # Remove samples with temperature = NA
  por.pcb.2 <- por.pcb.1[!is.na(por.pcb.1$temp), ]
  # Remove metadata
  por.pcb.3 <- subset(por.pcb.2, select = -c(SiteID:temp))
}

# LME for individual PCBs -------------------------------------------------
# Get covariates
time <- por.pcb.2$time
flow <- por.pcb.2$flow
temper <- por.pcb.2$temp
season <- por.pcb.2$season
site <- por.pcb.2$site.numb

# Create matrix to store results
lme.pcb <- matrix(nrow = length(por.pcb.3[1,]), ncol = 24)

# Perform LME
for (i in 1:length(por.pcb.3[1,])) {
  fit <- lmer(por.pcb.3[,i] ~ 1 + time + flow + temper + season + (1|site),
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
congeners <- colnames(por.pcb.3)
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
write.csv(lme.pcb, file = "Output/Data/csv/LmePorPCB.csv")

# Generate predictions
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(por.pcb.3[,1]),
                      ncol = length(por.pcb.3[1,]))

for (i in 1:length(por.pcb.3[1,])) {
  fit <- lmer(por.pcb.3[,i] ~ 1 + time + flow + temper + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  lme.fit.pcb[,i] <- fitted(fit)
}

# Estimate a factor of 2 between observations and predictions
factor2 <- 10^(por.pcb.3)/10^(lme.fit.pcb)
factor2.pcb <- sum(factor2 > 0.5 & factor2 < 2,
                   na.rm = TRUE)/(sum(!is.na(factor2)))*100

# Selected individual PCB regression --------------------------------------
# lme
lme.por.pcbi <- lmer(por.pcb.3$PCB17 ~ 1 + time + flow + temper + season +
                       (1|site), REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"),
                     na.action = na.exclude)

# See results
summary(lme.por.pcbi)
# Look at residuals
{
  res.lme.por.pcbi <- resid(lme.por.pcbi) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.lme.por.pcbi, main = expression(paste("Normal Q-Q Plot PCB 17")))
  # Add a straight diagonal line to the plot
  qqline(res.lme.por.pcbi)
}
# Shapiro test
shapiro.test(res.lme.por.pcbi)
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lme.por.pcbi))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lme.por.pcbi))[1, 'R2c']
# Extract coefficient values
time.coeff <- summary(lme.por.pcbi)$coef[2, "Estimate"]
time.coeff.ste <- summary(lme.por.pcbi)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -ln(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# (1) Get predicted values pcbi
date.pcbi <- format(por.pcb.2$SampleDate, "%Y-%m-%d")
obs <- por.pcb.3$PCB17
fox.pcbi <- cbind(date.pcbi, obs)
fit.lme.values.pcbi <- as.data.frame(fitted(lme.por.pcbi))
por.pcbi <- cbind(por.pcbi, fit.lme.values.pcbi)
colnames(por.pcbi) <- c("date", "obs", 'lme')
por.pcbi$date <- as.Date(por.pcbi$date)
por.pcbi$obs <- as.numeric(por.pcbi$obs)

# Plot residuals vs. predictions
{
  plot(por.pcbi$lme, res.lme.por.pcbi,
       points(por.pcbi$lme, res.lme.por.pcbi, pch = 16, 
              col = "#66ccff"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted concentration PCB 17 (pg/L)")),
       ylab = "Residual (lme)")
  abline(0, 0)
  abline(h = seq(-2, 2, 1), col = "grey")
  abline(v = seq(0.5, 2, 0.5), col = "grey")
}

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.por.pcbi <- as.data.frame(fitted(lme.por.pcbi))
# Add column name
colnames(fit.lme.values.por.pcbi) <- c("predicted")
# Add predicted values to data.frame
por.pcb.3$predicted <- 10^(fit.lme.values.por.pcbi$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(por.pcb.3, aes(x = 10^(PCB17), y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "#66ccff") +
  scale_y_log10(limits = c(0.1, 1000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.1, 1000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration PCB 17 (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration PCB 17 (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", linewidth = 1.3) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 0.7, y = 10^2.8,
           label = expression(atop("Portland Harbor (R"^2*"= 0.81)",
                                   paste("t"[1/2]*" = 14 ± 3"))),
           size = 3, fontface = 2)

# Modeling plots
# Change data.frame format to be plotted
por.pcbi.2 <- melt(por.pcbi, id.vars = c("date"))
# Plot
ggplot(por.pcbi.2, aes(x = date, y = 10^(value), group = variable)) +
  geom_point(aes(shape = variable, color = variable, size = variable,
                 fill = variable)) +
  scale_shape_manual(values = c(21, 3)) +
  scale_color_manual(values = c('black','#8856a7')) +
  scale_size_manual(values = c(2, 1)) +
  scale_fill_manual(values = c("#1b98e0", '#8856a7')) +
  scale_x_date(labels = date_format("%Y-%m")) +
  scale_y_log10(limits = c(1, 1000)) +
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
  annotate("text", x = as.Date("2018-06-01", format = "%Y-%m-%d"),
           y = 850, label = "Portland Harbor", size = 3.5)

