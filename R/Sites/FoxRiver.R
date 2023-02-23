## Water PCB concentrations data analysis
## Fox River 2005 - 2018
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

# Select Fox River data ---------------------------------------------------
fox.0 <- wdc[str_detect(wdc$LocationName, 'Fox River'),]
# Lake Winnebago is a background site.
# Data preparation --------------------------------------------------------
{
  # Remove samples (rows) with total PCBs  = 0
  fox.1 <- fox.0[!(rowSums(fox.0[, c(14:117)], na.rm = TRUE)==0),]
  # Calculate total PCB
  tpcb.fox <- rowSums(fox.1[, c(14:117)], na.rm = T)
  # Change date format
  fox.1$SampleDate <- as.Date(fox.1$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(fox.1$SampleDate) - min(as.Date(fox.1$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- fox.1$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(fox.1$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  fox.tpcb <- cbind(factor(fox.1$SiteID), fox.1$SampleDate,
                    fox.1$Latitude, fox.1$Longitude, as.matrix(tpcb.fox),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(fox.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
fox.location <- fox.tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
fox.location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                            data = fox.location, mean)

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(fox.tpcb$tPCB)
hist(log10(fox.tpcb$tPCB)) # Better approach

# (2) Time trend plots
ggplot(fox.tpcb, aes(y = tPCB,
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
  annotate("text", x = 5.8, y = 10^5, label = "Fox River",
           size = 4)
  
# (3) Seasonality
ggplot(fox.tpcb, aes(x = season, y = tPCB)) +
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
  annotate("text", x = 1, y = 10^5, label = "Fox River",
           size = 4)

# (4) Sites
ggplot(fox.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 7.8, y = 10^5, label = "Fox River",
           size = 4)

# Remove site -------------------------------------------------------------
# Remove site Lake Winnebago (background site)
fox.tpcb.2 <- subset(fox.tpcb, SiteID != c("WCPCB-FOX001"))

# Plots w/o Lake Winnebago ------------------------------------------------
# (1) Histograms
hist(fox.tpcb.2$tPCB)
hist(log10(fox.tpcb.2$tPCB))

# (2) Time trend plots
ggplot(fox.tpcb.2, aes(y = tPCB,
                     x = format(date,'%Y'))) +
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
  annotate("text", x = 5.8, y = 10^5, label = "Fox River",
           size = 4)

# (3) Seasonality
ggplot(fox.tpcb.2, aes(x = season, y = tPCB)) +
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
  annotate("text", x = 1, y = 10^5, label = "Fox River",
           size = 4)

# (4) Sites
ggplot(fox.tpcb.2, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 6.8, y = 10^5, label = "Fox River",
           size = 4)

# Include USGS flow data --------------------------------------------------
{
  # Include flow data from USGS station Fox River
  sitefoxN1 <- "04084445" # flow @ OX RIVER AT APPLETON, WI
  sitefoxN2 <- "040851385" # water temperature @ FOX RIVER AT OIL TANK DEPOT AT GREEN BAY, WI
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow <- readNWISdv(sitefoxN1, paramflow,
                     min(fox.tpcb.2$date), max(fox.tpcb.2$date))
  temp <- readNWISdv(sitefoxN2, paramtemp,
                     min(fox.tpcb.2$date), max(fox.tpcb.2$date))
  # Add USGS data to fox.tpcb.2, matching dates, conversion to m3/s
  fox.tpcb.2$flow <- 0.03*flow$X_.Primary.Stream.Flow._00060_00003[match(fox.tpcb.2$date,
                                                                         flow$Date)]
  fox.tpcb.2$temp <- 273.15 + temp$X_00010_00003[match(fox.tpcb.2$date,
                                                       temp$Date)]
  # Remove samples with temp = NA
  fox.tpcb.2 <- na.omit(fox.tpcb.2)
}

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- fox.tpcb.2$tPCB
time <- fox.tpcb.2$time
site <- fox.tpcb.2$site.code
season <- fox.tpcb.2$season
flow <- fox.tpcb.2$flow
tem <- fox.tpcb.2$temp
# tPCB vs. time + season + flow + temp + site
lmem.fox.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + flow + tem + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.fox.tpcb)
# Look at residuals
{
  res.fox.tpcb <- resid(lmem.fox.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.fox.tpcb, main = "log10(C)")
  qqnorm(res.fox.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.fox.tpcb)
}
# Shapiro test
shapiro.test(res.fox.tpcb)
# Random effect site Std Dev
RandonEffectSiteStdDev <- as.data.frame(VarCorr(lmem.fox.tpcb))[1,'sdcor']
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lmem.fox.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lmem.fox.tpcb))[1, 'R2c']

# Extract coefficient values
time.coeff <- summary(lmem.fox.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lmem.fox.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -ln(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.fox.tpcb <- as.data.frame(fitted(lmem.fox.tpcb))
# Add column name
colnames(fit.lme.values.fox.tpcb) <- c("predicted")
# Add predicted values to data.frame
fox.tpcb.2$predicted <- 10^(fit.lme.values.fox.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(fox.tpcb.2, aes(x = tPCB, y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "#66ccff") +
  scale_y_log10(limits = c(10, 10^5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10, 10^5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", linewidth = 1.3) +
  geom_abline(intercept = 0.3, slope = 1, col = "blue", linewidth = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = -0.3, slope = 1, col = "blue", linewidth = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 100, y = 10^4.9,
           label = expression("Fox River (R"^2*"= 0.78)"),
           size = 4, fontface = 2)

# Plot residuals vs. predictions
{
  plot(log10(fox.tpcb.2$predicted), res.fox.tpcb,
       points(log10(fox.tpcb.2$predicted), res.fox.tpcb, pch = 16, 
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
fox.tpcb.2$factor2 <- fox.tpcb.2$tPCB/fox.tpcb.2$predicted
factor2.tpcb <- nrow(fox.tpcb.2[fox.tpcb.2$factor2 > 0.5 & fox.tpcb.2$factor2 < 2,
                              ])/length(fox.tpcb.2[,1])*100

# Plot time series with lme predictions
# Create a data frame to storage data
{
  time.serie.tpcb <- as.data.frame(matrix(nrow = length(fox.tpcb.2[,1]),
                                          ncol = 3))
  # Add name to columns
  colnames(time.serie.tpcb) <- c('date', 'tPCB', 'lmetPCB')
  # Add data
  time.serie.tpcb$date <- fox.tpcb.2$date
  time.serie.tpcb$tPCB <- fox.tpcb.2$tPCB
  time.serie.tpcb$lmetPCB <- 10^(fit.lme.values.fox.tpcb)
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
  #scale_y_log10(limits = c(10, 10000)) +
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
  annotate("text", x = as.Date("2018-06-01", format = "%Y-%m-%d"),
           y = 10^3.8, label = "Fox River", size = 3.5)

# Individual PCB Analysis -------------------------------------------------
# Use fox.1 (no 0s samples)
# Prepare data.frame
{
  fox.pcb <- subset(fox.1, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  fox.pcb <- subset(fox.pcb, select = -c(A1016:A1260))
  # Log10 individual PCBs 
  fox.pcb <- log10(fox.pcb)
  # Replace -inf to NA
  fox.pcb <- do.call(data.frame,
                     lapply(fox.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  fox.pcb.1 <- fox.pcb[,
                       -which(colSums(is.na(fox.pcb))/nrow(fox.pcb) > 0.7)]
  # Add site ID
  fox.pcb.1$SiteID <- fox.1$SiteID
  # Change date format
  fox.pcb.1$SampleDate <- as.Date(fox.1$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  fox.pcb.1$time <- as.Date(fox.1$SampleDate) - min(as.Date(fox.1$SampleDate))
  # Create individual code for each site sampled
  fox.pcb.1$site.numb <- fox.1$SiteID %>% as.factor() %>% as.numeric
  # Include season
  fox.pcb.1$season <- factor(format(yq.s, "%q"), levels = 1:4,
                           labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Remove site Lake Winnebago (background site)
  fox.pcb.1 <- subset(fox.pcb.1, SiteID != c("WCPCB-FOX001"))
  # Include flow data from USGS station Fox River
  sitefoxN1 <- "04084445" # flow @ OX RIVER AT APPLETON, WI
  sitefoxN2 <- "040851385" # water temperature @ FOX RIVER AT OIL TANK DEPOT AT GREEN BAY, WI
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  paramtemp <- "00010" # water temperature, C
  # Retrieve USGS data
  flow <- readNWISdv(sitefoxN1, paramflow,
                     min(fox.pcb.1$SampleDate), max(fox.pcb.1$SampleDate))
  temp <- readNWISdv(sitefoxN2, paramtemp,
                     min(fox.pcb.1$SampleDate), max(fox.pcb.1$SampleDate))
  # Add USGS data to fox.pcb, matching dates, conversion to m3/s
  fox.pcb.1$flow <- 0.03*flow$X_.Primary.Stream.Flow._00060_00003[match(fox.pcb.1$SampleDate,
                                                                        flow$Date)]
  fox.pcb.1$temp <- 273.15 + temp$X_00010_00003[match(fox.pcb.1$SampleDate,
                                                      temp$Date)]
  # Remove samples with temperature = NA
  fox.pcb.2 <- fox.pcb.1[!is.na(fox.pcb.1$temp), ]
  # Remove metadata
  fox.pcb.3 <- subset(fox.pcb.2, select = -c(SiteID:temp))
}

# Get covariates
time <- fox.pcb.2$time
flow <- fox.pcb.2$flow
temper <- fox.pcb.2$temp
season <- fox.pcb.2$season
site <- fox.pcb.2$site.numb

# LME for individual PCBs -------------------------------------------------
# Create matrix to store results
lme.pcb <- matrix(nrow = length(fox.pcb.3[1,]), ncol = 24)

# Perform LME
for (i in 1:length(fox.pcb.3[1,])) {
    fit <- lmer(fox.pcb.3[,i] ~ 1 + time + flow + temper + season + (1|site),
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
congeners <- colnames(fox.pcb.3)
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
write.csv(lme.pcb, file = "Output/Data/csv/LmeFoxPCB.csv")

# Get predicted values for selected PCBs
# tPCB vs. time + season + flow + temp
# lme
lme.fox.pcbi <- lmer(fox.pcb.3$PCB17 ~ 1 + time + flow + temper + season +
                       (1|site), REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                  check.nobs.vs.rankZ = "ignore",
                                  check.nobs.vs.nRE="ignore"))

# See results
summary(lme.fox.pcbi)
# Look at residuals
{
  res.lme.fox.pcbi <- resid(lme.fox.pcbi) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.lme.fox.pcbi, main = expression(paste("Normal Q-Q Plot PCB 17")))
  # Add a straight diagonal line to the plot
  qqline(res.lme.fox.pcbi)
}
# Shapiro test
shapiro.test(res.lme.fox.pcbi)
# One-sample Kolmogorov-Smirnov test
ks.test(res.lme.fox.pcbi, 'pnorm')

# (1) Get predicted values pcbi
date.pcbi <- format(fox.pcb.2$SampleDate, "%Y-%m-%d")
obs <- fox.pcb.3$PCB17
fox.pcbi <- cbind(date.pcbi, obs)
# Remove NA value from observations
fox.pcbi <- na.omit(fox.pcbi)
fit.lme.values.pcbi <- as.data.frame(fitted(lme.fox.pcbi))
fox.pcbi <- cbind(fox.pcbi, fit.lme.values.pcbi)
colnames(fox.pcbi) <- c("date", "obs", 'lme')
fox.pcbi$date <- as.Date(fox.pcbi$date)
fox.pcbi$obs <- as.numeric(fox.pcbi$obs)

# Plot residuals vs. predictions
# lme
{
  plot(fox.pcbi$lme, res.lme.fox.pcbi,
       points(fox.pcbi$lme, res.lme.fox.pcbi, pch = 16, 
              col = "#66ccff"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted concentration PCB 17 (pg/L)")),
       ylab = "Residual (lme)")
  abline(0, 0)
  abline(h = seq(-2, 2, 1), col = "grey")
  abline(v = seq(0.5, 2, 0.5), col = "grey")
}

# Modeling plots
# Change data.frame format to be plotted
fox.pcbi.2 <- melt(fox.pcbi, id.vars = c("date"))
# Plot
ggplot(fox.pcbi.2, aes(x = date, y = 10^(value), group = variable)) +
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
           y = 850, label = "Fox River", size = 3.5)

