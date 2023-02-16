## Water PCB concentrations data analysis per site
## Kalamazoo River
## Aroclors 1242, 1254 and 1260

# Install packages
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

# Select Kalamazoo data ---------------------------------------------------
kal.0 <- wdc[str_detect(wdc$LocationName, 'Kalamazoo River'),]
# Superfund site from Morrow Dam (Kalamazoo River) to Lake Michigan
# and 30 miles of Portage Creek (south), Cork St and Portage Creek Cork St sites
# Dredging occurred at Plainwell Dam site.

# Data preparation --------------------------------------------------------
{
  # Remove samples (rows) with total PCBs  = 0
  kal.1 <- kal.0[!(rowSums(kal.0[, c(14:117)], na.rm = TRUE)==0),]
  # Calculate total PCB
  tpcb.kal <- rowSums(kal.1[, c(14:117)], na.rm = T)
  # Change date format
  kal.1$SampleDate <- as.Date(kal.1$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(kal.1$SampleDate) - min(as.Date(kal.1$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- kal.1$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(kal.1$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  kal.tpcb <- cbind(factor(kal.1$SiteID), kal.1$SampleDate,
                    kal.1$Latitude, kal.1$Longitude, as.matrix(tpcb.kal),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(kal.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
kal.location <- kal.tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
kal.location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                          data = kal.location, mean)

# General plots -------------------------------------------------------------------
# (1) Histograms
# tPCB
hist(kal.tpcb$tPCB)
hist(log10(kal.tpcb$tPCB))

# (2) Time trend plots
ggplot(kal.tpcb, aes(y = tPCB,
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
  annotate("text", x = 2, y = 100, label = "Kalamazoo River",
           size = 4)

# (3) Seasonality
ggplot(kal.tpcb, aes(x = season, y = tPCB)) +
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
  annotate("text", x = 1.2, y = 10^6.5, label = "Kalamazoo River",
           size = 4)

# (4) Sites
ggplot(kal.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 5, y = 10^6.5, label = "Kalamazoo River",
           size = 4)

# Remove site -------------------------------------------------------------
# Remove site PlainwellDam. Dredging = WCPCB-KAL023
kal.tpcb.1 <- subset(kal.tpcb, SiteID != c("WCPCB-KAL023"))

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
hist(kal.tpcb.2$tPCB)
hist(log10(kal.tpcb.2$tPCB))

# (2) Time trend plots
ggplot(kal.tpcb.1, aes(y = tPCB,
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
  annotate("text", x = 3, y = 10^5.2, label = "Kalamazoo River",
           size = 4)

# (3) Seasonality
ggplot(kal.tpcb.1, aes(x = season, y = tPCB)) +
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
  annotate("text", x = 1.1, y = 10^5.5, label = "Kalamazoo River",
           size = 4)

# (4) Sites
ggplot(kal.tpcb.1, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 5, y = 100, label = "Kalamazoo River",
           size = 4)

# Regressions -------------------------------------------------------------
# Perform Linear Mixed-Effects Model (LMEM)
tpcb <- kal.tpcb.2$tPCB
time <- kal.tpcb.2$time
site <- kal.tpcb.2$site.code
season <- kal.tpcb.2$season
flow <- kal.tpcb.2$flow.3

lmem.kal.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + flow + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.kal.tpcb)
# Look at residuals
{
  res.kal.tpcb <- resid(lmem.kal.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.kal.tpcb, main = "log10(C)")
  qqnorm(res.kal.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.kal.tpcb)
}
# Shapiro test
shapiro.test(res.kal.tpcb)
# One-sample Kolmogorov-Smirnov test
ks.test(res.kal.tpcb, 'pnorm')
# Random effect site Std Dev
RandonEffectSiteStdDev <- as.data.frame(VarCorr(lmem.kal.tpcb))[1,'sdcor']
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

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.kal.tpcb <- as.data.frame(fitted(lmem.kal.tpcb))
# Add column name
colnames(fit.lme.values.kal.tpcb) <- c("predicted")
# Add predicted values to data.frame
kal.tpcb.2$predicted <- 10^(fit.lme.values.kal.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(kal.tpcb.2, aes(x = tPCB, y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "#66ccff") +
  scale_y_log10(limits = c(10, 10^5.5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10, 10^5.5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", linewidth = 1.3) +
  geom_abline(intercept = 0.3, slope = 1, col = "blue", linewidth = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = -0.3, slope = 1, col = "blue", linewidth = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 200, y = 10^5.4,
           label = expression("Kalamazoo River (R"^2*"= 0.97)"),
           size = 4, fontface = 2)

# Plot residuals vs. predictions
{
  plot(log10(kal.tpcb.2$predicted), res.kal.tpcb,
       points(log10(kal.tpcb.2$predicted), res.kal.tpcb, pch = 16, 
              col = "#66ccff"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  abline(0, 0)
  abline(h = c(-1, 1), col = "grey")
  abline(v = seq(1, 5, 0.5), col = "grey")
  }

# Plot time series with lme predictions
# Create a data frame to storage data
{
  time.serie.tpcb <- as.data.frame(matrix(nrow = length(kal.tpcb.2[,1]),
                                          ncol = 3))
  # Add name to columns
  colnames(time.serie.tpcb) <- c('date', 'tPCB', 'lmetPCB')
  # Add data
  time.serie.tpcb$date <- kal.tpcb.2$date
  time.serie.tpcb$tPCB <- kal.tpcb.2$tPCB
  time.serie.tpcb$lmetPCB <- 10^(fit.lme.values.kal.tpcb)
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
  scale_y_log10(limits = c(10, 1000000)) +
  #xlim(as.Date("2008-11-01", format = "%Y-%m-%d"),
  #     x = as.Date("2009-02-01", format = "%Y-%m-%d")) +
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
  annotate("text", x = as.Date("2006-08-01", format = "%Y-%m-%d"),
           y = 10^5, label = "Kalamazoo River", size = 3.5)

# Individual PCB Analysis -------------------------------------------------
# Use fox.1 (no 0s samples)
# Prepare data.frame
{
  kal.pcb <- subset(kal.1, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  kal.pcb <- subset(kal.pcb, select = -c(A1016:A1260))
  # Log10 individual PCBs 
  kal.pcb <- log10(kal.pcb)
  # Replace -inf to NA
  kal.pcb <- do.call(data.frame,
                     lapply(kal.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  kal.pcb.1 <- kal.pcb[,
                       -which(colSums(is.na(kal.pcb))/nrow(kal.pcb) > 0.7)]
  # Add site ID
  kal.pcb.1$SiteID <- kal.1$SiteID
  # Change date format
  kal.pcb.1$SampleDate <- as.Date(kal.1$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  kal.pcb.1$time <- as.Date(kal.1$SampleDate) - min(as.Date(kal.1$SampleDate))
  # Create individual code for each site sampled
  kal.pcb.1$site.numb <- kal.1$SiteID %>% as.factor() %>% as.numeric
  # Include season
  kal.pcb.1$season <- factor(format(yq.s, "%q"), levels = 1:4,
                             labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Remove site PlainwellDam. Dredging = WCPCB-KAL023
  kal.pcb.1 <- subset(kal.pcb.1, SiteID != c("WCPCB-KAL023"))
  # Include flow data from USGS station Kalamazoo River
  siteKalN1 <- "04108660" # KALAMAZOO RIVER AT NEW RICHMOND, MI
  siteKalN2 <- "04106000" # KALAMAZOO RIVER AT COMSTOCK, MI
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  # Flow (ft3/s)
  flow.1 <- readNWISdv(siteKalN1, paramflow,
                       min(kal.pcb.1$SampleDate), max(kal.pcb.1$SampleDate))
  flow.2 <- readNWISdv(siteKalN2, paramflow,
                       min(kal.pcb.1$SampleDate), max(kal.pcb.1$SampleDate))
  
  kal.pcb.1$flow.1 <- 0.03*flow.1$X_00060_00003[match(kal.pcb.1$SampleDate,
                                                      flow.1$Date)]
  kal.pcb.1$flow.2 <- 0.03*flow.2$X_00060_00003[match(kal.pcb.1$SampleDate,
                                                      flow.2$Date)]
  # Create flow, flow.3
  kal.pcb.1$flow.3 <- ifelse(is.na(kal.pcb.1$flow.1) == TRUE,
                              kal.pcb.1$flow.2/0.46, kal.pcb.1$flow.1)
  # Remove samples with flow.1 = NA
  kal.pcb.2 <- na.omit(kal.pcb.1)
  # Remove metadata
  kal.pcb.3 <- subset(kal.pcb.2, select = -c(SiteID:flow.3))
}

# Something odd occurring.

# Get covariates
time <- kal.pcb.2$time
flow <- kal.pcb.2$flow.3
season <- kal.pcb.2$season
site <- kal.pcb.2$site.numb

# LME for individual PCBs -------------------------------------------------
# Create matrix to store results
lme.pcb <- matrix(nrow = length(kal.pcb.3[1,]), ncol = 24)

a <- lmer(10^(kal.pcb.3$PCB12.13) ~ 1 + time + (1|site),
     REML = FALSE,
     control = lmerControl(check.nobs.vs.nlev = "ignore",
                           check.nobs.vs.rankZ = "ignore",
                           check.nobs.vs.nRE="ignore"))


summary(a)

# Perform LME
for (i in 1:length(kal.pcb.3[1,])) {
  fit <- lmer(kal.pcb.3[,i] ~ 1 + time + flow + season + (1|site),
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
  lme.pcb[i,10] <- fixef(fit)[4] # season 1
  lme.pcb[i,11] <- summary(fit)$coef[4,"Std. Error"] # season 1 error
  lme.pcb[i,12] <- summary(fit)$coef[4,"Pr(>|t|)"] # season 1 p-value
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
congeners <- colnames(kal.pcb.3)
lme.pcb <- cbind(congeners, lme.pcb)
# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "flow", "flow.error", "flow.pv", "season2",
                       "season2.error", "season2, pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality")

# Export results
write.csv(lme.pcb, file = "Output/Data/csv/LmeKalmazooPCB.csv")

  
  
