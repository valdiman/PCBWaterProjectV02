## Water PCB concentrations data analysis per site
## Housatonic River
## Aroclors 1254 and 1260, no congener analysis
## GE facility map @https://semspub.epa.gov/work/01/574882.pdf

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

# Select Housatonic River data ---------------------------------------------------
hou.0 <- wdc[str_detect(wdc$LocationName, 'Housatonic River'),]

# Data preparation --------------------------------------------------------
{
  # Remove samples (rows) with total PCBs  = 0
  hou.1 <- hou.0[!(rowSums(hou.0[, c(14:117)], na.rm = TRUE)==0),]
  # Calculate total PCB
  tpcb.hou <- rowSums(hou.1[, c(14:117)], na.rm = T)
  # Change date format
  hou.1$SampleDate <- as.Date(hou.1$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(hou.1$SampleDate) - min(as.Date(hou.1$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- hou.1$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(hou.1$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  hou.tpcb <- cbind(factor(hou.1$SiteID), hou.1$SampleDate,
                    hou.1$Latitude, hou.1$Longitude, as.matrix(tpcb.hou),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(hou.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
hou.location <- hou.tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
hou.location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                          data = hou.location, mean)

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(hou.tpcb$tPCB)
hist(log10(hou.tpcb$tPCB))

# (2) Time trend plots
ggplot(hou.tpcb, aes(y = tPCB,
                     x = format(date,'%Y'))) +
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
  annotate("text", x = 23, y = 10^5.6, label = "Housotonic River",
           size = 3)

# (3) Seasonality
ggplot(hou.tpcb, aes(x = season, y = tPCB)) +
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
  annotate("text", x = 4, y = 10^5.6, label = "Housotonic River",
           size = 3)

# (4) Sites
ggplot(hou.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 32, y = 10^5.6, label = "Housotonic River",
           size = 3)

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Housatonic River
{
  siteHouN1 <- "01197000" # EAST BRANCH HOUSATONIC RIVER AT COLTSVILLE, MA
  siteHouN2 <- "01197500" # HOUSATONIC RIVER NEAR GREAT BARRINGTON, MA
  # Codes to retrieve data
  paramflow <- "00060" # discharge, ft3/s
  # Retrieve USGS data
  flow.1 <- readNWISdv(siteHouN1, paramflow,
                       min(hou.tpcb$date), max(hou.tpcb$date))
  flow.2 <- readNWISdv(siteHouN2, paramflow,
                       min(hou.tpcb$date), max(hou.tpcb$date))
  # Add USGS data to hou.tpcb, matching dates (m3/s, 0.03 conversion factor)
  hou.tpcb$flow.1 <- 0.03*flow.1$X_00060_00003[match(hou.tpcb$date,
                                                     flow.1$Date)]
  hou.tpcb$flow.2 <- 0.03*flow.2$X_00060_00003[match(hou.tpcb$date,
                                                     flow.2$Date)]
}

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- hou.tpcb$tPCB
time <- hou.tpcb$time
flow <- hou.tpcb$flow.1
site <- hou.tpcb$site.code
season <- hou.tpcb$season

# tPCB vs. time + flow + season + site
lme.hou.tpcb <- lmer(log10(tpcb) ~ 1 + time + flow + season + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lme.hou.tpcb)

# Look at residuals
{
  res.hou.tpcb <- resid(lme.hou.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.hou.tpcb, main = "log10(C)")
  qqnorm(res.hou.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.hou.tpcb)
}
# Shapiro test
shapiro.test(res.hou.tpcb)
# Lme does not provide a good model for all the data.
# If data are separated in ranges, lme models works.
# (i) > 982
hou.tpcb.1 <- subset(hou.tpcb, tPCB > 982)
tpcb <- hou.tpcb.1$tPCB
time <- hou.tpcb.1$time
flow <- hou.tpcb.1$flow.1
site <- hou.tpcb.1$site.code
season <- hou.tpcb.1$season

# tPCB vs. time + flow + season + site
lme.hou.tpcb.1 <- lmer(log10(tpcb) ~ 1 + time + flow + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.hou.tpcb.1)

{
  res.hou.tpcb.1 <- resid(lme.hou.tpcb.1) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.hou.tpcb.1, main = "log10(C)")
  qqnorm(res.hou.tpcb.1,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.hou.tpcb.1)
}
# Create matrix to store results
{
  lme.tpcb.1 <- matrix(nrow = 1, ncol = 24)
  lme.tpcb.1[1] <- fixef(lme.hou.tpcb.1)[1] # intercept
  lme.tpcb.1[2] <- summary(lme.hou.tpcb.1)$coef[1,"Std. Error"] # intercept error
  lme.tpcb.1[3] <- summary(lme.hou.tpcb.1)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb.1[4] <- fixef(lme.hou.tpcb.1)[2] # time
  lme.tpcb.1[5] <- summary(lme.hou.tpcb.1)$coef[2,"Std. Error"] # time error
  lme.tpcb.1[6] <- summary(lme.hou.tpcb.1)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb.1[7] <- fixef(lme.hou.tpcb.1)[3] # flow
  lme.tpcb.1[8] <- summary(lme.hou.tpcb.1)$coef[3,"Std. Error"] # flow error
  lme.tpcb.1[9] <- summary(lme.hou.tpcb.1)$coef[3,"Pr(>|t|)"] # flow p-value
  lme.tpcb.1[10] <- fixef(lme.hou.tpcb.1)[4] # season 1
  lme.tpcb.1[11] <- summary(lme.hou.tpcb.1)$coef[4,"Std. Error"] # season 1 error
  lme.tpcb.1[12] <- summary(lme.hou.tpcb.1)$coef[4,"Pr(>|t|)"] # season 1 p-value
  lme.tpcb.1[13] <- fixef(lme.hou.tpcb.1)[5] # season 2
  lme.tpcb.1[14] <- summary(lme.hou.tpcb.1)$coef[5,"Std. Error"] # season 2 error
  lme.tpcb.1[15] <- summary(lme.hou.tpcb.1)$coef[5,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb.1[16] <- fixef(lme.hou.tpcb.1)[6] # season 3
  lme.tpcb.1[17] <- summary(lme.hou.tpcb.1)$coef[6,"Std. Error"] # season 3 error
  lme.tpcb.1[18] <- summary(lme.hou.tpcb.1)$coef[6,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb.1[19] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb.1[20] <- abs(-log(2)/lme.tpcb.1[4]/365)*lme.tpcb.1[5]/abs(lme.tpcb.1[4]) # t0.5 error
  lme.tpcb.1[21] <- as.data.frame(VarCorr(lme.hou.tpcb.1))[1,'sdcor']
  lme.tpcb.1[22] <- as.data.frame(r.squaredGLMM(lme.hou.tpcb.1))[1, 'R2m']
  lme.tpcb.1[23] <- as.data.frame(r.squaredGLMM(lme.hou.tpcb.1))[1, 'R2c']
  lme.tpcb.1[24] <- shapiro.test(resid(lme.hou.tpcb.1))$p.value
}

# Just 3 significant figures
lme.tpcb.1 <- formatC(signif(lme.tpcb.1, digits = 3))
# Add column names
colnames(lme.tpcb.1) <- c("Intercept", "Intercept.error",
                        "Intercept.pv", "time", "time.error", "time.pv",
                        "flow", "flow.error", "flow.pv", "season1",
                        "season1.error", "season1.pv", "season2",
                        "season2.error", "season2, pv", "season3",
                        "season3.error", "season3.pv", "t05", "t05.error",
                        "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality")

# Export results
write.csv(lme.tpcb, file = "Output/Data/Sites/csv/HousatonicLmetPCBV01.csv")

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.hou.tpcb.1 <- as.data.frame(fitted(lme.hou.tpcb.1))
# Add column name
colnames(fit.lme.values.hou.tpcb.1) <- c("predicted")
# Add predicted values to data.frame
hou.tpcb.1$predicted <- 10^(fit.lme.values.hou.tpcb.1$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(hou.tpcb.1, aes(x = tPCB, y = predicted)) +
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
           label = expression("Housatonic River (R"^2*"= 0.75)"),
           size = 3, fontface = 2)

# Plot residuals vs. predictions
{
  plot(log10(hou.tpcb.1$predicted), res.hou.tpcb.1,
       points(log10(hou.tpcb.1$predicted), res.hou.tpcb.1, pch = 16, 
              col = "#66ccff"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  abline(0, 0)
  abline(h = c(-1, 1), col = "grey")
  abline(v = seq(2, 5, 0.5), col = "grey")
  }

# Estimate a factor of 2 between observations and predictions
hou.tpcb.1$factor2 <- hou.tpcb.1$tPCB/hou.tpcb.1$predicted
factor2.tpcb.1 <- nrow(hou.tpcb.1[hou.tpcb.1$factor2 > 0.5 & hou.tpcb.1$factor2 < 2,
                                  ])/length(hou.tpcb.1[,1])*100

# (ii) 750 - 1300
hou.tpcb.2 <- subset(hou.tpcb, tPCB < 1300 & tPCB > 750)
tpcb <- hou.tpcb.2$tPCB
time <- hou.tpcb.2$time
flow <- hou.tpcb.2$flow.1
site <- hou.tpcb.2$site.code
season <- hou.tpcb.2$season

# tPCB vs. time + flow + season + site
lme.hou.tpcb.2 <- lmer(log10(tpcb) ~ 1 + time + flow + season + (1|site),
                       REML = FALSE,
                       control = lmerControl(check.nobs.vs.nlev = "ignore",
                                             check.nobs.vs.rankZ = "ignore",
                                             check.nobs.vs.nRE="ignore"))

# See results
summary(lme.hou.tpcb.2)

{
  res.hou.tpcb.2 <- resid(lme.hou.tpcb.2) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.hou.tpcb.2, main = "log10(C)")
  qqnorm(res.hou.tpcb.2,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.hou.tpcb.2)
}

# Create matrix to store results
{
  lme.tpcb.2 <- matrix(nrow = 1, ncol = 24)
  lme.tpcb.2[1] <- fixef(lme.hou.tpcb.2)[1] # intercept
  lme.tpcb.2[2] <- summary(lme.hou.tpcb.2)$coef[1,"Std. Error"] # intercept error
  lme.tpcb.2[3] <- summary(lme.hou.tpcb.2)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb.2[4] <- fixef(lme.hou.tpcb.2)[2] # time
  lme.tpcb.2[5] <- summary(lme.hou.tpcb.2)$coef[2,"Std. Error"] # time error
  lme.tpcb.2[6] <- summary(lme.hou.tpcb.2)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb.2[7] <- fixef(lme.hou.tpcb.2)[3] # flow
  lme.tpcb.2[8] <- summary(lme.hou.tpcb.2)$coef[3,"Std. Error"] # flow error
  lme.tpcb.2[9] <- summary(lme.hou.tpcb.2)$coef[3,"Pr(>|t|)"] # flow p-value
  lme.tpcb.2[10] <- fixef(lme.hou.tpcb.2)[4] # season 1
  lme.tpcb.2[11] <- summary(lme.hou.tpcb.2)$coef[4,"Std. Error"] # season 1 error
  lme.tpcb.2[12] <- summary(lme.hou.tpcb.2)$coef[4,"Pr(>|t|)"] # season 1 p-value
  lme.tpcb.2[13] <- fixef(lme.hou.tpcb.2)[5] # season 2
  lme.tpcb.2[14] <- summary(lme.hou.tpcb.2)$coef[5,"Std. Error"] # season 2 error
  lme.tpcb.2[15] <- summary(lme.hou.tpcb.2)$coef[5,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb.2[16] <- fixef(lme.hou.tpcb.2)[6] # season 3
  lme.tpcb.2[17] <- summary(lme.hou.tpcb.2)$coef[6,"Std. Error"] # season 3 error
  lme.tpcb.2[18] <- summary(lme.hou.tpcb.2)$coef[6,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb.2[19] <- -log(2)/lme.tpcb.2[4]/365 # t0.5
  lme.tpcb.2[20] <- abs(-log(2)/lme.tpcb.2[4]/365)*lme.tpcb.2[5]/abs(lme.tpcb.2[4]) # t0.5 error
  lme.tpcb.2[21] <- as.data.frame(VarCorr(lme.hou.tpcb.2))[1,'sdcor']
  lme.tpcb.2[22] <- as.data.frame(r.squaredGLMM(lme.hou.tpcb.2))[1, 'R2m']
  lme.tpcb.2[23] <- as.data.frame(r.squaredGLMM(lme.hou.tpcb.2))[1, 'R2c']
  lme.tpcb.2[24] <- shapiro.test(resid(lme.hou.tpcb.2))$p.value
}

# Just 3 significant figures
lme.tpcb.2 <- formatC(signif(lme.tpcb.2, digits = 3))
# Add column names
colnames(lme.tpcb.2) <- c("Intercept", "Intercept.error",
                          "Intercept.pv", "time", "time.error", "time.pv",
                          "flow", "flow.error", "flow.pv", "season1",
                          "season1.error", "season1.pv", "season2",
                          "season2.error", "season2, pv", "season3",
                          "season3.error", "season3.pv", "t05", "t05.error",
                          "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality")

# Export results
write.csv(lme.tpcb, file = "Output/Data/Sites/csv/HousatonicLmetPCBV02.csv")

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.hou.tpcb.2 <- as.data.frame(fitted(lme.hou.tpcb.2))
# Add column name
colnames(fit.lme.values.hou.tpcb.2) <- c("predicted")
# Add predicted values to data.frame
hou.tpcb.2$predicted <- 10^(fit.lme.values.hou.tpcb.2$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(hou.tpcb.2, aes(x = tPCB, y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "#66ccff") +
  scale_y_log10(limits = c(10^2.6, 10^3.4), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10^2.5, 10^3.4), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", linewidth = 1.3) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 500, y = 10^3.3,
           label = expression(atop("Housatonic River (R"^2*"= 0.5)",
                                   paste("t"[1/2]*" = 175 ± 51 (yr)"))),
           size = 3, fontface = 2)

# Plot residuals vs. predictions
{
  plot(log10(hou.tpcb.2$predicted), res.hou.tpcb.2,
       points(log10(hou.tpcb.2$predicted), res.hou.tpcb.2, pch = 16, 
              col = "#66ccff"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  abline(0, 0)
  abline(h = c(-1, 1), col = "grey")
  abline(v = seq(2.8, 3.5, 0.02), col = "grey")
  }

# Estimate a factor of 2 between observations and predictions
hou.tpcb.2$factor2 <- hou.tpcb.2$tPCB/hou.tpcb.2$predicted
factor2.tpcb.2 <- nrow(hou.tpcb.2[hou.tpcb.2$factor2 > 0.5 & hou.tpcb.2$factor2 < 2,
                                  ])/length(hou.tpcb.2[,1])*100

# (iii) 79 - 100
hou.tpcb.3 <- subset(hou.tpcb, tPCB < 100 & tPCB > 79)
tpcb <- hou.tpcb.3$tPCB
time <- hou.tpcb.3$time
flow <- hou.tpcb.3$flow.1
site <- hou.tpcb.3$site.code
season <- hou.tpcb.3$season

# tPCB vs. time + flow + season + site
lme.hou.tpcb.3 <- lmer(log10(tpcb) ~ 1 + time + flow + season + (1|site),
                       REML = FALSE,
                       control = lmerControl(check.nobs.vs.nlev = "ignore",
                                             check.nobs.vs.rankZ = "ignore",
                                             check.nobs.vs.nRE="ignore"))

# See results
summary(lme.hou.tpcb.3)

{
  res.hou.tpcb.3 <- resid(lme.hou.tpcb.3) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.hou.tpcb.3, main = "log10(C)")
  qqnorm(res.hou.tpcb.3,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.hou.tpcb.3)
}

# Create matrix to store results
{
  lme.tpcb.3 <- matrix(nrow = 1, ncol = 24)
  lme.tpcb.3[1] <- fixef(lme.hou.tpcb.3)[1] # intercept
  lme.tpcb.3[2] <- summary(lme.hou.tpcb.3)$coef[1,"Std. Error"] # intercept error
  lme.tpcb.3[3] <- summary(lme.hou.tpcb.3)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb.3[4] <- fixef(lme.hou.tpcb.3)[2] # time
  lme.tpcb.3[5] <- summary(lme.hou.tpcb.3)$coef[2,"Std. Error"] # time error
  lme.tpcb.3[6] <- summary(lme.hou.tpcb.3)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb.3[7] <- fixef(lme.hou.tpcb.3)[3] # flow
  lme.tpcb.3[8] <- summary(lme.hou.tpcb.3)$coef[3,"Std. Error"] # flow error
  lme.tpcb.3[9] <- summary(lme.hou.tpcb.3)$coef[3,"Pr(>|t|)"] # flow p-value
  lme.tpcb.3[10] <- fixef(lme.hou.tpcb.3)[4] # season 1
  lme.tpcb.3[11] <- summary(lme.hou.tpcb.3)$coef[4,"Std. Error"] # season 1 error
  lme.tpcb.3[12] <- summary(lme.hou.tpcb.3)$coef[4,"Pr(>|t|)"] # season 1 p-value
  lme.tpcb.3[13] <- fixef(lme.hou.tpcb.3)[5] # season 2
  lme.tpcb.3[14] <- summary(lme.hou.tpcb.3)$coef[5,"Std. Error"] # season 2 error
  lme.tpcb.3[15] <- summary(lme.hou.tpcb.3)$coef[5,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb.3[16] <- fixef(lme.hou.tpcb.3)[6] # season 3
  lme.tpcb.3[17] <- summary(lme.hou.tpcb.3)$coef[6,"Std. Error"] # season 3 error
  lme.tpcb.3[18] <- summary(lme.hou.tpcb.3)$coef[6,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb.3[19] <- -log(2)/lme.tpcb.3[4]/365 # t0.5
  lme.tpcb.3[20] <- abs(-log(2)/lme.tpcb.3[4]/365)*lme.tpcb.3[5]/abs(lme.tpcb.3[4]) # t0.5 error
  lme.tpcb.3[21] <- as.data.frame(VarCorr(lme.hou.tpcb.3))[1,'sdcor']
  lme.tpcb.3[22] <- as.data.frame(r.squaredGLMM(lme.hou.tpcb.3))[1, 'R2m']
  lme.tpcb.3[23] <- as.data.frame(r.squaredGLMM(lme.hou.tpcb.3))[1, 'R2c']
  lme.tpcb.3[24] <- shapiro.test(resid(lme.hou.tpcb.3))$p.value
}

# Just 3 significant figures
lme.tpcb.3 <- formatC(signif(lme.tpcb.3, digits = 3))
# Add column names
colnames(lme.tpcb.3) <- c("Intercept", "Intercept.error",
                          "Intercept.pv", "time", "time.error", "time.pv",
                          "flow", "flow.error", "flow.pv", "season1",
                          "season1.error", "season1.pv", "season2",
                          "season2.error", "season2, pv", "season3",
                          "season3.error", "season3.pv", "t05", "t05.error",
                          "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality")

# Export results
write.csv(lme.tpcb, file = "Output/Data/Sites/csv/HousatonicLmetPCBV03.csv")

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.hou.tpcb.3 <- as.data.frame(fitted(lme.hou.tpcb.3))
# Add column name
colnames(fit.lme.values.hou.tpcb.3) <- c("predicted")
# Add predicted values to data.frame
hou.tpcb.3$predicted <- 10^(fit.lme.values.hou.tpcb.3$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(hou.tpcb.3, aes(x = tPCB, y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "#66ccff") +
  scale_y_log10(limits = c(50, 10^2.2), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(50, 10^2.2), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", linewidth = 1.3) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 65, y = 10^2.15,
           label = expression(atop("Housatonic River (R"^2*"= 0.11)",
                                   paste("t"[1/2]*" = 624 ± 230 (yr)"))),
           size = 3, fontface = 2)

# Plot residuals vs. predictions
{
  plot(log10(hou.tpcb.3$predicted), res.hou.tpcb.3,
       points(log10(hou.tpcb.3$predicted), res.hou.tpcb.3, pch = 16, 
              col = "#66ccff"),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  abline(0, 0)
  abline(h = c(-1, 1), col = "grey")
  abline(v = seq(1.92, 2, 0.01), col = "grey")
  }

# Estimate a factor of 2 between observations and predictions
hou.tpcb.3$factor2 <- hou.tpcb.3$tPCB/hou.tpcb.3$predicted
factor2.tpcb.3 <- nrow(hou.tpcb.3[hou.tpcb.3$factor2 > 0.5 & hou.tpcb.3$factor2 < 2,
                                ])/length(hou.tpcb.3[,1])*100

