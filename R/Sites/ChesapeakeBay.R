## Water PCB concentrations data analysis per site
## Chesapeake Bay & Delaware Canal

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

# Select Chesapeake Bay & Delaware Canal data ---------------------------------------------------
che.0 <- wdc[str_detect(wdc$LocationName, 'Chesapeake Bay'),]

# Data preparation --------------------------------------------------------
{
  # Remove samples (rows) with total PCBs  = 0
  che.1 <- che.0[!(rowSums(che.0[, c(14:117)], na.rm = TRUE)==0),]
  # Calculate total PCB
  tpcb.che <- rowSums(che.1[, c(14:117)], na.rm = T)
  # Change date format
  che.1$SampleDate <- as.Date(che.1$SampleDate, format = "%m/%d/%y")
  # Calculate sampling time
  time.day <- data.frame(as.Date(che.1$SampleDate) - min(as.Date(che.1$SampleDate)))
  # Create individual code for each site sampled
  site.numb <- che.1$SiteID %>% as.factor() %>% as.numeric
  # Include season
  yq.s <- as.yearqtr(as.yearmon(che.1$SampleDate, "%m/%d/%Y") + 1/12)
  season.s <- factor(format(yq.s, "%q"), levels = 1:4,
                     labels = c("0", "S-1", "S-2", "S-3")) # winter, spring, summer, fall
  # Create data frame
  che.tpcb <- cbind(factor(che.1$SiteID), che.1$SampleDate,
                    che.1$Latitude, che.1$Longitude, as.matrix(tpcb.che),
                    data.frame(time.day), site.numb, season.s)
  # Add column names
  colnames(che.tpcb) <- c("SiteID", "date", "Latitude", "Longitude",
                          "tPCB", "time", "site.code", "season")
}

# Get coordinates per site to plot in Google Earth
che.location <- che.tpcb[c('SiteID', 'Latitude', 'Longitude', 'tPCB')]
# Average tPCB per site
che.location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                          data = che.location, mean)

# General plots -------------------------------------------------------------------
# (1) Histograms
hist(che.tpcb$tPCB)
hist(log10(che.tpcb$tPCB))

# (2) Time trend plots
ggplot(che.tpcb, aes(y = tPCB,
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
  annotate("text", x = 2, y = 20, label = "Chesapeake Bay",
           size = 3)

# (3) Seasonality
ggplot(che.tpcb, aes(x = season, y = tPCB)) +
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
  annotate("text", x = 1, y = 20, label = "Chesapeake Bay",
           size = 3)

# (4) Sites
ggplot(che.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 5, y = 20, label = "Chesapeake Bay",
           size = 3)

# tPCB Regressions --------------------------------------------------------
# Perform Linear Mixed-Effects Model (lme)
# Get variables
tpcb <- che.tpcb$tPCB
time <- che.tpcb$time
site <- che.tpcb$site.code
season <- che.tpcb$season
# tPCB vs. time + season + flow + temp + site
lme.che.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + (1|site),
                     REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))

# See results
summary(lme.che.tpcb)
# Look at residuals
{
  res.che.tpcb <- resid(lme.che.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.che.tpcb, main = "log10(C)")
  qqnorm(res.che.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.che.tpcb)
}
# Shapiro test
shapiro.test(res.che.tpcb)
# Random effect site Std Dev
RandonEffectSiteStdDev <- as.data.frame(VarCorr(lme.che.tpcb))[1,'sdcor']
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lme.che.tpcb))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lme.che.tpcb))[1, 'R2c']
# Extract coefficient values
time.coeff <- summary(lme.che.tpcb)$coef[2, "Estimate"]
time.coeff.ste <- summary(lme.che.tpcb)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -ln(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.che.tpcb <- as.data.frame(fitted(lme.che.tpcb))
# Add column name
colnames(fit.lme.values.che.tpcb) <- c("predicted")
# Add predicted values to data.frame
che.tpcb$predicted <- 10^(fit.lme.values.che.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(che.tpcb, aes(x = tPCB, y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "#66ccff") +
  scale_y_log10(limits = c(5, 10^5.5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(5, 10^5.5), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", linewidth = 1.3) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 50, y = 10^5,
           label = expression(atop(" Chesapeake Bay (R"^2*"= 0.45)",
                                   paste("t"[1/2]*" = 14 ± 4 (yr)"))),
           size = 3, fontface = 2)

# Plot residuals vs. predictions
{
  plot(log10(che.tpcb$predicted), res.che.tpcb,
       points(log10(che.tpcb$predicted), res.che.tpcb, pch = 16, 
              col = "#66ccff"),
       xlim = c(2, 5),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  abline(0, 0)
  abline(h = c(-1, 1), col = "grey")
  abline(v = seq(2, 5, 0.5), col = "grey")
  }

# Estimate a factor of 2 between observations and predictions
che.tpcb$factor2 <- che.tpcb$tPCB/che.tpcb$predicted
factor2.tpcb <- nrow(che.tpcb[che.tpcb$factor2 > 0.5 & che.tpcb$factor2 < 2,
                                ])/length(che.tpcb[,1])*100

# Plot time series with lme predictions
# Create a data frame to storage data
{
  time.serie.tpcb <- as.data.frame(matrix(nrow = length(che.tpcb[,1]),
                                          ncol = 3))
  # Add name to columns
  colnames(time.serie.tpcb) <- c('date', 'tPCB', 'lmetPCB')
  # Add data
  time.serie.tpcb$date <- che.tpcb$date
  time.serie.tpcb$tPCB <- che.tpcb$tPCB
  time.serie.tpcb$lmetPCB <- 10^(fit.lme.values.che.tpcb)
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
  scale_y_log10(limits = c(1, 10^6), breaks = trans_breaks("log10", function(x) 10^x),
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
  annotate("text", x = as.Date("2004-05-01", format = "%Y-%m-%d"),
           y = 10, label = "Chesapeake Bay", size = 3)

# Individual PCB Analysis -------------------------------------------------
# Prepare data.frame
{
  # Values coming from Data preparation section (che.tpcb)
  # Remove metadata
  che.pcb <- subset(che.1, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  che.pcb <- subset(che.pcb, select = -c(A1016:A1260))
  # Log10 individual PCBs 
  che.pcb <- log10(che.pcb)
  # Replace -inf to NA
  che.pcb <- do.call(data.frame,
                     lapply(che.pcb,
                            function(x) replace(x, is.infinite(x), NA)))
  # Remove individual PCB that have 30% or less NA values
  che.pcb.1 <- che.pcb[,
                       -which(colSums(is.na(che.pcb))/nrow(che.pcb) > 0.7)]
  # Add site ID
  che.pcb.1$SiteID <- che.1$SiteID
  # Add SampleDate
  che.pcb.1$SampleDate <- che.1$SampleDate
  # Add sampling time
  che.pcb.1$time <- che.tpcb$time
  # Add sampling site code
  che.pcb.1$site.numb <- che.tpcb$site.code
  # Add season
  che.pcb.1$season <- che.tpcb$season
  # Remove metadata for analysis
  che.pcb.2 <- subset(che.pcb.1, select = -c(SiteID:season))
}

# LME for individual PCBs -------------------------------------------------
# Get covariates
time <- che.pcb.1$time
season <- che.pcb.1$season
site <- che.pcb.1$site.numb

# Create matrix to store results
lme.pcb <- matrix(nrow = length(che.pcb.2[1,]), ncol = 21)

# Perform LME
for (i in 1:length(che.pcb.2[1,])) {
  fit <- lmer(che.pcb.2[,i] ~ 1 + time + season + (1|site),
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
  lme.pcb[i,7] <- fixef(fit)[3] # # season 1
  lme.pcb[i,8] <- summary(fit)$coef[3,"Std. Error"] # season 1 error
  lme.pcb[i,9] <- summary(fit)$coef[3,"Pr(>|t|)"] # # season 1 p-value
  lme.pcb[i,10] <- fixef(fit)[4] # season 2
  lme.pcb[i,11] <- summary(fit)$coef[4,"Std. Error"] # season 2 error
  lme.pcb[i,12] <- summary(fit)$coef[4,"Pr(>|t|)"] # season 2 p-value
  lme.pcb[i,13] <- fixef(fit)[5] # season 3
  lme.pcb[i,14] <- summary(fit)$coef[5,"Std. Error"] # season 3 error
  lme.pcb[i,15] <- summary(fit)$coef[5,"Pr(>|t|)"] # season 3 p-value
  lme.pcb[i,16] <- -log(2)/lme.pcb[i,4]/365 # t0.5
  lme.pcb[i,17] <- abs(-log(2)/lme.pcb[i,4]/365)*lme.pcb[i,5]/abs(lme.pcb[i,4]) # t0.5 error
  lme.pcb[i,18] <- as.data.frame(VarCorr(fit))[1,'sdcor']
  lme.pcb[i,19] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2m']
  lme.pcb[i,20] <- as.data.frame(r.squaredGLMM(fit))[1, 'R2c']
  lme.pcb[i,21] <- shapiro.test(resid(fit))$p.value
}

# Just 3 significant figures
lme.pcb <- formatC(signif(lme.pcb, digits = 3))
# Add congener names
congeners <- colnames(che.pcb.2)
lme.pcb <- as.data.frame(cbind(congeners, lme.pcb))
# Add column names
colnames(lme.pcb) <- c("Congeners", "Intercept", "Intercept.error",
                       "Intercept.pv", "time", "time.error", "time.pv",
                       "season1", "season1.error", "season1.pv", "season2",
                       "season2.error", "season2, pv", "season3",
                       "season3.error", "season3.pv", "t05", "t05.error",
                       "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality")
# Remove congeners with no normal distribution
# Shapiro test p-value < 0.05
lme.pcb <- lme.pcb[lme.pcb$Normality > 0.05, ]

# Export results
write.csv(lme.pcb, file = "Output/Data/Sites/csv/ChesapeakeLmePCB.csv")

# Generate predictions
# Remove congeners with no Normality
che.pcb.3 <- select(che.pcb.2, -PCB20.21.28.31.33.50.53, -PCB40.41.64.71.72,
                    -PCB61.66.70.74.76.93.95.98.100.102, -PCB180.193)
# Create matrix to store results
lme.fit.pcb <- matrix(nrow = length(che.pcb.3[,1]),
                      ncol = length(che.pcb.3[1,]))

for (i in 1:length(che.pcb.3[1,])) {
  fit <- lmer(che.pcb.3[,i] ~ 1 + time + season + (1|site),
              REML = FALSE,
              control = lmerControl(check.nobs.vs.nlev = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE="ignore"),
              na.action = na.exclude)
  lme.fit.pcb[,i] <- fitted(fit)
}

# Estimate a factor of 2 between observations and predictions
factor2 <- 10^(che.pcb.3)/10^(lme.fit.pcb)
factor2.pcb <- sum(factor2 > 0.5 & factor2 < 2,
                   na.rm = TRUE)/(sum(!is.na(factor2)))*100

# Plot 1:1 for all congeners
# Add SiteID to both data.frames
che.pcb.4 <- cbind(che.pcb.1$SiteID, che.pcb.3)
lme.fit.pcb.2 <- as.data.frame(cbind(che.pcb.1$SiteID, lme.fit.pcb))
# Add/fix column name to SiteID
colnames(che.pcb.4)[1] <- 'SiteID'
colnames(lme.fit.pcb.2)[1] <- 'SiteID'
# Merge both new data.frames
pcb.plot <- merge(che.pcb.4, lme.fit.pcb.2, by = "SiteID")



# Selected individual PCB regression --------------------------------------
# lme
lme.che.pcbi <- lmer(che.pcb.2$PCB56.60 ~ 1 + time + season +
                       (1|site), REML = FALSE,
                     control = lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"),
                     na.action = na.exclude)

# See results
summary(lme.che.pcbi)
# Look at residuals
{
  res.lme.che.pcbi <- resid(lme.che.pcbi) # get list of residuals
  # Create Q-Q plot for residuals
  qqnorm(res.lme.che.pcbi, main = expression(paste("Normal Q-Q Plot PCBs 56+60")))
  # Add a straight diagonal line to the plot
  qqline(res.lme.che.pcbi)
}
# Shapiro test
shapiro.test(res.lme.che.pcbi)
# Extract R2 no random effect
R2.nre <- as.data.frame(r.squaredGLMM(lme.che.pcbi))[1, 'R2m']
# Extract R2 with random effect
R2.re <- as.data.frame(r.squaredGLMM(lme.che.pcbi))[1, 'R2c']
# Extract coefficient values
time.coeff <- summary(lme.che.pcbi)$coef[2, "Estimate"]
time.coeff.ste <- summary(lme.che.pcbi)$coef[2, "Std. Error"]
# Calculate half-life tPCB in yr (-log(2)/slope/365)
t0.5 <- -log(2)/time.coeff/365 # half-life tPCB in yr = -ln(2)/slope/365
# Calculate error
t0.5.error <- abs(t0.5)*time.coeff.ste/abs(time.coeff)

# (1) Get predicted values pcbi
date.pcbi <- format(che.pcb.1$SampleDate, "%Y-%m-%d")
obs <- che.pcb.2$PCB56.60
che.pcbi <- cbind(date.pcbi, obs)
fit.lme.values.pcbi <- as.data.frame(fitted(lme.che.pcbi))
che.pcbi <- cbind(che.pcbi, fit.lme.values.pcbi)
colnames(che.pcbi) <- c("date", "obs", 'lme')
che.pcbi$date <- as.Date(che.pcbi$date)
che.pcbi$obs <- as.numeric(che.pcbi$obs)

# Plot residuals vs. predictions
{
  plot(che.pcbi$lme, res.lme.che.pcbi,
       points(che.pcbi$lme, res.lme.che.pcbi, pch = 16, 
              col = "#66ccff"),
       xlim = c(0, 4),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted concentration PCBs 56+60 (pg/L)")),
       ylab = "Residual (lme)")
  abline(0, 0)
  abline(h = seq(-2, 2, 1), col = "grey")
  abline(v = seq(0, 4, 0.5), col = "grey")
}

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.che.pcbi <- as.data.frame(fitted(lme.che.pcbi))
# Add column name
colnames(fit.lme.values.che.pcbi) <- c("predicted")
# Add predicted values to data.frame
che.pcb.2$predicted <- 10^(fit.lme.values.che.pcbi$predicted)

# Plot prediction vs. observations, 1:1 line
ggplot(che.pcb.2, aes(x = 10^(PCB56.60), y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "#66ccff") +
  scale_y_log10(limits = c(0.1, 100000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(0.1, 100000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration PCBs 56+60 (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration PCBs 56+60 (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", linewidth = 1.3) +
  geom_abline(intercept = log10(2), slope = 1, col = "blue", linewidth = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = log10(0.5), slope = 1, col = "blue", linewidth = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 5, y = 10^4.5,
           label = expression(atop("Chesapeake Bay (R"^2*" = 0.73)",
                                   paste("t"[1/2]*" = 0 ± 5.7x10"^-0.5*" (yr)"))),
           size = 3, fontface = 2)

# Modeling plots
# Change data.frame format to be plotted
che.pcbi.2 <- melt(che.pcbi, id.vars = c("date"))
# Plot
ggplot(che.pcbi.2, aes(x = date, y = 10^(value), group = variable)) +
  geom_point(aes(shape = variable, color = variable, size = variable,
                 fill = variable)) +
  scale_shape_manual(values = c(21, 3)) +
  scale_color_manual(values = c('black','#8856a7')) +
  scale_size_manual(values = c(2, 1)) +
  scale_fill_manual(values = c("#1b98e0", '#8856a7')) +
  scale_x_date(labels = date_format("%Y-%m")) +
  scale_y_log10(limits = c(0.1, 100000)) +
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
  annotate("text", x = as.Date("2004-06-01", format = "%Y-%m-%d"),
           y = 1, label = "Chesapeake Bay", size = 3)

