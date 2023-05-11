## Water PCB concentrations data analysis per site
## Kalamazoo River
## Aroclors 1242, 1254 and 1260, no congener analysis

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
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotate("text", x = 2, y = 100, label = "Kalamazoo River",
           size = 3)

# (3) Seasonality
ggplot(kal.tpcb, aes(x = season, y = tPCB)) +
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
  annotate("text", x = 1.2, y = 10^6.5, label = "Kalamazoo River",
           size = 3)

# (4) Sites
ggplot(kal.tpcb, aes(x = factor(SiteID), y = tPCB)) + 
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
  annotate("text", x = 5, y = 10^6.5, label = "Kalamazoo River",
           size = 3)

# Remove site -------------------------------------------------------------
# Remove site PlainwellDam. Dredging = WCPCB-KAL023
kal.tpcb.1 <- subset(kal.tpcb, SiteID != c("WCPCB-KAL023"))

# Include USGS flow data --------------------------------------------------
# Include flow data from USGS station Kalamazoo River, no temperature available
{
  siteKalN1 <- "04108660" # KALAMAZOO RIVER AT NEW RICHMOND, MI
  siteKalN2 <- "04106000" # KALAMAZOO RIVER AT COMSTOCK, MI
  # Codes to retrieve data
  paramflow <- "00060" # discharge
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
}

# 2nd plots -----------------------------------------------------------
# (1.1) tPCB
hist(kal.tpcb.1$tPCB)
hist(log10(kal.tpcb.1$tPCB))
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
  ylab(expression(bold(atop("Water Concentration",
                            paste(Sigma*"PCB (pg/L)"))))) +
  theme(axis.text.y = element_text(face = "bold", size = 9),
        axis.title.y = element_text(face = "bold", size = 10)) +
  theme(axis.text.x = element_text(face = "bold", size = 9,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 9)) +
  annotate("text", x = 2, y = 100, label = "Kalamazoo River",
           size = 3)

ggplot(kal.tpcb.2, aes(y = tPCB,
                       x = format(date,'%Y'))) +
  geom_point(shape = 21, fill = "#66ccff") +
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
  annotate("text", x = 1, y = 100, label = "Kalamazoo River",
           size = 3)

# Regressions -------------------------------------------------------------
# Perform Linear Mixed-Effects Model (LMEM)
# Use kal.tpcb.2
tpcb <- kal.tpcb.2$tPCB
time <- kal.tpcb.2$time
site <- kal.tpcb.2$site.code
season <- kal.tpcb.2$season
flow <- kal.tpcb.2$flow.1

lme.kal.tpcb <- lmer(log10(tpcb) ~ 1 + time + season + flow + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"),
                  na.action = na.exclude)

# See results
summary(lme.kal.tpcb)
# Look at residuals
{
  res.kal.tpcb <- resid(lme.kal.tpcb) # get list of residuals
  # Create Q-Q plot for residuals
  # Create pdf file
  pdf("Output/Plots/Sites/Q-Q/KalamazooRiverQ-QtPCB.pdf")
  qqnorm(res.kal.tpcb,
         main = expression(paste("Normal Q-Q Plot (log"[10]* Sigma,
                                 "PCB)")))
  # Add a straight diagonal line to the plot
  qqline(res.kal.tpcb)
  dev.off()
}

# Create matrix to store results
{
  lme.tpcb <- matrix(nrow = 1, ncol = 24)
  lme.tpcb[1] <- fixef(lme.kal.tpcb)[1] # intercept
  lme.tpcb[2] <- summary(lme.kal.tpcb)$coef[1,"Std. Error"] # intercept error
  lme.tpcb[3] <- summary(lme.kal.tpcb)$coef[1,"Pr(>|t|)"] # intercept p-value
  lme.tpcb[4] <- fixef(lme.kal.tpcb)[2] # time
  lme.tpcb[5] <- summary(lme.kal.tpcb)$coef[2,"Std. Error"] # time error
  lme.tpcb[6] <- summary(lme.kal.tpcb)$coef[2,"Pr(>|t|)"] # time p-value
  lme.tpcb[7] <- fixef(lme.kal.tpcb)[3] # flow
  lme.tpcb[8] <- summary(lme.kal.tpcb)$coef[3,"Std. Error"] # flow error
  lme.tpcb[9] <- summary(lme.kal.tpcb)$coef[3,"Pr(>|t|)"] # flow p-value
  lme.tpcb[10] <- fixef(lme.kal.tpcb)[4] # temperature
  lme.tpcb[11] <- summary(lme.kal.tpcb)$coef[4,"Std. Error"] # temperature error
  lme.tpcb[12] <- summary(lme.kal.tpcb)$coef[4,"Pr(>|t|)"] # temperature p-value
  lme.tpcb[13] <- fixef(lme.kal.tpcb)[5] # season 2
  lme.tpcb[14] <- summary(lme.kal.tpcb)$coef[5,"Std. Error"] # season 2 error
  lme.tpcb[15] <- summary(lme.kal.tpcb)$coef[5,"Pr(>|t|)"] # season 2 p-value
  lme.tpcb[16] <- fixef(lme.kal.tpcb)[6] # season 3
  lme.tpcb[17] <- summary(lme.kal.tpcb)$coef[6,"Std. Error"] # season 3 error
  lme.tpcb[18] <- summary(lme.kal.tpcb)$coef[6,"Pr(>|t|)"] # season 3 p-value
  lme.tpcb[19] <- -log(2)/lme.tpcb[4]/365 # t0.5
  lme.tpcb[20] <- abs(-log(2)/lme.tpcb[4]/365)*lme.tpcb[5]/abs(lme.tpcb[4]) # t0.5 error
  lme.tpcb[21] <- as.data.frame(VarCorr(lme.kal.tpcb))[1,'sdcor']
  lme.tpcb[22] <- as.data.frame(r.squaredGLMM(lme.kal.tpcb))[1, 'R2m']
  lme.tpcb[23] <- as.data.frame(r.squaredGLMM(lme.kal.tpcb))[1, 'R2c']
  lme.tpcb[24] <- shapiro.test(resid(lme.kal.tpcb))$p.value
}

# Just 3 significant figures
lme.tpcb <- formatC(signif(lme.tpcb, digits = 3))
# Add column names
colnames(lme.tpcb) <- c("Intercept", "Intercept.error",
                        "Intercept.pv", "time", "time.error", "time.pv",
                        "flow", "flow.error", "flow.pv", "temperature",
                        "temperature.error", "temperature.pv", "season2",
                        "season2.error", "season2, pv", "season3",
                        "season3.error", "season3.pv", "t05", "t05.error",
                        "RandonEffectSiteStdDev", "R2nR", "R2R", "Normality")

# Export results
write.csv(lme.tpcb, file = "Output/Data/Sites/csv/KalamazooRiver/KalamazooLmetPCB.csv")

# Modeling plots
# (1) Get predicted values tpcb
fit.lme.values.kal.tpcb <- as.data.frame(fitted(lme.kal.tpcb))
# Add column name
colnames(fit.lme.values.kal.tpcb) <- c("predicted")
# Add predicted values to data.frame
kal.tpcb.2$predicted <- 10^(fit.lme.values.kal.tpcb$predicted)

# Plot prediction vs. observations, 1:1 line
p <- ggplot(kal.tpcb.2, aes(x = tPCB, y = predicted)) +
  geom_point(shape = 21, size = 3, fill = "#66ccff") +
  scale_y_log10(limits = c(10, 10^6), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(limits = c(10, 10^6), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab(expression(bold("Observed concentration " *Sigma*"PCB (pg/L)"))) +
  ylab(expression(bold("Predicted lme concentration " *Sigma*"PCB (pg/L)"))) +
  geom_abline(intercept = 0, slope = 1, col = "red", linewidth = 1.3) +
  geom_abline(intercept = 0.3, slope = 1, col = "blue", linewidth = 0.8) + # 1:2 line (factor of 2)
  geom_abline(intercept = -0.3, slope = 1, col = "blue", linewidth = 0.8) + # 2:1 line (factor of 2)
  theme_bw() +
  theme(aspect.ratio = 15/15) +
  annotation_logticks(sides = "bl") +
  annotate('text', x = 200, y = 10^5.8,
           label = expression(atop("Kalamazoo River (R"^2*"= 0.97)",
                                   paste("t"[1/2]*" = 3 ± 0.1 (yr)"))),
           size = 3, fontface = 2)
# See plot
print(p)
# Save plot
ggsave(filename = "Output/Plots/Sites/ObsPred/KalamazooRiver/KalamazooRiverObsPredtPCB.pdf",
       plot = p, device = "pdf")

# Plot residuals vs. predictions
{
  # Create pdf file
  pdf("Output/Plots/Sites/Residual/KalamazooRiverResidualtPCB.pdf")
  plot(log10(kal.tpcb.2$predicted), res.kal.tpcb,
       points(log10(kal.tpcb.2$predicted), res.kal.tpcb, pch = 16, 
              col = "#66ccff"),
       xlim = c(1, 6),
       ylim = c(-2, 2),
       xlab = expression(paste("Predicted lme concentration ",
                               Sigma, "PCB (pg/L)")),
       ylab = "Residual")
  abline(0, 0)
  abline(h = c(-1, 1), col = "grey")
  abline(v = seq(1, 6, 0.5), col = "grey")
  dev.off()
  }

# Estimate a factor of 2 between observations and predictions
kal.tpcb.2$factor2 <- kal.tpcb.2$tPCB/kal.tpcb.2$predicted
factor2.tpcb <- nrow(kal.tpcb.2[kal.tpcb.2$factor2 > 0.5 & kal.tpcb.2$factor2 < 2,
                                ])/length(kal.tpcb.2[,1])*100



