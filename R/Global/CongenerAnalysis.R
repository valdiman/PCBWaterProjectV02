## Water PCB concentrations analysis.
## Data were obtained from EPA and contractors from PCB Superfund
## sites in USA. Data only AroclorCongener = Congener

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

# Data preparation --------------------------------------------------------
# # Only consider congener data 
{
  cong <- subset(wdc, AroclorCongener == "Congener")
  # Remove samples with only 0s
  cong <- cong[!(rowSums(cong[, c(14:117)], na.rm = TRUE)==0), ]
  # Remove metadata
  cong.1 <- subset(cong, select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  cong.1 <- subset(cong.1, select = -c(A1016:A1260))
}

# Create PCB profile distribution
{
  # Generate PCB profile for individual samples
  tmp <- rowSums(cong.1, na.rm = TRUE)
  prof <- sweep(cong.1, 1, tmp, FUN = "/")
  # Generate average PCB profile
  prof.ave <- data.frame(colMeans(prof, na.rm = TRUE))
  colnames(prof.ave) <- c("mean")
  prof.sd <- data.frame(apply(prof, 2, sd, na.rm = TRUE))
  colnames(prof.sd) <- c("sd")
  congener <- row.names(prof.ave)
  prof.ave <- cbind(congener, prof.ave$mean, prof.sd$sd)
  colnames(prof.ave) <- c("congener", "mean", "sd")
  prof.ave <- data.frame(prof.ave)
  prof.ave$mean <- as.numeric(as.character(prof.ave$mean))
  prof.ave$sd <- as.numeric(as.character(prof.ave$sd))
  prof.ave$congener <- as.character(prof.ave$congener)
  # Then turn it back into a factor with the levels in the correct order
  prof.ave$congener <- factor(prof.ave$congener,
                              levels = unique(prof.ave$congener))
}

# Plot average PCB profile
ggplot(prof.ave, aes(x = congener, y = mean)) +
  geom_bar(position = position_dodge(), stat = "identity",
           fill = "black") +
  geom_errorbar(aes(ymin = mean, ymax = (mean+sd)), width = 0.9,
                position = position_dodge(0.9)) +
  xlab("") +
  ylim(0, 0.5) +
  theme_bw() +
  theme(aspect.ratio = 4/12) +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 8),
        axis.title.y = element_text(face = "bold", size = 9)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  annotate("text", x = 5.4, y = 0.4, label = "PCBs 4+10", size = 2,
           fontface = 1, angle = 90) +
  annotate("text", x = 8, y = 0.36, label = "PCB 11", size = 2,
           fontface = 1, angle = 90) +
  annotate("text", x = 14, y = 0.32, label = "PCBs 18+30", size = 2,
           fontface = 1, angle = 90) +
  annotate("text", x = 18, y = 0.33, label = "PCBs 20+21+28+\n31+33+50+53",
           size = 2, fontface = 1, angle = 90) +
  annotate("text", x = 30.5, y = 0.25, label = "PCBs 43+49+\n52+69+73",
           size = 2, fontface = 1, angle = 90) +
  annotate("text", x = 39.3, y = 0.30, label = "PCBs 61+66+70+74+\n76+93+95+98+100+102",
           size = 2, fontface = 1, angle = 90) +
  annotate("text", x = 70.5, y = 0.33, label = "PCBs 132+153+\n161+168",
           size = 2, fontface = 1, angle = 90)

# Prepare data for PCA
t.prof <- data.frame(t(prof))

# Add column names with samples name
colnames(t.prof) <- cong$SiteID
t.prof <- t(t.prof)

is.na(t.prof)
is.infinite(t.prof)


# Perform PCA all samples
PCA <- prcomp(t.prof, na.action = na.pass)
summary(PCA)

# Add Aroclor data

# PCA analysis

# Cosine analysis


