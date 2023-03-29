## Water PCB concentrations analysis.
## Data were obtained from EPA and contractors from PCB Superfund
## sites in USA. Data only AroclorCongener = Congener

# Install packages
install.packages("ggplot2")
install.packages('FactoMineR')
install.packages('factoextra')

# Load libraries
{
  library(ggplot2)
  library(FactoMineR) # Perform PCA
  library(factoextra) # Plot result from PCA
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

# PCB congener analysis ---------------------------------------------------
# Prepare data for PCA
# Remove congeners with < 50% detection frequency
prof.2 <- prof[, colMeans(!is.na(prof)) >= 0.5]
# Add SampleID names to row name
rownames(prof.2) <- cong$SampleID

# Perform PCA all samples
PCA <- PCA(prof.2, graph = FALSE)
fviz_eig(PCA, addlabels = TRUE, ylim = c(0, 100))
fviz_pca_var(PCA, col.var = "cos2",
             repel = TRUE) 
fviz_pca_ind(PCA, geom.ind = "point", pointshape = 21, 
             pointsize = 2, col.ind = "black", palette = "jco", 
             addEllipses = TRUE, label = "var",
             col.var = "black", repel = TRUE)

# Remove samples with less than x% of congeners
prof <- cbind(cong$SampleID, prof)

prop_non_na <- rowMeans(!is.na(prof))

# subset data frame to rows with >= 50% non-NA values
prof.2 <- prof[prop_non_na >= 0.8, ]


# Just looking at samples with Method 1668
{
  cong.1668 <- subset(wdc, EPAMethod == "M1668")
  # Remove samples with only 0s
  cong.1668 <- cong.1668[!(rowSums(cong.1668[, c(14:117)],
                                   na.rm = TRUE)==0), ]
  sampleID <- cong.1668$SampleID
  # Remove metadata
  cong.1668 <- subset(cong.1668,
                      select = -c(SampleID:AroclorCongener))
  # Remove Aroclor data
  cong.1668 <- subset(cong.1668, select = -c(A1016:A1260))
}
tmp <- rowSums(cong.1668, na.rm = TRUE)
prof.1668 <- sweep(cong.1668, 1, tmp, FUN = "/")
# Remove congeners with < 50% detection frequency
prof.1668.2 <- prof.1668[, colMeans(!is.na(prof.1668)) >= 0.5]
# Add SampleID names to row name
rownames(prof.1668.2) <- sampleID

PCA <- PCA(prof.1668.2, graph = FALSE)
fviz_eig(PCA, addlabels = TRUE, ylim = c(0, 30))
fviz_pca_ind(PCA, geom.ind = "point", pointshape = 21, 
             pointsize = 2, col.ind = "black", palette = "jco", 
             addEllipses = TRUE, label = "var",
             col.var = "black", repel = TRUE)




