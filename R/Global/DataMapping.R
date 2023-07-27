## Water PCB concentrations mapping.
# Data were obtained from EPA and contractors from PCB Superfund
# sites in USA

# Install packages
install.packages("ggplot2")
install.packages("devtools")
install.packages("dplyr")
install.packages("stringr")
install.packages("maps")
install.packages("mapdata")
install.packages("ggmap")
install.packages("usethis")
install.packages("GISTools")
install.packages("rgeos")
install.packages("ggsn")
install.packages("sf")
install.packages("ggrepel")
install.packages("ggpp")
install.packages("raster")
install.packages("grid")

# Load libraries
{
  library(dplyr)
  library(usethis)
  library(devtools)
  library(ggplot2)
  library(ggmap) # function map_data
  library(maps)
  library(leaflet)
  library(rgeos)
  library(ggsn)
  library(ggrepel)
  library(reshape2)
  library(ggpmisc)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc.0 <- read.csv("Data/WaterDataCongenerAroclor08052022.csv")
# Extract sample site locations -------------------------------------------
# Calculate total PCB
# Data preparation
{
  # Remove samples (rows) with total PCBs  = 0
  wdc.2 <- wdc.0[!(rowSums(wdc.0[, c(14:117)], na.rm = TRUE)==0),]
  # Calculate total PCB
  tpcb <- rowSums(wdc.2[, c(14:117)], na.rm = T)
  # Select and combine sample sites anf tPCB
  wdc.location <- cbind.data.frame(wdc.2$SiteID, wdc.2$Latitude, wdc.2$Longitude,
                                   tpcb)
  # Name the columns
  colnames(wdc.location) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Average tPCB per sample site
  wdc.location <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                            data = wdc.location, mean)
}

# Global maps -------------------------------------------------------------
us <- map_data("usa")
states <- map_data("state")

# Find number of samples per state to be included as table in maps
wdc.1 <- wdc.0 %>%
  group_by(StateSampled) %>%
  summarise(n = n())
#colnames(wdc.1) <- c("State", "# samples")
wdc.2 <- data.frame(t(wdc.1))
name <- c('State', '# samples')
wdc.2 <- data.frame(col1 = name, wdc.2)
names(wdc.2) <- NULL

# (1) Map of US with locations
ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = "lightblue") +
  coord_fixed(1.3) +
  theme_nothing() +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_path(data = states, aes(x = long, y = lat, group = group),
             colour = "white") +
  geom_polygon(color = "black", fill = NA) +
  geom_point(data = wdc.0, aes(x = Longitude, y = Latitude),
             color = "black",
             size = 1.2, shape = 20) +
  annotate(geom = 'table', x = -65, y = 53,
           label = list(wdc.2), size = 2.9) # add table with info

# (2) Map + tPCB
# Cannot include legend
ggplot() +
  geom_polygon(data = us, aes(x = long, y = lat, group = group),
               color = "black", fill = "lightblue") +
  coord_fixed(1.3) +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_path(data = states, aes(x = long, y = lat, group = group),
            colour = "white") +
  geom_polygon(color = "black", fill = NA) +
  geom_point(data = wdc.location, aes(x = Longitude, y = Latitude,
                                       size = tPCB),
             color = "red") +
  theme(legend.position = "right") +
  scale_size_area(breaks = c(1000, 50*1000, 500*1000, 1000*1000, 1500*1000,
                             2000*1000), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (SiteID mean)",
                                              paste("1990-2020 (pg/L)")))),
                  max_size = 5)
  
# Specific locations ------------------------------------------------------
# Portland Harbor ---------------------------------------------------------
# Select only from Portland Harbor
wdc.PO <- subset(wdc.0, LocationName == "Portland Harbor")

# Create general map
PO.box <- make_bbox(lon = wdc.PO$Longitude, lat = wdc.PO$Latitude, f = 0.8)
PO.map <- get_stamenmap(bbox = PO.box, zoom = 10)

# Plot map with sites
# Prepare data
# Get tPCB and coordinates
tPCB.PO <- data.frame(cbind(wdc.PO$LocationName, wdc.PO$Latitude,
                            wdc.PO$Longitude, rowSums(wdc.PO[, c(14:117)],
                                                      na.rm = TRUE)))
# Name the columns
colnames(tPCB.PO) <- c("LocationID", "Latitude", "Longitude", "tPCB")
# Change no numeric to numeric
tPCB.PO$Latitude <- as.numeric(tPCB.PO$Latitude)
tPCB.PO$Longitude <- as.numeric(tPCB.PO$Longitude)
tPCB.PO$tPCB <- as.numeric(tPCB.PO$tPCB)
# Average tPCB per site
tPCB.PO.mean <- aggregate(tPCB ~ LocationID + Latitude + Longitude,
                       data = tPCB.PO, FUN = mean)

# (1) Plot map + locations
ggmap(PO.map) +
  geom_point(data = tPCB.PO.mean, aes(x = Longitude, y = Latitude),
             shape = 21, color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = LocationID),
                   data = tPCB.PO.mean, family = 'Times New Roman', size = 3, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(PO.map) +
  geom_point(data = tPCB.PO.mean, aes(x = Longitude, y = Latitude,
                                      size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_size_area(breaks = c(100, 250, 500, 1000, 1500), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean)",
                                              paste("2007-2019 (pg/L)")))),
                  max_size = 8)
  

# Fox River ---------------------------------------------------------------
# Select only from Fox River
wdc.Fox <- subset(wdc.0, SiteName == "Fox River")

# Create general map
Fox.box <- make_bbox(lon = wdc.Fox$Longitude, lat = wdc.Fox$Latitude,
                    f = 0.4)
Fox.map <- get_stamenmap(bbox = Fox.box, zoom = 10)

# Plot map with sites
# Prepare data
# Remove samples (rows) with total PCBs  = 0
wdc.Fox.1 <- wdc.Fox[!(rowSums(wdc.Fox[, c(14:117)],
                             na.rm = TRUE)==0),] # sum of PCB1 to PCB209
# Get tPCB and coordinates
tPCB.Fox <- data.frame(cbind(wdc.Fox.1$LocationID, wdc.Fox.1$Latitude,
                             wdc.Fox.1$Longitude, rowSums(wdc.Fox.1[, c(14:117)],
                                                          na.rm = TRUE)))
# Name the columns
colnames(tPCB.Fox) <- c("LocationID", "Latitude", "Longitude", "tPCB")
# Change no numeric to numeric
tPCB.Fox$Latitude <- as.numeric(tPCB.Fox$Latitude)
tPCB.Fox$Longitude <- as.numeric(tPCB.Fox$Longitude)
tPCB.Fox$tPCB <- as.numeric(tPCB.Fox$tPCB)
# Average tPCB per site
tPCB.Fox.mean <- aggregate(tPCB ~ LocationID + Latitude + Longitude,
                       data = tPCB.Fox, FUN = mean)

# (1) Plot map + locations
ggmap(Fox.map) +
  geom_point(data = tPCB.Fox.mean, aes(x = Longitude, y = Latitude),
             shape = 21, color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotate('text', x = -87.5, y = 45,
           label = 'Green Bay', colour = 'black', size = 3,
           fontface = 2) +
  annotate('text', x = -88.2, y = 44.435,
           label = 'Fox River', colour = 'black', size = 3,
           fontface = 2) +
  annotate('text', x = -87.85, y = 44.7,
           label = 'Green Bay', colour = 'black', size = 3,
           fontface = 2) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = LocationID),
                   data = tPCB.mean, family = 'Times New Roman', size = 2.8, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(Fox.map) +
  geom_point(data = tPCB.Fox.mean, aes(x = Longitude, y = Latitude,
                                      size = tPCB), alpha = 1,
             color  = "red") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_size_area(breaks = c(100, 1500, 30000, 65000),
                  name = expression(bold(atop(Sigma*"PCBs (mean)",
                                              paste("2005-2018 (ng/L)")))),
                  max_size = 8)

# Hudson River ------------------------------------------------------------
# Select only from Hudson River
wdc.Hud <- subset(wdc.0, SiteName == "Hudson River")

# Create general map
Hud.box <- make_bbox(lon = wdc.Hud$Longitude, lat = wdc.Hud$Latitude, f = 0.3)
Hud.map <- get_stamenmap(bbox = Hud.box, zoom = 8)

# Plot map with sites
# Prepare data
# Get tPCB and coordinates
tPCB.Hud <- data.frame(cbind(wdc.Hud$LocationID, wdc.Hud$Latitude,
                             wdc.Hud$Longitude, rowSums(wdc.Hud[, c(14:117)],
                                                        na.rm = TRUE)))
# Name the columns
colnames(tPCB.Hud) <- c("LocationID", "Latitude", "Longitude", "tPCB")
# Change no numeric to numeric
tPCB.Hud$Latitude <- as.numeric(tPCB.Hud$Latitude)
tPCB.Hud$Longitude <- as.numeric(tPCB.Hud$Longitude)
tPCB.Hud$tPCB <- as.numeric(tPCB.Hud$tPCB)
# Average tPCB per site
tPCB.Hud.mean <- aggregate(tPCB ~ LocationID + Latitude + Longitude,
                       data = tPCB.Hud, FUN = mean)

# (1) Plot map + locations
ggmap(Hud.map) +
  geom_point(data = tPCB.Hud.mean, aes(x = Longitude, y = Latitude),
             shape = 21, color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = LocationID),
                   data = tPCB.mean, family = 'Times', size = 1.8, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(Hud.map) +
  geom_point(data = tPCB.Hud.mean, aes(x = Longitude, y = Latitude,
                                       size = tPCB), alpha = 1,
             color  = "red") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_size_area(breaks = c(1000, 10000, 30000, 40000),
                  name = expression(bold(atop(Sigma*"PCBs (mean)",
                                              paste("2005-2017 (ng/L)")))),
                  max_size = 8)


# Housatonic River --------------------------------------------------------
# Select only from Housotonic River
wdc.Hou <- subset(wdc.0, SiteName == "Housatonic River")

# Create general map
Hou.box <- make_bbox(lon = wdc.Hou$Longitude, lat = wdc.Hou$Latitude, f = 0.2)
Hou.map <- get_stamenmap(bbox = Hou.box, zoom = 10)

# Plot map with sites
# Prepare data
# Get tPCB and coordinates
tPCB.Hou <- data.frame(cbind(wdc.Hou$LocationID, wdc.Hou$Latitude,
                             wdc.Hou$Longitude, rowSums(wdc.Hou[, c(14:117)],
                                    na.rm = TRUE)))
# Name the columns
colnames(tPCB.Hou) <- c("LocationID", "Latitude", "Longitude", "tPCB")
# Change no numeric to numeric
tPCB.Hou$Latitude <- as.numeric(tPCB.Hou$Latitude)
tPCB.Hou$Longitude <- as.numeric(tPCB.Hou$Longitude)
tPCB.Hou$tPCB <- as.numeric(tPCB.Hou$tPCB)
# Average tPCB per site
tPCB.mean <- aggregate(tPCB ~ LocationID + Latitude + Longitude,
                       data = tPCB.Hou, FUN = mean)

# (1) Plot map + locations
ggmap(Hou.map) +
  geom_point(data = tPCB.mean, aes(x = Longitude, y = Latitude),
             shape = 21, color = "red",
             fill = "white", size = 1.75, stroke = 0.75)
  #geom_label_repel(aes(x = Longitude, y = Latitude, label = Site),
  #                 data = tPCB.mean, family = 'Times', size = 1.8, 
  #                 box.padding = 0.2, point.padding = 0.3,
  #                 segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(HR.map) +
  geom_point(data = tPCB.mean, aes(x = Longitude, y = Latitude,
                                   size = tPCB), alpha = 1, color  = "red") +
  scale_size_area(breaks = c(100, 125, 150, 175, 200),
                  name = "Ave. PCBs \n2018-2019 (pg/L)") +
  xlab("Longitude") +
  ylab("Latitude")

# Kalamazoo River ---------------------------------------------------------
# Select only from Kalamazoo  River
wdc.Kal <- subset(wdc.0, SiteName == "KalamazooRiver")

# Create general map
Kal.box <- make_bbox(lon = wdc.Kal$Longitude, lat = wdc.Kal$Latitude, f = 0.2)
Kal.map <- get_stamenmap(bbox = Kal.box, zoom = 10)

# Plot map with sites
# Prepare data
# Remove samples (rows) with total PCBs  = 0
wdc.Kal.1 <- wdc.Kal[!(rowSums(wdc.Kal[, c(12:115)],
                             na.rm = TRUE)==0),] # sum of PCB1 to PCB209
# Get tPCB and coordinates
tPCB.Kal <- data.frame(cbind(wdc.Kal.1$Latitude, wdc.Kal.1$Longitude,
                            wdc.Kal.1$SiteSampled,
                            rowSums(wdc.Kal.1[, c(12:115)],
                                    na.rm = TRUE)))
# Name the columns
colnames(tPCB.Kal) <- c("Latitude", "Longitude", "Site", "tPCB")
# Change no numeric to numeric
tPCB.Kal$Latitude <- as.numeric(tPCB.Kal$Latitude)
tPCB.Kal$Longitude <- as.numeric(tPCB.Kal$Longitude)
tPCB.Kal$tPCB <- as.numeric(tPCB.Kal$tPCB)
# Average tPCB per site
tPCB.mean <- aggregate(tPCB ~ Site + Latitude + Longitude,
                       data = tPCB.Kal, FUN = mean)

# (1) Plot map + locations
ggmap(Kal.map) +
  geom_point(data = tPCB.mean, aes(x = Longitude, y = Latitude), shape = 21,
             color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
geom_label_repel(aes(x = Longitude, y = Latitude, label = Site),
                 data = tPCB.mean, family = 'Times', size = 1.8, 
                 box.padding = 0.2, point.padding = 0.3,
                 segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(Kal.map) +
  geom_point(data = tPCB.mean, aes(x = Longitude, y = Latitude,
                                   size = tPCB), alpha = 1, color  = "red") +
  scale_size_area(breaks = c(100, 125, 150, 175, 200),
                  name = "Ave. PCBs \n2018-2019 (pg/L)") +
  xlab("Longitude") +
  ylab("Latitude")

# New Bedford -------------------------------------------------------------
# Select only MA
wdc.MA <- subset(wdc.0, StateSampled == "MA")
# Select only from Hudson River
wdc.NB <- subset(wdc.MA, SiteName == "NewBedford")

# Create general map
NB.box <- make_bbox(lon = wdc.NB$Longitude, lat = wdc.NB$Latitude, f = 0.3)
NB.map <- get_stamenmap(bbox = NB.box, zoom = 13)

# Plot map with sites
# Prepare data
# Remove samples (rows) with total PCBs  = 0
wdc.NB.1 <- wdc.NB[!(rowSums(wdc.NB[, c(12:115)],
                             na.rm = TRUE)==0),] # sum of PCB1 to PCB209
# Get tPCB and coordinates
tPCB.NB <- data.frame(cbind(wdc.NB.1$Latitude, wdc.NB.1$Longitude,
                            wdc.NB.1$SiteSampled,
                            rowSums(wdc.NB.1[, c(12:115)],
                                    na.rm = TRUE)))
# Name the columns
colnames(tPCB.NB) <- c("Latitude", "Longitude", "Site", "tPCB")
# Change no numeric to numeric
tPCB.NB$Latitude <- as.numeric(tPCB.NB$Latitude)
tPCB.NB$Longitude <- as.numeric(tPCB.NB$Longitude)
tPCB.NB$tPCB <- as.numeric(tPCB.NB$tPCB)
# Average tPCB per site
tPCB.mean <- aggregate(tPCB ~ Latitude + Longitude + Site,
                       data = tPCB.NB, FUN = mean)

# (1) Plot map + locations
ggmap(NB.map) +
  geom_point(data = tPCB.mean, aes(x = Longitude, y = Latitude), shape = 21,
             color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
geom_label_repel(aes(x = Longitude, y = Latitude, label = Site),
                 data = tPCB.mean, family = 'Times', size = 1.8, 
                 box.padding = 0.2, point.padding = 0.3,
                 segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(NB.map) +
  geom_point(data = tPCB.mean, aes(x = Longitude, y = Latitude,
                                   size = tPCB), alpha = 1, color  = "red") +
  #scale_size_area(breaks = c(100, 125, 150, 175, 200),
  #                name = "Ave. PCBs \n2018-2019 (pg/L)") +
  xlab("Longitude") +
  ylab("Latitude")

# Spokane River -----------------------------------------------------------
# Select only ID
wdc.ID <- subset(wdc.0, StateSampled == "ID")
# Select only from Hudson River
wdc.Spo <- subset(wdc.0, SiteName == "SpokaneRiver")

# Create general map
Spo.box <- make_bbox(lon = wdc.Spo$Longitude, lat = wdc.Spo$Latitude, f = 0.6)
Spo.map <- get_stamenmap(bbox = Spo.box, zoom = 10)

# Plot map with sites
# Prepare data
# Remove samples (rows) with total PCBs  = 0
wdc.Spo.1 <- wdc.Spo[!(rowSums(wdc.Spo[, c(12:115)],
                             na.rm = TRUE)==0),] # sum of PCB1 to PCB209
# Get tPCB and coordinates
tPCB.Spo <- data.frame(cbind(wdc.Spo.1$Latitude, wdc.Spo.1$Longitude,
                            wdc.Spo.1$SiteSampled,
                            rowSums(wdc.Spo.1[, c(12:115)],
                                    na.rm = TRUE)))
# Name the columns
colnames(tPCB.Spo) <- c("Latitude", "Longitude", "Site", "tPCB")
# Change no numeric to numeric
tPCB.Spo$Latitude <- as.numeric(tPCB.Spo$Latitude)
tPCB.Spo$Longitude <- as.numeric(tPCB.Spo$Longitude)
tPCB.Spo$tPCB <- as.numeric(tPCB.Spo$tPCB)
# Average tPCB per site
tPCB.mean <- aggregate(tPCB ~ Latitude + Longitude + Site,
                       data = tPCB.Spo, FUN = mean)

# (1) Plot map + locations
ggmap(Spo.map) +
  geom_point(data = tPCB.mean, aes(x = Longitude, y = Latitude), shape = 21,
             color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = Site),
                   data = tPCB.mean, family = 'Times', size = 1.8, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(Spo.map) +
  geom_point(data = tPCB.mean, aes(x = Longitude, y = Latitude,
                                   size = tPCB), alpha = 1, color  = "red") +
  #scale_size_area(breaks = c(100, 125, 150, 175, 200),
  #                name = "Ave. PCBs \n2018-2019 (pg/L)") +
  xlab("Longitude") +
  ylab("Latitude")

# Blue River --------------------------------------------------------------
# Select only from Blue River
wdc.Blu <- subset(wdc.0, SiteName == "BlueRiver")

# Create general map
Blu.box <- make_bbox(lon = wdc.Blu$Longitude, lat = wdc.Blu$Latitude, f = 1.0)
Blu.map <- get_stamenmap(bbox = Blu.box, zoom = 14)

# Plot map with sites
# Prepare data
# Remove samples (rows) with total PCBs  = 0
wdc.Blu.1 <- wdc.Blu[!(rowSums(wdc.Blu[, c(12:115)],
                               na.rm = TRUE)==0),] # sum of PCB1 to PCB209
# Get tPCB and coordinates
tPCB.Blu <- data.frame(cbind(wdc.Blu.1$Latitude, wdc.Blu.1$Longitude,
                             wdc.Blu.1$SiteSampled,
                             rowSums(wdc.Blu.1[, c(12:115)],
                                     na.rm = TRUE)))
# Name the columns
colnames(tPCB.Blu) <- c("Latitude", "Longitude", "Site", "tPCB")
# Change no numeric to numeric
tPCB.Blu$Latitude <- as.numeric(tPCB.Blu$Latitude)
tPCB.Blu$Longitude <- as.numeric(tPCB.Blu$Longitude)
tPCB.Blu$tPCB <- as.numeric(tPCB.Blu$tPCB)
# Average tPCB per site
tPCB.mean <- aggregate(tPCB ~ Latitude + Longitude + Site,
                       data = tPCB.Blu, FUN = mean)

# (1) Plot map + locations
ggmap(Blu.map) +
  geom_point(data = tPCB.mean, aes(x = Longitude, y = Latitude), shape = 21,
             color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = Site),
                   data = tPCB.mean, family = 'Times', size = 1.8, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(Blu.map) +
  geom_point(data = tPCB.mean, aes(x = Longitude, y = Latitude,
                                   size = tPCB), alpha = 1, color  = "red") +
  #scale_size_area(breaks = c(100, 125, 150, 175, 200),
  #                name = "Ave. PCBs \n2018-2019 (pg/L)") +
  xlab("Longitude") +
  ylab("Latitude")


# Extra -------------------------------------------------------------------
# Prepare congener data for plotting
# Get coordinates per site
LW <- subset(w.WI, w.WI$SiteSampled == 'LakeWinnebago')
LW <- data.frame(c(LW[1,6], LW[1,7]))
OU1 <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit1')
OU1 <- data.frame(c(OU1[1,6], OU1[1,7]))
OU2A <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit2A')
OU2A <- data.frame(c(OU2A[1,6], OU2A[1,7]))
OU2B <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit2B')
OU2B <- data.frame(c(OU2B[1,6], OU2B[1,7]))
OU2C <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit2C')
OU2C <- data.frame(c(OU2C[1,6], OU2C[1,7]))
OU3 <- subset(w.WI, w.WI$SiteSampled == 'OperableUnit3')
OU3 <- data.frame(c(OU3[1,6], OU3[1,7]))
wi.coord <- rbind(LW, OU1, OU2A, OU2B, OU2C, OU3)

# Total PCBs
# # remove samples (rows) with total PCBs  = 0
w.WI.t <- w.WI[!(rowSums(w.WI[,
                           c(12:115)],
                         na.rm = TRUE)==0),] # sum of PCB1 to PCB209
site.sampled <- w.WI.t$SiteSampled
w.WI.t <- subset(w.WI.t, select = -c(ID:AroclorCongener))
w.WI.t <- subset(w.WI.t, select = -c(AroclorA1016:AroclorA1260))
# Get mean congener per site, excluding zeros
tPCB <- rowSums(w.WI.t, na.rm = TRUE)
tPCB <- data.frame(cbind(site.sampled, tPCB))
tPCB$tPCB <- as.numeric(as.character(tPCB$tPCB))
tPCB.mean <- aggregate(tPCB ~ site.sampled, data = tPCB, mean)
# add coordinates
tPCB.mean <- data.frame(c(tPCB.mean, wi.coord))

# (3) Plot map + tPCB
ggmap(wi.map) +
  geom_point(data = tPCB.mean, aes(x = Long, y = Lat,
                              size = tPCB), alpha = 0.5) +
  scale_size_area(breaks = c(250, 500, 750, 1000, 1500),
                  labels = c(250, 500, 750, 1000, 1500),
                  name = "PCBs ng/L") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(title = "Fox River PCBs water concentration (pg/L) 2010-2018")

# Congener maps
# Select congener and remove samples with = 0 and NA for selected congener
w.WI.2 <- subset(w.WI, w.WI$PCB1 != 0 & w.WI$PCB1 != "NA")
# Get mean congener per site, excluding zeros
PCB1 <- aggregate(PCB1 ~ SiteSampled, data = w.WI.2, mean)
PCB1 <- data.frame(c(PCB1, wi.coord))

# (4) Plot map + congener
ggmap(wi.map) +
  geom_point(data = PCB1, aes(x = Long, y = Lat,
                              size = PCB1), alpha = 0.5) +
  scale_size_area(breaks = c(0.1, 1, 2, 4, 6),
                 labels = c(0.1, 1, 2, 4, 6),
                 name = "PCB 1 ng/L") +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(title = "Fox River PCB 1 water concentration (pg/L) 2010-2018")
  #geom_label_repel(data = PCB1, aes(x = Long, y = Lat, label = SiteSampled),
  #                 fill = "white", box.padding = unit(0.3, "lines"),
  #                 label.padding = unit(0.15, "lines"),
  #                 segment.color = "black", segment.size = 1)
                   
                   
