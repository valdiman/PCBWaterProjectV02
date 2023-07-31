## Water PCB concentrations mapping.

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
install.packages("ggrepel")
install.packages("ggpp")
install.packages("scales")

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
  library(scales) # add commas in legend in maps
  library(cowplot)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor08052022.csv")
# Extract sample site locations -------------------------------------------
# Calculate total PCB
# Data preparation
{
  # Remove samples (rows) with total PCBs  = 0
  wdc.1 <- wdc[!(rowSums(wdc[, c(14:117)], na.rm = TRUE)==0),]
  # Calculate total PCB
  tpcb <- rowSums(wdc.1[, c(14:117)], na.rm = T)
  # Select and combine sample sites anf tPCB
  location <- cbind.data.frame(wdc.1$SiteID, wdc.1$Latitude, wdc.1$Longitude,
                               tpcb)
  # Name the columns
  colnames(location) <- c("SiteID", "Latitude", "Longitude", "tPCB.ave")
  # Average tPCB per sample site
  tpcb.ave.loc <- aggregate(tPCB.ave ~ SiteID + Latitude + Longitude,
                           data = location, mean)
}

# USA/State maps -------------------------------------------------------------
us <- map_data("usa")
states <- map_data("state")

# Find number of samples per state to be included as table in maps
{
  wdc.2 <- wdc %>%
    group_by(StateSampled) %>%
    summarise(n = n())
  #colnames(wdc.1) <- c("State", "# samples")
  wdc.3 <- data.frame(t(wdc.2))
  name <- c('State', '# samples')
  wdc.3 <- data.frame(col1 = name, wdc.3)
  names(wdc.3) <- NULL
}

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
  geom_point(data = wdc, aes(x = Longitude, y = Latitude),
             color = "black",
             size = 1.2, shape = 20) +
  annotate(geom = 'table', x = -65, y = 53,
           label = list(wdc.3), size = 2.9) # add table with info

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
  geom_point(data = tpcb.ave.loc, aes(x = Longitude, y = Latitude,
                                       size = tPCB.ave), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  theme(legend.position = "right") +
  scale_size_area(breaks = c(1000, 50*1000, 500*1000, 1000*1000, 1500*1000,
                             2000*1000), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (SiteID mean)",
                                              paste("1990-2020 (pg/L)")))),
                  max_size = 5)

# Specific locations ------------------------------------------------------
# Portland Harbor ---------------------------------------------------------
{
  # Select only from Portland Harbor
  wdc.PO <- subset(wdc, LocationName == "Portland Harbor")
  # Create general map
  PO.box <- make_bbox(lon = wdc.PO$Longitude, lat = wdc.PO$Latitude, f = 0.8)
  PO.map <- get_stamenmap(bbox = PO.box, zoom = 10)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.PO <- data.frame(cbind(wdc.PO$SiteID, wdc.PO$Latitude,
                              wdc.PO$Longitude, rowSums(wdc.PO[, c(14:117)],
                                                        na.rm = TRUE)))
  # Name the columns
  colnames(tPCB.PO) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.PO$Latitude <- as.numeric(tPCB.PO$Latitude)
  tPCB.PO$Longitude <- as.numeric(tPCB.PO$Longitude)
  tPCB.PO$tPCB <- as.numeric(tPCB.PO$tPCB)
  # Average tPCB per site
  tPCB.PO.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                           data = tPCB.PO, FUN = mean)
}

# (1) Plot map + locations
ggmap(PO.map) +
  geom_point(data = tPCB.PO.ave, aes(x = Longitude, y = Latitude),
             shape = 21, color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.PO.ave, family = 'Times New Roman', size = 3, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(PO.map) +
  geom_point(data = tPCB.PO.ave, aes(x = Longitude, y = Latitude,
                                     size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_size_area(breaks = c(100, 250, 500, 1000, 1500), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 2007-2019 (pg/L)"))),
                  max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(1.6, 0.8),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

# Fox River ---------------------------------------------------------------
{
  # Select only from Fox River
  wdc.Fox <- subset(wdc, LocationName == "Fox River")
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
  tPCB.Fox <- data.frame(cbind(wdc.Fox.1$SiteID, wdc.Fox.1$Latitude,
                               wdc.Fox.1$Longitude,
                               rowSums(wdc.Fox.1[, c(14:117)], na.rm = TRUE)))
  # Name the columns
  colnames(tPCB.Fox) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.Fox$Latitude <- as.numeric(tPCB.Fox$Latitude)
  tPCB.Fox$Longitude <- as.numeric(tPCB.Fox$Longitude)
  tPCB.Fox$tPCB <- as.numeric(tPCB.Fox$tPCB)
  # Average tPCB per site
  tPCB.Fox.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                             data = tPCB.Fox, FUN = mean)
}

# (1) Plot map + locations
ggmap(Fox.map) +
  geom_point(data = tPCB.Fox.ave, aes(x = Longitude, y = Latitude),
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
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.Fox.ave, family = 'Times New Roman', size = 2.8, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(Fox.map) +
  geom_point(data = tPCB.Fox.ave, aes(x = Longitude, y = Latitude,
                                      size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_size_area(breaks = c(100, 500, 1500, 30000, 65000), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 2005-2018 (pg/L)"))),
                  max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(1.35, 0.8),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

# Hudson River ------------------------------------------------------------
{
  # Select only from Hudson River
  wdc.Hud <- subset(wdc, LocationName == "Hudson River")
  # Create 2 locations to increase with of the map
  #l <- data.frame(SiteID = "WCPCB-HUD015", Latitude = 43.29781, Longitude = -75, tPCB = 0)
  #tPCB.Hud.ave <- rbind(tPCB.Hud.ave, l)
  
  # Increase the longitude range to make the map wider
  lon_range <- 0.5  # Modify this value to control the width
  
  # Create a new bounding box with the adjusted longitude range
  Hud.box <- make_bbox(
    lon = c(min(wdc.Hud$Longitude) - lon_range, max(wdc.Hud$Longitude) + lon_range),
    lat = wdc.Hud$Latitude,
    f = 0.1
  )
  
  # Fetch the map using the new bounding box
  Hud.map <- get_stamenmap(bbox = Hud.box, zoom = 8)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.Hud <- data.frame(cbind(wdc.Hud$SiteID, wdc.Hud$Latitude,
                               wdc.Hud$Longitude, rowSums(wdc.Hud[, c(14:117)],
                                                          na.rm = TRUE)))
  # Name the columns
  colnames(tPCB.Hud) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.Hud$Latitude <- as.numeric(tPCB.Hud$Latitude)
  tPCB.Hud$Longitude <- as.numeric(tPCB.Hud$Longitude)
  tPCB.Hud$tPCB <- as.numeric(tPCB.Hud$tPCB)
  # Average tPCB per site
  tPCB.Hud.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                             data = tPCB.Hud, FUN = mean)
}

# (1) Plot map + locations
ggmap(Hud.map) +
  geom_point(data = tPCB.Hud.ave, aes(x = Longitude, y = Latitude),
             shape = 21, color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.Hud.ave, family = 'Times', size = 1.8, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# Plot the map with points and legend
ggmap(Hud.map) +
  geom_point(data = tPCB.Hud.ave, aes(x = Longitude, y = Latitude, size = tPCB),
             alpha = 1, color = "black", shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_size_area(
    breaks = c(2000, 5000, 10000, 20000, 30000, 40000),
    labels = comma,
    name = expression(bold(atop(Sigma*"PCBs (mean) 2005-2017 (pg/L)"))),
    max_size = 8
  ) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = "right",  # Move the legend to the bottom
        legend.justification = c(0.3, 1),  # Center the legend horizontally
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

# Housatonic River --------------------------------------------------------
{
  # Select only from Housotonic River
  wdc.Hou <- subset(wdc, LocationName == "Housatonic River")
  
  # Create general map
  Hou.box <- make_bbox(lon = wdc.Hou$Longitude, lat = wdc.Hou$Latitude, f = 0.2)
  Hou.map <- get_stamenmap(bbox = Hou.box, zoom = 10)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.Hou <- data.frame(cbind(wdc.Hou$SiteID, wdc.Hou$Latitude,
                               wdc.Hou$Longitude, rowSums(wdc.Hou[, c(14:117)],
                                                          na.rm = TRUE)))
  # Name the columns
  colnames(tPCB.Hou) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.Hou$Latitude <- as.numeric(tPCB.Hou$Latitude)
  tPCB.Hou$Longitude <- as.numeric(tPCB.Hou$Longitude)
  tPCB.Hou$tPCB <- as.numeric(tPCB.Hou$tPCB)
  # Average tPCB per site
  tPCB.Hou.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                            data = tPCB.Hou, FUN = mean)
}

# (1) Plot map + locations
ggmap(Hou.map) +
  geom_point(data = tPCB.Hou.ave, aes(x = Longitude, y = Latitude),
             shape = 21, color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.Hou.ave, family = 'Times New Roman', size = 3, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(Hou.map) +
  geom_point(data = tPCB.Hou.ave, aes(x = Longitude, y = Latitude,
                                      size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_size_area(breaks = c(50, 100, 500, 1000, 1500), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 1998-2020 (pg/L)"))),
                  max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(2.2, 0.8),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

# Kalamazoo River ---------------------------------------------------------
{
  # Select only from Kalamazoo  River
  wdc.Kal <- subset(wdc, LocationName == "Kalamazoo River")
  
  # Create general map
  Kal.box <- make_bbox(lon = wdc.Kal$Longitude, lat = wdc.Kal$Latitude, f = 0.2)
  Kal.map <- get_stamenmap(bbox = Kal.box, zoom = 10)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.Kal <- data.frame(cbind(wdc.Kal$SiteID, wdc.Kal$Latitude,
                               wdc.Kal$Longitude,
                               rowSums(wdc.Kal[, c(14:117)],
                                       na.rm = TRUE)))
  # Name the columns
  colnames(tPCB.Kal) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.Kal$Latitude <- as.numeric(tPCB.Kal$Latitude)
  tPCB.Kal$Longitude <- as.numeric(tPCB.Kal$Longitude)
  tPCB.Kal$tPCB <- as.numeric(tPCB.Kal$tPCB)
  # Average tPCB per site
  tPCB.Kal.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                         data = tPCB.Kal, FUN = mean)
}

# (1) Plot map + locations
ggmap(Kal.map) +
  geom_point(data = tPCB.Kal.ave, aes(x = Longitude, y = Latitude),
             shape = 21, color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.Kal.ave, family = 'Times New Roman', size = 3, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(Kal.map) +
  geom_point(data = tPCB.Kal.ave, aes(x = Longitude, y = Latitude,
                                      size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_size_area(breaks = c(100, 1000, 10000, 50000, 700000), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 1994-2010 (pg/L)"))),
                  max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(1.27, 0.78),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

# New Bedford -------------------------------------------------------------
{
  # Select only from Hudson River
  wdc.NB <- subset(wdc, LocationName == "New Bedford Harbor")
  
  # Create general map
  NB.box <- make_bbox(lon = wdc.NB$Longitude, lat = wdc.NB$Latitude, f = 0.2)
  NB.map <- get_stamenmap(bbox = NB.box, zoom = 12)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.NB <- data.frame(cbind(wdc.NB$SiteID, wdc.NB$Latitude, wdc.NB$Longitude,
                              rowSums(wdc.NB[, c(14:117)],
                                      na.rm = TRUE)))
  # Name the columns
  colnames(tPCB.NB) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.NB$Latitude <- as.numeric(tPCB.NB$Latitude)
  tPCB.NB$Longitude <- as.numeric(tPCB.NB$Longitude)
  tPCB.NB$tPCB <- as.numeric(tPCB.NB$tPCB)
  # Average tPCB per site
  tPCB.NB.ave <- aggregate(tPCB ~ SiteID+ Latitude + Longitude,
                         data = tPCB.NB, FUN = mean)
}

# (1) Plot map + locations
ggmap(NB.map) +
  geom_point(data = tPCB.NB.ave, aes(x = Longitude, y = Latitude), shape = 21,
             color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                 data = tPCB.NB.ave, family = 'Times', size = 1.8, 
                 box.padding = 0.2, point.padding = 0.3,
                 segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(NB.map) +
  geom_point(data = tPCB.NB.ave, aes(x = Longitude, y = Latitude,
                                      size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_size_area(breaks = c(10000, 500000, 1000000, 3000000), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 1994-2010 (pg/L)"))),
                  max_size = 8) +
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(4, 0.78),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

# Spokane River -----------------------------------------------------------
{
  # Select only from Hudson River
  wdc.Spo <- subset(wdc, LocationName == "Spokane River")
  
  # Create general map
  Spo.box <- make_bbox(lon = wdc.Spo$Longitude, lat = wdc.Spo$Latitude, f = 0.6)
  Spo.map <- get_stamenmap(bbox = Spo.box, zoom = 10)
  
  # Plot map with sites
  # Prepare data
  # Get tPCB and coordinates
  tPCB.Spo <- data.frame(cbind(wdc.Spo.1$SiteID, wdc.Spo.1$Latitude,
                               wdc.Spo.1$Longitude,
                               rowSums(wdc.Spo.1[, c(14:117)],
                                       na.rm = TRUE)))
  # Name the columns
  colnames(tPCB.Spo) <- c("SiteID", "Latitude", "Longitude", "tPCB")
  # Change no numeric to numeric
  tPCB.Spo$Latitude <- as.numeric(tPCB.Spo$Latitude)
  tPCB.Spo$Longitude <- as.numeric(tPCB.Spo$Longitude)
  tPCB.Spo$tPCB <- as.numeric(tPCB.Spo$tPCB)
  # Average tPCB per site
  tPCB.Spo.ave <- aggregate(tPCB ~ SiteID + Latitude + Longitude,
                         data = tPCB.Spo, FUN = mean)
}

# (1) Plot map + locations
ggmap(Spo.map) +
  geom_point(data = tPCB.Spo.ave, aes(x = Longitude, y = Latitude), shape = 21,
             color = "red",
             fill = "white", size = 1.75, stroke = 0.75) +
  geom_label_repel(aes(x = Longitude, y = Latitude, label = SiteID),
                   data = tPCB.Spo.ave, family = 'Times', size = 1.8, 
                   box.padding = 0.2, point.padding = 0.3,
                   segment.color = 'grey50')

# (2) Plot map + tPCB
ggmap(Spo.map) +
  geom_point(data = tPCB.Spo.ave, aes(x = Longitude, y = Latitude,
                                     size = tPCB), alpha = 1, color  = "black",
             shape = 21, fill = "white", stroke = 0.75) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_size_area(breaks = c(120, 1000, 3000, 5000, 8000), labels = comma,
                  name = expression(bold(atop(Sigma*"PCBs (mean) 2014-2016 (pg/L)"))),
                  max_size = 8)
  guides(size = guide_legend(label.hjust = 0.5)) +
  theme(legend.position = c(4, 0.78),  # Adjust the legend.position values
        legend.title = element_text(margin = margin(b = -16, unit = "pt")))

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
                   
                   
