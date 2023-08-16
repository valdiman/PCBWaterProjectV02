library(cowplot)
library(grid)
library(png)
library(gridExtra)

# Read the first plot from a folder and convert to ggplot-compatible object
map <- rasterGrob(readPNG("Output/Maps/Sites/maptPCBFoxRiverAveV01.png"))

# Read the second plot from another folder and convert to ggplot-compatible object
plot <- rasterGrob(readPNG("Output/Plots/Sites/Temporal/plotFoxRiverTimeV01.png"))

# Create a grid arrangement with the map and plot positioned next to each other
grid_arrange <- arrangeGrob(
  map, plot,
  layout_matrix = matrix(c(1, 2), ncol = 2, byrow = TRUE),
  widths = unit(c(0.7, 0.3), c("npc", "npc")),
  heights = unit(c(1), "null")
)

# Display or save the grid arrangement
grid.arrange(grid_arrange)
