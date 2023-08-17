
library(ggplot2)
library(png)

# Read the PNG images
map <- readPNG("Output/Maps/Sites/maptPCBFoxRiverAveV01.png")
plot <- readPNG("Output/Plots/Sites/Temporal/plotFoxRiverTimeV01.png")

# Create a plot with a smaller canvas
p <- ggplot() +
  xlim(0, 6) + ylim(0, 6) +  # Adjust the limits as needed
  theme_void()

# Add the images at specific positions with reduced size
p <- p +
  annotation_custom(rasterGrob(map),
                    xmin = 0, xmax = 5, ymin = 1, ymax = 6) +
  annotation_custom(rasterGrob(plot, width = unit(1.5*.8, "in"), height = unit(0.75*.8, "in")),
                    xmin = 3.3, xmax = 5.3, ymin = 1.8, ymax = 3.8)

# Print the plot
print(p)

# Save plot in folder
ggsave("Output/Figures/Sites/FoxRiverV01.png",
       plot = p, width = 6, height = 4, dpi = 300)
