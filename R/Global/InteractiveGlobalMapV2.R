
# Install packages
install.packages("dplyr")
install.packages("ggplot2")
install.packages("leaflet")
install.packages("shiny")
install.packages("tidyverse")

# Load libraries
{
  library(dplyr)
  library(ggplot2)
  library(leaflet)
  library(shiny)
  library(stringr) # str_detect
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor08052022.csv")

# Prepare data ------------------------------------------------------------
datasets <- list()

# Create datasets based on LocationName
unique_locations <- unique(trimws(wdc$LocationName))
for (location in unique_locations) {
  dataset_name <- gsub(" ", "_", location)  # Generate dataset name based on LocationName
  
  datasets[[dataset_name]] <- wdc[trimws(wdc$LocationName) == location, ]
}

# Remove the "All_Locations" dataset from the list
datasets[["All_Locations"]] <- NULL

# Process and combine the datasets
combinedData <- NULL  # Initialize an empty data frame

for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  
  # Remove samples (rows) with total PCBs = 0
  dataset <- dataset[!(rowSums(dataset[, 14:117], na.rm = TRUE) == 0), ]
  
  # Calculate total PCB
  tpcb <- rowSums(dataset[, 14:117], na.rm = TRUE)
  
  # Create data.frame
  processed_data <- data.frame(
    SampleID = with(dataset, SampleID),
    date = with(dataset, as.Date(SampleDate, format = "%m/%d/%y")),
    Latitude = with(dataset, as.numeric(Latitude)),
    Longitude = with(dataset, as.numeric(Longitude)),
    tPCB = as.numeric(tpcb),
    dataset = dataset_name
  )
  
  colnames(processed_data) <- c("SampleID", "date", "Latitude",
                                "Longitude", "tPCB", "dataset")
  
  # Append the processed data to the combinedData
  combinedData <- rbind(combinedData, processed_data)
}

# Remove duplicated rows based on SampleID column
combinedData <- distinct(combinedData, SampleID, .keep_all = TRUE)

# Shiny app ---------------------------------------------------------------
# Define the Shiny UI
ui <- fluidPage(
  titlePanel("PCB Water Concentration Data Visualization"),
  selectInput("dataset", "Select Dataset:", choices = c("All", names(datasets))),
  leafletOutput("map"),
  splitLayout(
    tableOutput("data"),
    plotOutput("plot", height = "300px", width = "700px"),
    cellWidths = c("35%", "65%")
  ),
  verbatimTextOutput("plot_text")
)

# Define the Shiny server
server <- function(input, output, session) {
  # Subset data based on selected dataset
  filtered_data <- reactive({
    if (input$dataset == "All") {
      return(combinedData)
    } else {
      return(combinedData[combinedData$dataset == input$dataset, ])
    }
  })
  
  # Render the map
  output$map <- renderLeaflet({
    leaflet(filtered_data()) %>%
      addTiles() %>%
      addMarkers(
        lng = ~Longitude,
        lat = ~Latitude,
        label = ~as.character(SampleID),
        labelOptions = labelOptions(noHide = TRUE)
      )
  })
  
  # Render the table
  output$data <- renderTable({
    filtered_data <- filtered_data()
    
    if (!is.null(input$map_marker_click)) {
      sampleid <- input$map_marker_click$id
      subset_data <- subset(filtered_data, SampleID == sampleid)[, c("SampleID", "date", "tPCB")]
      subset_data$date <- format(as.Date(subset_data$date, format = "%m/%d/%y"), "%m-%d-%Y")
      colnames(subset_data)[3] <- paste("\u03A3", "PCB ", "(ng/L)", sep = "")
      
      # Sort the data by date
      subset_data <- subset_data[order(as.Date(subset_data$date, format = "%m-%d-%Y")), ]
      
      # Get the number of samples
      num_samples <- nrow(subset_data)
      
      colnames(subset_data)[1] <- paste("SampleID (n =", num_samples, ")")
      
      return(subset_data)
    } else {
      return(data.frame())
    }
  })
  
  output$plot <- renderPlot({
    filtered_data <- filtered_data()
    
    if (!is.null(input$map_marker_click)) {
      sampleid <- input$map_marker_click$id
      subset_data <- subset(filtered_data, SampleID == sampleid)
      
      if (nrow(subset_data) == 0) {
        # No data available for the selected SampleID
        return(NULL)
      }
      
      # Aggregate data by week and calculate the average of PCB values
      subset_data$week <- cut(subset_data$date, breaks = "week")
      data_agg <- aggregate(tPCB ~ week, data = subset_data, mean)
      
      num_values <- nrow(data_agg)
      width <- ifelse(num_values > 5, 0.8, 0.2 + (num_values * 0.1))
      
      p <- ggplot(data_agg, aes(x = week, y = tPCB)) +
        geom_col(fill = "steelblue", width = width) +
        labs(x = NULL, y = paste("\u03A3", "PCB (ng/L)", sep = "")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size = 12),  # Increase the font size of y-axis numbers
          axis.title.y = element_text(size = 14),  # Increase the font size of y-axis title
          axis.text = element_text(size = 8)
        )
      
      # Check if y-values are too large for linear scale
      max_value <- max(data_agg$tPCB)
      if (max_value >= 50000) {
        p <- p + scale_y_log10()
      }
      
      return(p)
    } else {
      return(NULL)
    }
  })
  
  output$plot_text <- renderPrint({
    if (!is.null(input$map_marker_click)) {
      cat("Plots are showing the aggregated data per week.\n")
      cat("If the maximum tPCB is too large (>80000 ng/L), the y-axis changes to log10 scale.\n")
      cat("Source:")
    }
  })  
  
  observe({
    leafletProxy("map", data = filtered_data()) %>%
      clearMarkers() %>%
      addMarkers(
        lng = ~Longitude,
        lat = ~Latitude,
        layerId = ~SampleID,
        popup = ~paste(
          "SampleID: ", SampleID, "<br>",
          "Latitude: ", Latitude, "<br>",
          "Longitude: ", Longitude, "<br>",
          "Dataset: ", dataset
        )
      )
  })
}

# Run the Shiny app
shinyApp(ui, server)
