
# Install packages
install.packages("leaflet")
install.packages("shiny")


# Load libraries
{
  library(dplyr)
  library(shiny)
  library(leaflet)
  library(stringr) # str_detect
  library(ggplot2)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor08052022.csv")

datasets <- list(
  "Chesapeake Bay" = wdc[str_detect(wdc$LocationName, 'Chesapeake Bay'),],
  "Fox River" = wdc[str_detect(wdc$LocationName, 'Fox River'),],
  "Housatonic River" = wdc[str_detect(wdc$LocationName, 'Housatonic River'),],
  "Hudson River" = wdc[str_detect(wdc$LocationName, 'Hudson River'),],
  "Kalamazoo River" = wdc[str_detect(wdc$LocationName, 'Kalamazoo River'),],
  "New Bedford" = wdc[str_detect(wdc$LocationName, 'New Bedford'),],
  "Portland Harbor" = wdc[str_detect(wdc$LocationName, 'Portland Harbor'),],
  "Spokane River" = wdc[str_detect(wdc$LocationName, 'Spokane River'),]
)

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
    SiteID = with(dataset, SiteID),
    date = with(dataset, as.Date(SampleDate, format = "%m/%d/%y")),
    Latitude = with(dataset, as.numeric(Latitude)),
    Longitude = with(dataset, as.numeric(Longitude)),
    tPCB = as.numeric(tpcb),
    dataset = dataset_name
  )
  
  colnames(processed_data) <- c("SiteID", "date", "Latitude", "Longitude", "tPCB", "dataset")
  
  # Append the processed data to the combinedData
  combinedData <- rbind(combinedData, processed_data)
}

# Define the Shiny UI
ui <- fluidPage(
  titlePanel("PCB Water Concentration Data Visualization"),
  selectInput("dataset", label = "Select dataset:", choices = unique(combinedData$dataset)),
  leafletOutput("map"),
  splitLayout(
    tableOutput("data"),
    plotOutput("plot", height = "300px"),
    cellWidths = c("70%", "30%")
  )
)

# Define the Shiny server
server <- function(input, output, session) {
  
  data <- reactive({
    filter(combinedData, dataset == input$dataset)
  })
  
  output$map <- renderLeaflet({
    leaflet(data()) %>%
      addTiles() %>%
      fitBounds(lng1 = min(data()$Longitude),
                lat1 = min(data()$Latitude),
                lng2 = max(data()$Longitude),
                lat2 = max(data()$Latitude)) %>%
      setView(lng = mean(data()$Longitude), lat = mean(data()$Latitude), zoom = 8)
  })
  
  output$data <- renderTable({
    if (!is.null(input$map_marker_click)) {
      siteid <- input$map_marker_click$id
      data <- subset(data(), SiteID == siteid)[, c("SiteID", "date", "tPCB")]
      data$date <- format(as.Date(data$date), "%m-%d-%Y")
      colnames(data)[3] <- paste("\u03A3", "PCB ", "(ng/L)", sep = "")
      return(data)
    } else {
      return(data.frame())
    }
  })
  
  output$plot <- renderPlot({
    if (!is.null(input$map_marker_click)) {
      siteid <- input$map_marker_click$id
      data <- subset(data(), SiteID == siteid)[, c("date", "tPCB")]
      data$date <- as.Date(data$date)
      
      if (nrow(data) == 0) {
        # No data available for aggregation
        return(NULL)
      }
      
      # Aggregate data by week and calculate the average of PCB values
      data$week <- cut(data$date, breaks = "week")
      data_agg <- aggregate(tPCB ~ week, data, mean)
      
      p <- ggplot(data_agg, aes(x = week, y = tPCB)) +
        geom_col(fill = "steelblue", width = 0.8) +
        labs(x = NULL, y = paste("\u03A3", "PCB (ng/L)", sep = "")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text = element_text(size = 12)
        )
      return(p)
    } else {
      return(NULL)
    }
  })
  
  
  observe({
    leafletProxy("map", data = data()) %>%
      clearMarkers() %>%
      addMarkers(
        lng = ~Longitude,
        lat = ~Latitude,
        layerId = ~SiteID,
        popup = ~paste(
          "SiteID: ", SiteID, "<br>",
          "Latitude: ", Latitude, "<br>",
          "Longitude: ", Longitude
        )
      )
  })
  
}

shinyApp(ui, server)
