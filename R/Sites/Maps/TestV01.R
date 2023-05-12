
# Install packages
install.packages("leaflet")
install.packages("shiny")
install.packages("sp")

# Load libraries
{
  library(dplyr)
  library(shiny)
  library(leaflet)
  library(sp)
  library(stringr) # str_detect
  library(ggplot2)
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor08052022.csv")

# Select Chesapeake Bay & Delaware Canal data ---------------------------------------------------
che.0 <- wdc[str_detect(wdc$LocationName, 'Chesapeake Bay'),]
fox.0 <- wdc[str_detect(wdc$LocationName, 'Fox River'),]

# Data preparation --------------------------------------------------------
{
  # Remove samples (rows) with total PCBs  = 0
  che.1 <- che.0[!(rowSums(che.0[, c(14:117)], na.rm = TRUE)==0),]
  # Calculate total PCB
  tpcb.che <- rowSums(che.1[, c(14:117)], na.rm = T)
  # Create data.frame
  che.tpcb <- data.frame(
    SiteID = with(che.1, SiteID),
    date = with(che.1, as.Date(SampleDate, format = "%m/%d/%y")),
    Latitude = with(che.1, as.numeric(Latitude)),
    Longitude = with(che.1, as.numeric(Longitude)),
    tPCB = as.numeric(tpcb.che))
  colnames(che.tpcb) <- c("SiteID", "date", "Latitude", "Longitude", "tPCB")
}

{
  # Remove samples (rows) with total PCBs  = 0
  fox.1 <- fox.0[!(rowSums(fox.0[, c(14:117)], na.rm = TRUE)==0),]
  # Calculate total PCB
  tpcb.fox <- rowSums(fox.1[, c(14:117)], na.rm = T)
  # Create data.frame
  fox.tpcb <- data.frame(
    SiteID = with(fox.1, SiteID),
    date = with(fox.1, as.Date(SampleDate, format = "%m/%d/%y")),
    Latitude = with(fox.1, as.numeric(Latitude)),
    Longitude = with(fox.1, as.numeric(Longitude)),
    tPCB = as.numeric(tpcb.fox))
  colnames(fox.tpcb) <- c("SiteID", "date", "Latitude", "Longitude", "tPCB")
}


# (1)
ui <- fluidPage(
  leafletOutput("map"),
  tableOutput("data")
)

server <- function(input, output, session) {
  output$map <- renderLeaflet({
    leaflet(che.tpcb) %>%
      addTiles() %>%
      setView(lng = mean(che.tpcb$Longitude), lat = mean(che.tpcb$Latitude), zoom = 8)
  })
  
  output$data <- renderTable({
    if (!is.null(input$map_marker_click)) {
      siteid <- input$map_marker_click$id
      data <- subset(che.tpcb, SiteID == siteid)[, c("SiteID", "date", "tPCB")]
      data$date <- format(as.Date(data$date), "%m-%d-%Y")
      colnames(data)[3] <- paste0("\u03A3", "PCB (ng/L)")
      return(data)
    } else {
      return(data.frame())
    }
  })
  
  
  observe({
    leafletProxy("map", data = che.tpcb) %>%
      clearMarkers() %>%
      addMarkers(lng = ~Longitude, lat = ~Latitude,
                 layerId = ~SiteID,
                 popup = ~paste("SiteID: ", SiteID, "<br>",
                                "Latitude: ", Latitude, "<br>",
                                "Longitude: ", Longitude))
  })
  
}

shinyApp(ui, server)


# (2)
ui <- fluidPage(
  leafletOutput("map"),
  splitLayout(tableOutput("data"), plotOutput("plot", height = "300px"), 
              cellWidths = c("70%", "30%"))
)

server <- function(input, output, session) {
  output$map <- renderLeaflet({
    leaflet(che.tpcb) %>%
      addTiles() %>%
      setView(lng = mean(che.tpcb$Longitude), lat = mean(che.tpcb$Latitude), zoom = 8)
  })
  
  output$data <- renderTable({
    if (!is.null(input$map_marker_click)) {
      siteid <- input$map_marker_click$id
      data <- subset(che.tpcb, SiteID == siteid)[, c("SiteID", "date", "tPCB")]
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
      data <- subset(che.tpcb, SiteID == siteid)[, c("date", "tPCB")]
      data$date <- as.Date(data$date)
      p <- ggplot(data, aes(x = date, y = tPCB)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        labs(x = NULL, y = paste("\u03A3", "PCB (ng/L)", sep = "")) +
        scale_x_date(date_labels = "%b-%Y") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              axis.text = element_text(size = 12),
              plot.margin = margin(l = 10))
      return(p)
    } else {
      return(NULL)
    }
  })
  

  
  observe({
    leafletProxy("map", data = che.tpcb) %>%
      clearMarkers() %>%
      addMarkers(lng = ~Longitude, lat = ~Latitude,
                 layerId = ~SiteID,
                 popup = ~paste("SiteID: ", SiteID, "<br>",
                                "Latitude: ", Latitude, "<br>",
                                "Longitude: ", Longitude))
  })
}

shinyApp(ui, server)

# (3)
ui <- fluidPage(
  selectInput("dataset", label = "Select dataset:", choices = c("Chesapeake Bay" = "che.tpcb", "Fox River" = "fox.tpcb")),
  leafletOutput("map"),
  splitLayout(tableOutput("data"), plotOutput("plot", height = "300px"), 
              cellWidths = c("70%", "30%"))
)

server <- function(input, output, session) {
  
  data <- reactive({
    if(input$dataset == "che.tpcb") {
      che.tpcb
    } else {
      fox.tpcb
    }
  })
  
  output$map <- renderLeaflet({
    leaflet(data()) %>%
      addTiles() %>%
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
      p <- ggplot(data, aes(x = date, y = tPCB)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        labs(x = NULL, y = paste("\u03A3", "PCB (ng/L)", sep = "")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              axis.text = element_text(size = 12),
              plot.margin = margin(l = 10))
      return(p)
    } else {
      return(NULL)
    }
  })
  
  observe({
    leafletProxy("map", data = data()) %>%
      clearMarkers() %>%
      addMarkers(lng = ~Longitude, lat = ~Latitude,
                 layerId = ~SiteID,
                 popup = ~paste("SiteID: ", SiteID, "<br>",
                                "Latitude: ", Latitude, "<br>",
                                "Longitude: ", Longitude))
  })
  
}

shinyApp(ui, server)

