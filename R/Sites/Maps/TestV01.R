
# Install packages
install.packages("leaflet")
install.packages("shiny")
install.packages("sp")
install.packages("plotly")

# Load libraries
{
  library(dplyr)
  library(shiny)
  library(leaflet)
  library(sp)
  
}

# Read data ---------------------------------------------------------------
# Data in pg/L
wdc <- read.csv("Data/WaterDataCongenerAroclor08052022.csv")

# Select Chesapeake Bay & Delaware Canal data ---------------------------------------------------
che.0 <- wdc[str_detect(wdc$LocationName, 'Chesapeake Bay'),]

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
      data <- subset(che.tpcb, SiteID == siteid)[, c("SiteID", "date", "Latitude", "Longitude", "tPCB")]
      data$date <- format(as.Date(data$date), "%m-%d-%Y")
      colnames(data)[5] <- paste0("\u03A3", "PCB (ng/L)")
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
                                "Date: ", format(as.Date(date), "%m-%d-%Y"), "<br>",
                                "tPCB: ", tPCB, " ppb", "<br>",
                                tableOutput("data")))
  })
  
}

shinyApp(ui, server)


