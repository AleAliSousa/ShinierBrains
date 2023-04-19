library(shiny)
library(ggplot2)
library(caper)
library(dplyr)

# Define the user interface
ui <- fluidPage(
  titlePanel("PGLS Mammal Brain Size Prediction (Beta)"),
  sidebarLayout(
    sidebarPanel(
      selectInput("select_order", "Select Order:",
                  choices = c("Afrosoricida", "Artiodactyla", "Carnivora",  "Cetacea",  "Chiroptera","Cingulata","Dasyuromorphia","Diprotodontia", "Erinaceomorpha","Hyracoidea", "Lagomorpha", "Macroscelidea", "Monotremata","Notoryctemorphia", 
                              "Paucituberculata", "Peramelemorphia","Perissodactyla", "Pholidota", "Pilosa",  "Primates",  "Proboscidea",  "Rodentia", "Scandentia",  "Sirenia",  "Soricomorpha","Tubulidentata")),
      textInput("body_size_text",
                "enter body size (g)"),
      textInput("brain_size_text",
                "if known, enter brain size (g) (appears as red point)"),
      #actionButton("my_prediction", "Predict Brain Size")
    ),
    mainPanel(
      plotOutput("plot"),
      verbatimTextOutput("result"),
      textOutput("citation")
    )
  )
)
# Define the server logic
server <- function(input, output) {
  
  # Load the data
  MammalData <- read.csv("gyz043_suppl_Supplement_Data.csv")
  MammalData$log_brain_mass_g <- log(MammalData$Mean_brain_mass_g)
  MammalData$log_body_mass_g <- log(MammalData$Mean_body_mass_g)
  
  # Load the tree
  MammalTree <- read.tree("species.nwk")
  
  # Create a reactive object to filter the data based on the user's selected order 
  OrderData <- reactive({
    filter(MammalData, order %in% input$select_order)
  })
  
  # Create the comparative data object
  MammalOrder <- reactive({
    comparative.data(phy = MammalTree, data = OrderData(), names.col = Binomial, vcv = TRUE, na.omit = FALSE, warn.dropped = FALSE)
  })
  
  # Define the PGLS model
  model.pgls<-reactive({
    pgls(log(Mean_brain_mass_g) ~ log(Mean_body_mass_g), data = MammalOrder())
  })
  
  # Predict the brain size from the PGLS model based on the user's input
  predicted_brain <- reactive({
    exp(predict(model.pgls(), data.frame(Mean_body_mass_g = as.numeric(input$body_size_text))))
  })
  
  # Display the prediction result
  output$result <- renderText({
    paste("Based on a body size of", input$body_size_text, "g, the predicted brain size for", input$select_order, "is", round(predicted_brain(), 2), "g (blue point).")
  })
  
  # Generate the plot
  output$plot <- renderPlot({
    req(input$select_order, input$body_size_text)
    
    predicted_brain_body_size <- exp(predict(model.pgls(), newdata = data.frame(Mean_body_mass_g = as.numeric(input$body_size_text))))
    
    # Get the slope and intercept of the PGLS regression line
    slope <- coef(model.pgls())[2]
    intercept <- coef(model.pgls())[1]
    
    ggplot(OrderData(), aes(x = log(Mean_body_mass_g), y = log(Mean_brain_mass_g))) +
      geom_point() +
      geom_abline(intercept = intercept, slope = slope, color = "blue") +
      labs(x = "Body Size (log g)", y = "Brain Size (log g)") +
      geom_point(aes(x = log(Mean_body_mass_g), y = log(Mean_brain_mass_g)), color = "black", size = 3) +
      geom_point(aes(x = log(as.numeric(input$body_size_text)), y = log(as.numeric(input$brain_size_text))), color = "red", size = 3) +
      geom_point(aes(x = log(as.numeric(input$body_size_text)), y = log(predicted_brain_body_size)), color = "blue", size = 3) +
      theme(plot.margin = margin(10, 10, 30, 10))
  })
  
  # Display the citation
  output$citation <- renderText({
      "Data: Burger et al., 2019; Tree: Kumar et al., 2022; See https://github.com/AleAliSousa/ShinierBrains"
  })
}

  # Run the app
  shinyApp(ui = ui, server = server)

