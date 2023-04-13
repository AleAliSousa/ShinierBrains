library(shiny)
library(ggplot2)
library(caper)
library(dplyr)

# Load the data
MammalData <- read.csv("gyz043_suppl_Supplement_Data.csv")
MammalData$log_brain_mass_g <- log(MammalData$Mean_brain_mass_g)
MammalData$log_body_mass_g <- log(MammalData$Mean_body_mass_g)

# Load the tree
MammalTree <- read.tree("species.nwk")

# comparative data object
Mammals <- comparative.data(phy = MammalTree, data = MammalData, names.col = Binomial, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

# Define the user interface
ui <- fluidPage(
  titlePanel("PGLS Mammal Brain Size Prediction"),
  sidebarLayout(
    sidebarPanel(
      textInput("body_size_text",
                "enter body size (g)"),
      textInput("brain_size_text",
                "if known, enter brain size (g)"),
      actionButton("predict", "Predict Brain Size")
    ),
    mainPanel(
      plotOutput("plot"),
      verbatimTextOutput("result")
    )
  )
)

# Define the server logic
server <- function(input, output) {

  # Define the PGLS model
  model.pgls<-pgls(log(Mean_brain_mass_g) ~ log(Mean_body_mass_g), data = Mammals)

  # Generate the plot
  output$plot <- renderPlot({
    new_data <- data.frame(Mean_body_mass_g = as.numeric(input$body_size_text))
    predicted_brain_body_size <- exp(predict(model.pgls, newdata = new_data))

    # Get the slope and intercept of the PGLS regression line
    slope <- coef(model.pgls)[2]
    intercept <- coef(model.pgls)[1]
    
    ggplot(MammalData, aes(x = log(Mean_body_mass_g), y = log(Mean_brain_mass_g))) +
      geom_point() +
      geom_abline(intercept = intercept, slope = slope, color = "blue") +
      labs(x = "Body Size (log g)", y = "Brain Size (log g)") +
      geom_point(aes(x = log(Mean_body_mass_g), y = log(Mean_brain_mass_g)), color = "black", size = 3) +
      geom_point(aes(x = log(as.numeric(input$body_size_text)), y = log(as.numeric(input$brain_size_text))), color = "red", size = 3) +
      geom_point(aes(x = log(as.numeric(input$body_size_text)), y = log(predicted_brain_body_size)), color = "blue", size = 3) +
      theme(plot.margin = margin(10, 10, 30, 10))
  })

  # Predict the brain size based on the user's input
  predicted_brain <- eventReactive(input$predict, {
    new_data <- data.frame(Mean_body_mass_g = as.numeric(input$body_size_text))
    exp(predict(model.pgls, newdata = new_data))
  })

    # Display the prediction result
    output$result <- renderText({
      paste("Based on a body size of", input$body_size_text, "g, the predicted brain size is", round(predicted_brain(), 2), "g.")
    })

  }

  # Run the app
  shinyApp(ui = ui, server = server)
