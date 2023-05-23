library(shiny)
library(ggplot2)
library(caper)
library(dplyr)

# Define the user interface
ui <- fluidPage(
  titlePanel("PGLS Mammal Brain Size Prediction (Beta)"),
  tags$small("Patience! This is a large dataset."),  # Add the subheading using tags$small
  sidebarLayout(
    sidebarPanel(
      selectInput("select_order", "\U0001F984 select Order:",
                  choices = c("all mammals", "Afrosoricida", "Artiodactyla", "Carnivora",  "Cetacea",  "Chiroptera","Cingulata","Dasyuromorphia","Diprotodontia", "Erinaceomorpha","Hyracoidea", "Lagomorpha", "Macroscelidea", "Monotremata","Notoryctemorphia",
                              "Paucituberculata", "Peramelemorphia","Perissodactyla", "Pholidota", "Pilosa",  "Primates",  "Proboscidea",  "Rodentia", "Scandentia",  "Sirenia",  "Soricomorpha","Tubulidentata")),
      textInput("body_size_text",
                "enter Body Size (g)"),
      textInput("brain_size_text",
                "to calculate Residual Index, enter Brain Size (g) if known \U0001F534"),
      checkboxInput("lambda_checkbox", "estimate phylogenetic signal (λ) using maximum likelihood (takes up to 8 min.)", value = FALSE),
      actionButton("my_prediction", "\U0001F535 Predict Brain Size ")
    ),
    mainPanel(
      plotOutput("plot"),
      verbatimTextOutput("result"),
      uiOutput("citation")
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
    if (input$select_order == "all mammals") {
      MammalData
    } else {
      filter(MammalData, order %in% input$select_order)
    }
  })

  # Create the comparative data object
  MammalOrder <- reactive({
    comparative.data(phy = MammalTree, data = OrderData(), names.col = Binomial, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
  })

  # Define the PGLS model
  model.pgls <- reactive({
    if (input$lambda_checkbox) {
      pgls(log(Mean_brain_mass_g) ~ log(Mean_body_mass_g), data = MammalOrder(), lambda = "ML")
    } else {
      pgls(log(Mean_brain_mass_g) ~ log(Mean_body_mass_g), data = MammalOrder())
    }
  })

  # Predict the brain size from the PGLS model based on the user's input
  predicted_brain <- reactive({
    exp(predict(model.pgls(), data.frame(Mean_body_mass_g = as.numeric(input$body_size_text))))
  })

  # Determine the brain size Residual Index from the PGLS model based on the user's input
  residual_brain <- reactive({
    logobserved <- log(as.numeric(input$brain_size_text))  # LOG Observed brain size
    logpredicted <- log(predicted_brain())  # LOG Predicted brain size
    residual_brain <- (logobserved - logpredicted) / logpredicted # Residual Index calculation
  })

  # Display the prediction result
  output$result <- renderText({
    if (input$lambda_checkbox) {
      lambda_estimate <- model.pgls()$param["lambda"][1]
      lambda_text <- paste("Phylogenetic signal estimated by maximum likelihood, lambda =", lambda_estimate)
    } else {
      lambda_text <- paste("Brownian motion assumed, fixed to lambda = 1")
    }
    
    prediction_text <- paste("Based on a body size of", input$body_size_text, "g, for", input$select_order, "the predicted brain size is", round(predicted_brain(), 2), "g \U0001F535.")
    
    if (!is.null(input$brain_size_text) && input$brain_size_text != "") {
      residual_text <- paste("It has a brain Residual Index value of", round(residual_brain(), 2), "based on an", input$select_order, "regression.")
    } else {
      residual_text <- ""
    }
    
    paste(prediction_text, "\n", residual_text, "\n", lambda_text)
  })
  
  # Generate the plot
  output$plot <- renderPlot({
    req(input$select_order, input$body_size_text, input$my_prediction)

    # Get the slope and intercept of the PGLS regression line
    slope <- coef(model.pgls())[2]
    intercept <- coef(model.pgls())[1]

    ggplot(OrderData(), aes(x = log(Mean_body_mass_g), y = log(Mean_brain_mass_g))) +
      geom_point() +
      geom_abline(intercept = intercept, slope = slope, color = "blue") +
      labs(x = "log Body Size (g)", y = "log Brain Size (g)") +
      geom_point(aes(x = log(Mean_body_mass_g), y = log(Mean_brain_mass_g)), color = "black", size = 3) +
      geom_point(aes(x = log(as.numeric(input$body_size_text)), y = log(as.numeric(input$brain_size_text))), color = "red", size = 3) +
      geom_point(aes(x = log(as.numeric(input$body_size_text)), y = log(predicted_brain())), color = "blue", size = 3) +
      theme(plot.margin = margin(10, 10, 30, 10))
  })

  # Display the citation
  output$citation <- renderText({
    citation_text <- "Data: Burger et al., 2019; Tree: Kumar et al., 2022; See"
    citation_link <- "https://github.com/AleAliSousa/ShinierBrains"
    paste(citation_text, "<a href='", citation_link, "'>", citation_link, "</a>")
  })
}

# Run the app
shinyApp(ui = ui, server = server)



