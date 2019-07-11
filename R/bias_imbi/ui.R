library(shiny)
library(markdown)
library(plotly)

# Define UI for slider demo application
shinyUI(fluidPage(
     
     #  Application title
   titlePanel("drugdevelopR: bias"),
     # Sidebar with sliders that demonstrate various available
     # options  
  sidebarLayout(
    sidebarPanel(
       selectInput("Select", label = h3("Select adjustment method"), 
                   choices = list("Multiplicative adjustment " = 1,
                                  "Additive adjustment " = 2,
                                  "Multiplicative vs. additive adjustment " = 3), 
                   selected = 1),

          sliderInput("D2", HTML("Total sample size for phase II d<sub>2</sub>:"), 
                             min = 10, max = 750, value = c(250,350), step=2),
       sliderInput("stepD2", HTML("Step size for d<sub>2</sub>:"), 
                   min = 2, max = 50, value = 10, step=5),
          sliderInput("HRgo", HTML("Threshold value for decision rule HR<sub>go</sub>:"), 
                      min = 0.7, max = 0.95, value = c(0.82,0.9), step= 0.01),
       sliderInput("stepHRgo", HTML("Step size for HR<sub>go</sub>:"), 
                   min = 0.01, max = 0.1, value = 0.02, step=0.01),
       
       conditionalPanel("input.Select==1|input.Select==3", 
       sliderInput("R", HTML("Multiplicative adjustment parameter &lambda;:"), 
                   min = 0.2, max = 1, value = c(0.7,0.9), step=0.05),
       sliderInput("stepR", HTML("Step size for &lambda;:"), 
                   min = 0.01, max = 0.1, value = 0.05, step=0.01)),
       

       
       
       conditionalPanel("input.Select==2|input.Select==3",
                        sliderInput("alphaL", HTML("Multiplicative adjustment parameter &alpha;<sub>CI</sub>:"), 
                                    min = 0.025, max = 0.5, value = c(0.3,0.35), step=0.025),
       sliderInput("stepalphaL", HTML("Step size for &alpha;<sub>CI</sub>:"), 
                   min = 0.025, max = 0.1, value = 0.025, step=0.025)),
       
          sliderInput("alpha", HTML("Significance level &alpha;:"), 
                      min = 0.005, max = 0.2, value = 0.05, step= 0.005),
          sliderInput("beta", HTML("1 - Power = &beta;:"), 
                      min = 0.05, max = 0.3, value = 0.1, step= 0.01),
          sliderInput("HR", HTML("Hazard ratio HR=exp(-&beta;):"), 
                      min = 0.4, max = 1.2, value = 0.75, step= 0.01),
          sliderInput("p2", HTML("Event rate for phase II &xi;<sub>2</sub>:"), 
                      min = 0, max = 2, value = 0.45, step= 0.01),
       sliderInput("p3", HTML("Event rate for phase phase III &xi;<sub>3</sub>:"), 
                   min = 0, max = 2, value = 0.45, step= 0.01),
          numericInput("c02", HTML("Fixed costs for phase II (in 10<sup>5</sup>) c<sub>02</sub>:"), 
                      min = 0, max = 200, value = 100, step= 1),
          numericInput("c2", HTML("Per patient costs for phase II (in 10<sup>5</sup>) c<sub>2</sub>:"), 
                      min = 0, max = 2, value = 0.75, step= 0.01),
          numericInput("c03", HTML("Fixed costs for phase III (in 10<sup>5</sup>) c<sub>03</sub>:"), 
                      min = 0, max = 300, value = 150, step= 1),
          numericInput("c3", HTML("Per patient costs for phase III (in 10<sup>5</sup>) c<sub>3</sub>:"), 
                      min = 0, max = 3, value = 1, step= 0.01),
          numericInput("b1", HTML("Amount of benefit for small effect size (in 10<sup>5</sup>) b<sub>1</sub>:"), 
                      min = 1000, max = 6000, value = 1000, step=100),
          numericInput("b2", HTML("Amount of benefit for medium effect size (in 10<sup>5</sup>) b<sub>2</sub>:"), 
                      min = 1000, max = 6000, value = 3000, step=100),
          numericInput("b3", HTML("Amount of benefit for large effect size (in 10<sup>5</sup>) b<sub>3</sub>:"), 
                      min = 1000, max = 6000, value = 5000, step=100), 
       selectInput("refresh", label = h3("Refresh table?"), 
                   choices = list("No" = 0,
                                  "Yes" = 1), 
                   selected = 1),
       selectInput("Plot", label = h3("Plot optimization region?"), 
                   choices = list("No" = 0,
                                  "Yes" = 1), 
                   selected = 1),
       actionButton("go", "Go"),
               
   tags$head(tags$style("#plot{height:200vh !important;}"))),
          
          
          # Show a table summarizing the values entered
          mainPanel(
          includeMarkdown("help51.md"),
          img(src = "biasdesign.png",width = 800),
          includeMarkdown("help52.md"),
          tableOutput("table"),
          includeMarkdown("help53.md"),
          plotlyOutput("plot", height = "30px")
          )
     )
))
