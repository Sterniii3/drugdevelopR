library(shiny)
library(markdown)
library(plotly)

# Define UI for slider demo application
shinyUI(fluidPage(
     
     #  Application title
   titlePanel("drugdevelopR: multitrial"),
     # Sidebar with sliders that demonstrate various available
     # options  
  sidebarLayout(
    sidebarPanel(
       selectInput("Case", label = h3("Select case"), 
                   choices = list("Case 1: At least one significant trial" = 1,
                                  "Case 2: At least two significant trials" = 2,
                                  "Case 3: At least three significant trials" = 3), 
                   selected = 1),

          sliderInput("D2", HTML("Total number of events for phase II d<sub>2</sub>:"), 
                             min = 10, max = 750, value = c(250,350), step=2),
       sliderInput("stepD2", HTML("Step size for d<sub>2</sub>:"), 
                   min = 2, max = 50, value = 10, step=2),
          sliderInput("HRgo", HTML("Threshold value for decision rule HR<sub>go</sub>:"), 
                      min = 0.7, max = 0.95, value = c(0.82,0.9), step= 0.01),
       sliderInput("stepHRgo", HTML("Step size for HR<sub>go</sub>:"), 
                   min = 0.01, max = 0.1, value = 0.02, step=0.01),
          sliderInput("alpha", HTML("Significance level &alpha;:"), 
                      min = 0.005, max = 0.2, value = 0.025, step= 0.005),
          sliderInput("beta", HTML("1 - Power = &beta;:"), 
                      min = 0.05, max = 0.3, value = 0.1, step= 0.01),
          sliderInput("HR", HTML("Hazard ratio HR=exp(-&beta;):"), 
                      min = 0.4, max = 1.2, value = 0.75, step= 0.01),
          sliderInput("p2", HTML("Event rate for phase II &xi;<sub>2</sub>:"), 
                      min = 0, max = 2, value = 0.7, step= 0.01),
       sliderInput("p3", HTML("Event rate for phase II &xi;<sub>3</sub>:"), 
                   min = 0, max = 2, value = 0.7, step= 0.01),
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
          img(src = "cases_strategies1.png",width = 900),
          includeMarkdown("help62.md"),
          img(src = "vs.png",width = 1150),
          includeMarkdown("help63.md"),
          tableOutput("table"),
          includeMarkdown("help64.md"),
          plotlyOutput("plot", height = "30px")
          )
     )
))
