library(shiny)
library(markdown)
library(plotly)
library(viridis)
library(shinyWidgets)
# Define UI for slider demo application
shinyUI(fluidPage(
  chooseSliderSkin(skin = "Shiny", color = viridis(3,alpha=0.9)[2]),
     #  Application title
   titlePanel("drugdevelopR: multiarm"),
     # Sidebar with sliders that demonstrate various available
     # options  
  sidebarLayout(
    sidebarPanel(
       selectInput("Select", label = h3("Select strategy"), 
                   choices = list("Strategy 1" = 1,
                                  "Strategy 2" = 2,
                                  "Strategy 1 vs. 2" = 3), 
                   selected = 1),

          sliderInput("N2", HTML("Total sample size for (three-arm) phase II n<sub>2</sub>:"), 
                             min = 15, max = 750, value = c(60,210), step=3),
       sliderInput("stepN2", HTML("Step size for n<sub>2</sub>:"), 
                   min = 3, max = 60, value = 30, step=3),
          sliderInput("HRgo", HTML("Threshold value for decision rule HR<sub>go</sub>:"), 
                      min = 0.7, max = 0.95, value = c(0.82,0.88), step= 0.01),
       sliderInput("stepHRgo", HTML("Step size for HR<sub>go</sub>:"), 
                   min = 0.01, max = 0.1, value = 0.02, step=0.01),
          sliderInput("alpha", HTML("One-sided significance level &alpha;:"), 
                      min = 0.005, max = 0.2, value = 0.025, step= 0.005),
          sliderInput("beta", HTML("1 - Power = &beta;:"), 
                      min = 0.05, max = 0.3, value = 0.1, step= 0.01),
          sliderInput("HR1", HTML("Assumed true hazard ratio HR<sub>1</sub>=exp(-&beta;<sub>1</sub>):"), 
                      min = 0.4, max = 1.2, value = 0.8, step= 0.01),
          sliderInput("HR2", HTML("Assumed true hazard ratio HR<sub>2</sub>=exp(-&beta;<sub>2</sub>):"), 
                      min = 0.4, max = 1.2, value = 0.8, step= 0.01),
          sliderInput("ec", HTML("Control arm event rate for phase II and III e<sub>c</sub>:"), 
                      min = 0, max = 2, value = 0.6, step= 0.01),
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
       numericInput("Nc", HTML("Maximal expected sample size for phase II and III N:"), 
                    min = 100, max = 100000, value = 5000, step=100), 
       selectInput("refresh", label = h3("Refresh table?"), 
                   choices = list("No" = 0,
                                  "Yes" = 1), 
                   selected = 1),
       selectInput("Plot", label = h3("Plot optimization region?"), 
                   choices = list("No" = 0,
                                  "Yes" = 1), 
                   selected = 1),
       actionButton("go", "Go"),
               
   tags$head(tags$style("#plot{height:75vh !important;}"))),
          
          
          # Show a table summarizing the values entered
          mainPanel(
          includeMarkdown("help51.md"),
          img(src = "Trial Design_small5.png",width = 1200),
          includeMarkdown("help52.md"),
          tableOutput("table"),
          includeMarkdown("help53.md"),
          includeMarkdown("help54.md"),
          plotlyOutput("plot", height = "30px")
          )
     )
))
