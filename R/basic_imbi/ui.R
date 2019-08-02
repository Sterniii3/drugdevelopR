library(shiny)
library(markdown)
library(plotly)

shinyUI(
   fluidPage(
      titlePanel("drugdevelopR: basic"),
      sidebarLayout(
         sidebarPanel(

   selectInput("Select", 
               label = h3("Select type of endpoint"),
               choices = list("Time-to-event " = 1,
                              "Binary " = 2,
                              "Normally distributed " = 3),
          selected = 1),

   conditionalPanel("input.Select==1",
   sliderInput("HR", HTML("Hazard ratio HR=exp(-&beta;):"),
      min = 0.4, max = 1.2, value = 0.7, step= 0.01),
   sliderInput("HRgo", 
      HTML("Threshold value for decision rule 
           HR<sub>go</sub>:"),
      min = 0.65, max = 0.95, value = c(0.82,0.94), step= 0.01),
   sliderInput("stepHRgo", HTML("Step size for 
                                HR<sub>go</sub>:"),
      min = 0.01, max = 0.1, value = 0.02, step=0.01),
   sliderInput("Num1", 
      HTML("Total number of events for phase II 
           d<sub>2</sub>:"),
      min = 10, max = 500, value = c(250,350), step=5),
   sliderInput("stepNum1", HTML("Step size for 
                                d<sub>2</sub>:"),
      min = 1, max = 50, value = 5, step=2),
   sliderInput("p2", 
      HTML("Event rate for phase II &xi;<sub>2</sub>"),
      min = 0, max = 2, value = 0.7, step= 0.01),
   sliderInput("p3", 
      HTML("Event rate for phase III &xi;<sub>3</sub>"),
      min = 0, max = 2, value = 0.7, step= 0.01),
   sliderInput("alpha1", HTML("Significance level &alpha;:"),
      min = 0.005, max = 0.2, value = 0.05, step= 0.005),
   sliderInput("beta1", HTML("1 - Power = &beta;:"),
      min = 0.05, max = 0.3, value = 0.1, step= 0.01),
   numericInput("c021", 
      HTML("Fixed costs for phase II 
            (in 10<sup>5</sup>) c<sub>02</sub>:"),
      min = 0, max = 200, value = 100, step= 1),
   numericInput("c21", 
      HTML("Per patient costs for phase II 
            (in 10<sup>5</sup>) c<sub>2</sub>:"),
      min = 0, max = 2, value = 0.75, step= 0.01),
   numericInput("c031", 
      HTML("Fixed costs for phase III 
            (in 10<sup>5</sup>) c<sub>03</sub>:"),
      min = 0, max = 300, value = 150, step= 1),
   numericInput("c31", 
      HTML("Per patient costs for phase III 
            (in 10<sup>5</sup>) c<sub>3</sub>:"),
      min = 0, max = 3, value = 1, step= 0.01),
   sliderInput("small1", 
      HTML("Lower boundary for small effect 
            size category in HR scale:"),
      min = 0, max = 1.2, value = 1, step= 0.01),
   sliderInput("medium1", 
      HTML("Lower boundary for medium effect 
            size category in HR scale:"),
      min = 0, max = 1.2, value = 0.95, step= 0.01),
   sliderInput("large1", 
      HTML("Lower boundary for large effect 
            size category in HR scale:"),
      min = 0, max = 1.2, value = 0.85, step=0.01),
   numericInput("b11", 
      HTML("Amount of benefit for small effect size 
           (in 10<sup>5</sup>) b<sub>1</sub>:"),
      min = 1000, max = 6000, value = 1000, step=100),
   numericInput("b21", 
      HTML("Amount of benefit for medium effect size 
          (in 10<sup>5</sup>) b<sub>2</sub>:"),
      min = 1000, max = 6000, value = 3000, step=100),
   numericInput("b31", 
      HTML("Amount of benefit for large effect size 
           (in 10<sup>5</sup>) b<sub>3</sub>:"),
      min = 1000, max = 6000, value = 5000, step=100),
   numericInput("K1", 
      HTML("Maximal costs for phase II and III 
           (in 10<sup>5</sup>) K:"),
      min = 100, max = 100000, value = 5000, step=100),
   selectInput("refresh1", label = h3("Refresh table?"),
      choices = list("No" = 0,
                "Yes" = 1),selected = 1),
   selectInput("Plot1", 
               label = h3("Plot optimization region?"),
      choices = list("No" = 0,
                "Yes" = 1),
      selected = 1),
   actionButton("go1", "Go")
   ),

   conditionalPanel("input.Select==2",
   sliderInput("p0", 
      HTML("Assumed control rate p<sub>0</sub>:"),
      min = 0.01, max = 1, value = 0.5, step= 0.01),
   sliderInput("p1", 
      HTML("Assumed treatment rate p<sub>1</sub>:"),
      min = 0.01, max = 1, value = 0.3, step= 0.01),
   sliderInput("Num2", 
      HTML("Total sample size in phase II 
           n<sub>2</sub>:"),
      min = 10, max = 500, value = c(250,350), step=5),
   sliderInput("stepNum2", 
      HTML("Step size for n<sub>2</sub>:"),
      min = 2, max = 50, value = 10, step=2),
   sliderInput("RRgo", 
      HTML("Threshold value for decision rule 
           RR<sub>go</sub>:"),
      min = 0.65, max = 0.95, value = c(0.82,0.94), step= 0.01),
   sliderInput("stepRRgo", 
      HTML("Step size for RR<sub>go</sub>:"),
      min = 0.01, max = 0.1, value = 0.02, step=0.01),
   sliderInput("alpha2", HTML("Significance level &alpha;:"),
      min = 0.005, max = 0.2, value = 0.05, step= 0.005),
   sliderInput("beta2", HTML("1 - Power = &beta;:"),
      min = 0.05, max = 0.3, value = 0.1, step= 0.01),
   numericInput("c022", 
      HTML("Fixed costs for phase II 
           (in 10<sup>5</sup>) c<sub>02</sub>:"),
      min = 0, max = 200, value = 100, step= 1),
   numericInput("c22", 
      HTML("Per patient costs for phase II 
           (in 10<sup>5</sup>) c<sub>2</sub>:"),
      min = 0, max = 2, value = 0.75, step= 0.01),
   numericInput("c032", 
      HTML("Fixed costs for phase III 
           (in 10<sup>5</sup>) c<sub>03</sub>:"),
      min = 0, max = 300, value = 150, step= 1),
   numericInput("c32", 
      HTML("Per patient costs for phase III 
           (in 10<sup>5</sup>) c<sub>3</sub>:"),
      min = 0, max = 3, value = 1, step= 0.01),
   sliderInput("small2", 
      HTML("Lower boundary for small effect 
           size category in RR scale:"),
      min = 0, max = 1.2, value = 1, step= 0.01),
   sliderInput("medium2", 
      HTML("Lower boundary for medium effect 
            size category in RR scale:"),
      min = 0, max = 1.2, value = 0.95, step= 0.01),
   sliderInput("large2", 
      HTML("Lower boundary for large effect 
            size category in RR scale:"),
      min = 0, max = 1.2, value = 0.85, step=0.01),
   numericInput("b12", 
      HTML("Amount of benefit for small effect 
            size (in 10<sup>5</sup>) b<sub>1</sub>:"),
      min = 1000, max = 6000, value = 1000, step=100),
   numericInput("b22", 
      HTML("Amount of benefit for medium effect size 
           (in 10<sup>5</sup>) b<sub>2</sub>:"),
      min = 1000, max = 6000, value = 3000, step=100),
   numericInput("b32", 
      HTML("Amount of benefit for large effect size 
           (in 10<sup>5</sup>) b<sub>3</sub>:"),
      min = 1000, max = 6000, value = 5000, step=100),
   numericInput("K2", 
      HTML("Maximal costs for phase II and III 
           (in 10<sup>5</sup>) K:"),
      min = 100, max = 100000, value = 5000, step=100),
   selectInput("refresh2", label = h3("Refresh table?"),
      choices = list("No" = 0, "Yes" = 1), 
      selected = 1),
   selectInput("Plot2", 
               label = h3("Plot optimization region?"),
      choices = list("No" = 0,"Yes" = 1),
      selected = 0),
   actionButton("go2", "Go")),

   conditionalPanel("input.Select==3",
   sliderInput("Delta", 
      HTML("Assumed standardized 
            difference in means &Delta;:"),
      min = 0.01, max = 1, value = 0.7, step= 0.01),
   sliderInput("kappa", 
      HTML("Threshold value for decision rule &kappa;:"),
      min = 0.01, max = 1, value = c(0.05,0.3), step= 0.01),
   sliderInput("stepkappa", HTML("Step size for &kappa;:"),
      min = 0.01, max = 0.1, value = 0.02, step=0.01),
   sliderInput("Num3", 
      HTML("Total sample size in phase II n<sub>2</sub>:"),
      min = 10, max = 500, value = c(20,150), step=5),
   sliderInput("stepNum3", 
      HTML("Step size for n<sub>2</sub>:"),
      min = 2, max = 50, value = 10, step=2),
   sliderInput("alpha3", HTML("Significance level &alpha;:"),
      min = 0.005, max = 0.2, value = 0.05, step= 0.005),
   sliderInput("beta3", HTML("1 - Power = &beta;:"),
      min = 0.05, max = 0.3, value = 0.1, step= 0.01),
   numericInput("c023", 
      HTML("Fixed costs for phase II 
           (in 10<sup>5</sup>) c<sub>02</sub>:"),
      min = 0, max = 200, value = 100, step= 1),
   numericInput("c23", 
      HTML("Per patient costs for phase II 
           (in 10<sup>5</sup>) c<sub>2</sub>:"),
      min = 0, max = 2, value = 0.75, step= 0.01),
   numericInput("c033", 
      HTML("Fixed costs for phase III 
           (in 10<sup>5</sup>) c<sub>03</sub>:"),
      min = 0, max = 300, value = 150, step= 1),
   numericInput("c33", 
      HTML("Per patient costs for phase III 
           (in 10<sup>5</sup>) c<sub>3</sub>:"),
      min = 0, max = 3, value = 1, step= 0.01),
   sliderInput("small3", 
      HTML("Lower boundary for small effect 
           size category:"),
      min = 0, max = 1, value = 0, step= 0.005),
   sliderInput("medium3", 
      HTML("Lower boundary for medium effect 
           size category:"),
      min = 0, max = 1, value = 0.5, step= 0.005),
   sliderInput("large3", 
      HTML("Lower boundary for large effect 
           size category:"),
      min = 0, max = 1, value = 0.8, step=0.005),
   numericInput("b13", 
      HTML("Amount of benefit for small effect size 
           (in 10<sup>5</sup>) b<sub>1</sub>:"),
      min = 1000, max = 6000, value = 1000, step=100),
   numericInput("b23", 
      HTML("Amount of benefit for medium effect size 
           (in 10<sup>5</sup>) b<sub>2</sub>:"),
      min = 1000, max = 6000, value = 3000, step=100),
   numericInput("b33", 
      HTML("Amount of benefit for large effect size 
           (in 10<sup>5</sup>) b<sub>3</sub>:"),
      min = 1000, max = 6000, value = 5000, step=100),
   numericInput("K3", 
      HTML("Maximal costs for phase II and III 
           (in 10<sup>5</sup>) K:"),
      min = 100, max = 100000, value = 5000, step=100),
   selectInput("refresh3", label = h3("Refresh table?"),
      choices = list("No" = 0, "Yes" = 1),
            selected = 1),
   selectInput("Plot3", 
               label = h3("Plot optimization region?"),
            choices = list("No" = 0, "Yes" = 1),
            selected = 0),
   actionButton("go3", "Go")),
   
   tags$head(tags$style("#plot{height:200vh !important;}"))),

   mainPanel(
      conditionalPanel(
         condition="input.Select == 1",
         includeMarkdown("help11.md"),
         img(src = "trialdesign1.png",width = 900),
         includeMarkdown("help12.md"),
         img(src = "CI1.png",width = 450),
         includeMarkdown("help122.md")),

      conditionalPanel(
         condition="input.Select == 2",
         includeMarkdown("help21.md"),
         img(src = "trialdesign2.png",width = 900),
         includeMarkdown("help22.md"),
         img(src = "CI2.png",width = 450),
         includeMarkdown("help222.md")),

      conditionalPanel(
         condition="input.Select == 3",
         includeMarkdown("help31.md"),
         img(src = "trialdesign3.png",width = 900),
         includeMarkdown("help32.md"),
         img(src = "CI3.png",width = 450),
         includeMarkdown("help322.md")),

         tableOutput("table"),

      conditionalPanel(
         condition="input.Select == 1",
         includeMarkdown("help13.md")),

      conditionalPanel(
         condition="input.Select == 2",
         includeMarkdown("help23.md")),

      conditionalPanel(
         condition="input.Select == 3",
         includeMarkdown("help33.md")),

         plotlyOutput("plot", height = "20px")
      )
  ))
)
