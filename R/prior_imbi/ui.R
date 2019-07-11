library(shiny)
library(markdown)
library(plotly)
library(viridis)
library(msm)

# Define UI for slider demo application
shinyUI(fluidPage(

     #  Application title
   titlePanel("drugdevelopR: prior"),
     # Sidebar with sliders that demonstrate various available
     # options
  sidebarLayout(
    sidebarPanel(
       selectInput("Select", label = h3("Select type of endpoint"),
                   choices = list("Time-to-event " = 1,
                                  "Binary " = 2,
                                  "Normally distributed " = 3),
                   selected = 1),

conditionalPanel("input.Select==1",
          sliderInput("HR", HTML("Assumed true treatment effects on HR scale (HR<sub>1</sub>, HR<sub>2</sub>):"),
                             min = 0.4, max = 1.2, value = c(0.69, 0.88), step=0.01),
       sliderInput("n1", HTML("Amount of information for HR<sub>1</sub> and HR<sub>2</sub> in terms of number of events (i<sup>d</sup><sub>1</sub>, i<sup>d</sup><sub>2</sub>):"),
                   min = 10, max = 1000, value = c(210, 420), step=10),
       sliderInput("w1", HTML("Weights w:"),
                   min = 0, max = 1, value = c(0.3,0.9), step=0.05)),



conditionalPanel("input.Select==2",
                 sliderInput("p0", HTML("Assumed true rate of control group p<sub>0</sub>:"),
                             min = 0.1, max = 1, value = 0.6, step=0.01),
                 sliderInput("p1", HTML("Assumed true rates of treatment group (p<sub>1,1</sub>, p<sub>1,2</sub>):"),
                             min = 0, max = 1, value = c(0.3, 0.5), step=0.02),
                 sliderInput("n2", HTML("Amount of information for p<sub>1,2</sub> and p<sub>1,2</sub> in terms of sample size (i<sup>n</sup><sub>1</sub>, i<sup>n</sup><sub>2</sub>):"),
                             min = 10, max = 1000, value = c(30, 60), step=10),
                 sliderInput("w2", HTML("Weights w:"),
                             min = 0, max = 1, value = c(0.1,0.7), step=0.05)),


conditionalPanel("input.Select==3",
                 sliderInput("delta", HTML("Assumed true treatment effects for standardized difference in means (&Delta;<sub>1</sub>, &Delta;<sub>2</sub>):"),
                             min = 0.01, max = 1, value = c(0.25, 0.625), step=0.005),
                 sliderInput("n3", HTML("Amount of information for &Delta;<sub>2</sub> and &Delta;<sub>1</sub> in terms of sample size (i<sup>n</sup><sub>1</sub>, i<sup>n</sup><sub>2</sub>):"),
                             min = 10, max = 1000, value = c(300, 600), step=10),
                 sliderInput("b", HTML("Lower and upper boundary for truncation (a, b):"),
                             min = 0, max = 1, value = c(0, 0.75), step=0.01),
                 sliderInput("w3", HTML("Weights w:"),
                             min = 0, max = 1, value = c(0.2,0.8), step=0.05)),

conditionalPanel(condition="input.Select==1",img(src = "prior1.png",width = 425),
                 sliderInput("gamma1", HTML("Assume more pessimistic (&gamma;<0)/ optimistic (&gamma;>0) view about treatment effect in phase III (&theta;<sub>3</sub>) compared to phase II (&theta;<sub>2</sub>):"),
                             min = -0.05, max = 0.05, value = 0, step=0.0125),
                 img(src = "pop11.png",width = 425)
                 ),
conditionalPanel(condition="input.Select==2",img(src = "prior2.png",width = 425),
                 sliderInput("gamma2", HTML("Assume more pessimistic (&gamma;<0)/ optimistic (&gamma;>0) view about treatment effect in phase III (p<sub>1,3</sub>) compared to phase II (p<sub>1,2</sub>):"),
                             min = -0.05, max = 0.05, value = 0, step=0.01),
                 img(src = "pop12.png",width = 425)),
conditionalPanel(condition="input.Select==3",img(src = "prior3.png",width = 425),
                 sliderInput("gamma3", HTML("Assume more pessimistic (&gamma;<0)/ optimistic (&gamma;>0) view about treatment effect in phase III (&Delta;<sub>3</sub>) compared to phase II (&Delta;<sub>2</sub>):"),
                             min = -0.05, max = 0.05, value = 0, step=0.005),
                 img(src = "pop13.png",width = 425))

),


          # Show a table summarizing the values entered
          mainPanel(
          #includeMarkdown("help51.md"),
          #img(src = "biasdesign.png",width = 800),
          #includeMarkdown("help52.md"),
          #includeMarkdown("help53.md"),
         # plotOutput("plot")
             fluidRow(
                column(8, plotOutput("plot1")),
                column(12, plotOutput("plot"))
             ),
             conditionalPanel(condition="input.Select==1",includeMarkdown("help1.md")),
             conditionalPanel(condition="input.Select==2",includeMarkdown("help2.md")),
             conditionalPanel(condition="input.Select==3",includeMarkdown("help3.md"))
          )
     )
))
