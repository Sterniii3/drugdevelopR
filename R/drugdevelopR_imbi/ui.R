shinyUI(fluidPage(
   
   
   
   
   
   
   titlePanel("drugdevelopR: Optimal Phase II/III Drug Development Planning"),
   includeMarkdown("shortexpl.md"),
   img(src = "basic1.png",width = 490),
   uiOutput("tabbasic"),
   uiOutput("tabbasic1"),
   includeMarkdown("blank.md"),
   #includeMarkdown("shortexpl2.md"),
   img(src = "bias1.png",width = 490),
   uiOutput("tabbias"),
   uiOutput("tabbias1"),
   img(src = "multitrial1.png",width = 490),
   uiOutput("tabmultitrial"),
   uiOutput("tabmultitrial1"),
   img(src = "multiarm1.png",width = 490),
   uiOutput("tabmultiarm"),
   uiOutput("tabmultiarm1"),
   includeMarkdown("maintainer.md")
   
   
   
))