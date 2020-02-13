library("shiny")
shinyServer(function(input, output,session) {
   
   urlbasic <- a("basic R Shiny App", href="https://web.imbi.uni-heidelberg.de/basic/")
   prior <- a("prior", href="https://web.imbi.uni-heidelberg.de/prior/")
   
   
   urlmultiarm <- a("multiarm R Shiny App", href="https://web.imbi.uni-heidelberg.de/multiarm/")

   
   urlmultitrial <- a("multitrial R Shiny App", href="https://web.imbi.uni-heidelberg.de/multitrial/")

   
   urlbias <- a("bias R Shiny App", href="https://web.imbi.uni-heidelberg.de/bias/")

   urlpackage <- a("drugdevelopR package", href="https://github.com/Sterniii3/drugdevelopR")
   
   output$tabbasic <- renderUI({
      tagList("- with fixed assumed treatment effects:", urlbasic)
   })

 
   
   output$tabbasic1 <- renderUI({
     tagList("- with", prior, "distribution for the assumed treatment effects:", urlpackage)
   })
   

   output$tabmultiarm <- renderUI({
      tagList("- with fixed assumed treatment effects:", urlmultiarm, ",", urlpackage)
   })
   
 
   output$tabmultitrial <- renderUI({
      tagList("- with fixed assumed treatment effects:", urlmultitrial)
   })
   output$tabmultitrial1 <- renderUI({
     tagList("- with", prior, "distribution for the assumed treatment effects:", urlpackage) })  

   
   
   output$tabbias <- renderUI({
      tagList("- with fixed assumed treatment effects:", urlbias)
   })
   
   output$tabbias1 <- renderUI({
     tagList("- with", prior, "distribution for the assumed treatment effects:", urlpackage) })
   
})