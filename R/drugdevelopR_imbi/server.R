library("shiny")
shinyServer(function(input, output,session) {
   
   urlbasic <- a("basic R Shiny App", href="https://web.imbi.uni-heidelberg.de/basic/")
   urlbasic1 <- a("basic R package", href="https://web.imbi.uni-heidelberg.de/basic_package/")
   prior <- a("prior", href="https://web.imbi.uni-heidelberg.de/prior/")
   
   
   urlmultiarm <- a("multiarm R Shiny App", href="https://web.imbi.uni-heidelberg.de/multiarm/")
   urlmultiarm1 <- a("multiarm R package", href="https://web.imbi.uni-heidelberg.de/multiarm_package/")

   
   urlmultitrial <- a("multitrial R Shiny App", href="https://web.imbi.uni-heidelberg.de/multitrial/")
   urlmultitrial1 <- a("multitrial R package", href="https://web.imbi.uni-heidelberg.de/multitrial_package/")

   
   urlbias <- a("bias R Shiny App", href="https://web.imbi.uni-heidelberg.de/bias/")
   urlbias1 <- a("bias R package", href="https://web.imbi.uni-heidelberg.de/multiarm_package/")

   
   output$tabbasic <- renderUI({
      tagList("- with fixed assumed treatment effects:", urlbasic)
   })

   #output$tabbasic1 <- renderUI({
   #   tagList("- with", prior, "distribution for the assumed treatment effects:", urlbasic1)
   #})
   
   output$tabbasic1 <- renderUI({
     tagList("- with", prior, "distribution for the assumed treatment effects")
   })
   

   output$tabmultiarm <- renderUI({
      tagList("- with fixed assumed treatment effects:", urlmultiarm)
   })
   
  # output$tabmultiarm1 <- renderUI({
   #   tagList("- with", prior, "distribution for the assumed treatment effects:", urlmultiarm1)})
   
   output$tabmultitrial <- renderUI({
      tagList("- with fixed assumed treatment effects:", urlmultitrial)
   })
   
   #output$tabmultitrial1 <- renderUI({
  #    tagList("- with", prior, "distribution for the assumed treatment effects:", urlmultitrial1)
   #})
   
   output$tabmultitrial1 <- renderUI({
     tagList("- with", prior, "distribution for the assumed treatment effects")
   })
   
   
   output$tabbias <- renderUI({
      tagList("- with fixed assumed treatment effects:", urlbias)
   })
   
   #output$tabbias1 <- renderUI({
   #   tagList("- with", prior, "distribution for the assumed treatment effects:", urlbias1) })
   
})