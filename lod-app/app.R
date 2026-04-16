library(shiny)
library(tidyverse)

ui <- fluidPage(
  
  titlePanel("Statistical limit of detection for next generation sequencing"),
  
  "Check out ",
    tags$a(href= "https://joe-m-shaw.github.io/professional_website/posts/coverage_and_confidence.html",
           " this blogpost"),
  " for the full explantion. The code is on ",
  tags$a(href = "https://github.com/joe-m-shaw/limit-of-detection-app",
        "Github."),

  sidebarLayout(
    sidebarPanel(
      numericInput("true_vaf",
                   "True variant allele frequency (%)",
                   min = 0,
                   max = 100,
                   value = 10),
      numericInput("vaf_threshold",
                   "Pipeline variant allele frequency threshold (%)",
                   value = 4, 
                   min = 0, 
                   max = 100),
      numericInput("coverage", 
                   "Coverage (X)", 
                   value = 138, 
                   min = 0, 
                   max = 5000),
      numericInput("x_max", 
                   "X axis max", 
                   value = 50, 
                   min = 0, 
                   max = 5000),
    ),
    mainPanel(
      plotOutput("lod_plot")
    )
  )
)

server <- function(input, output) {
  
  output$lod_plot <- renderPlot({
    
    number_of_reads_needed_to_call <- round(input$coverage*(input$vaf_threshold/100), 0)
    
    number_of_reads_which_wont_call <- number_of_reads_needed_to_call-1
    
    probability_of_not_calling_variant <- pbinom(q = number_of_reads_which_wont_call,
                                                 size = input$coverage,
                                                 prob = input$true_vaf/100)
    
    probability_of_calling_variant <- 1-probability_of_not_calling_variant
    
    df_data_for_plot <- data.frame(
      n_reads = seq(0, input$x_max)) |> 
      mutate(y_prob = dbinom(x = n_reads, 
                             size = input$coverage, 
                             prob = input$true_vaf/100),
             status = case_when(
               n_reads >= number_of_reads_needed_to_call ~"Detected",
               n_reads < number_of_reads_needed_to_call ~"Not detected"))
    
    ggplot(df_data_for_plot, aes(x = n_reads, y = y_prob)) +
      geom_col(aes(fill = status), 
               colour = "black") +
      scale_fill_manual(values = c("blue", "red")) +
      geom_vline(xintercept = number_of_reads_needed_to_call, 
                 linetype = "dashed") +
      annotate(geom = "text", 
               x = number_of_reads_needed_to_call - 1,
               y = 0.05, 
               label = "VAF threshold", 
               angle = 90) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_x_continuous(breaks = c(seq(0, 
                                        input$x_max, 
                                        by = 10),
                                    number_of_reads_needed_to_call)) +
      labs(x = "Number of reads",
           y = "Probability of detecting variant in number of reads",
           subtitle = paste0("If coverage is ", 
                             input$coverage, 
                             "X and true VAF is ", 
                             input$true_vaf,
                             "% and VAF detection threshold is ",
                             input$vaf_threshold, 
                             "% (",
                             number_of_reads_needed_to_call,
                             " reads)"),
           title = paste0("Probability of detecting variant: ",
                          round(probability_of_calling_variant*100, 2), 
                          "%") ,
           fill = "Variant status") 
    
  })
}

shinyApp(ui = ui, server = server)
