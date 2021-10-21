library(dplyr)
library(tibble)
library(shiny)
library(shinydashboard)
library(ggplot2)
library(plotly)

ui <- shinyUI(
  dashboardPage(
    dashboardHeader(title = "Lipid QC"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Files", selected = TRUE,
                  fileInput('file1',
                           'lipids LC-MS data',
                           accept=c('text/csv'))

        ),
        actionButton("analyze", "QC Analyze")
    )
    ),
    dashboardBody(

      fluidRow(
        h1(strong("QC plots")),
        tabBox(width = 12,
                     tabPanel(title = "plot1",
                               plotlyOutput("plot1"),
                               downloadButton('downloadPlot1', 'Save plot')),
                      tabPanel(title = "plot2",
                               plotlyOutput("plot2"),
                              downloadButton('downloadPlot2', 'Save plot')),
                      tabPanel(title = "plot3",
                              plotlyOutput("plot3"),
                              downloadButton('downloadPlot3', 'Save plot'))
                  )

    ),

      fluidRow(
        h1("Data"),
        DT::dataTableOutput("table")

      )
      )
    )
  )


server <- shinyServer(function(input, output) {
  options(shiny.maxRequestSize=200*1024^2)


  ### Reactive functions ### --------------------------------------------------
  data <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header = TRUE,
             sep = ",", stringsAsFactors = FALSE) %>%
      mutate(id = row_number())
  })




  ### All object and functions upon 'Analyze' input  ### ----------------------

  observeEvent(input$analyze, {



    ### Reactive functions ### ------------------------------------------------


   plot1_input <- reactive({
      withProgress(message = 'Plotting', value = 0.66, {
        p <- ggplot(data(), aes(Petal.Length, Petal.Width))+geom_point()
        ggplotly(p)
      })
    })

    plot2_input <- reactive({
      p <- ggplot(data(), aes(Species, Sepal.Width))+geom_boxplot()
      ggplotly(p)
    })

    plot3_input <- reactive({
        p <- ggplot(data(), aes(Species, Sepal.Width))+geom_bar(stat='identity')
        ggplotly(p)
    })




    ### Output functions ### --------------------------------------------------
    output$table <- DT::renderDataTable({
      data()
    }, options = list(pageLength = 5, scrollX = T),
    selection = list(selected = c(1)))



    output$plot1 <- renderPlotly({
      plot1_input()
    })

    output$plot2 <- renderPlotly({
      plot2_input()
    })

    output$plot3 <- renderPlotly({
      plot3_input()
    })


    ### Download objects and functions ### ------------------------------------
    datasetInput <- reactive({
      switch(input$dataset,
             "results" = get_results(dep()),
             "significant_proteins" = get_results(dep()) %>%
               filter(significant) %>%
               select(-significant),
             "displayed_subset" = res() %>%
               filter(significant) %>%
               select(-significant),
             "full_dataset" = get_df_wide(dep()))
    })

    output$downloadData <- downloadHandler(
      filename = function() { paste(input$dataset, ".txt", sep = "") },
      content = function(file) {
        write.table(datasetInput(),
                    file,
                    col.names = TRUE,
                    row.names = FALSE,
                    sep ="\t") }
    )

    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0("Barplot_", table()[input$table_rows_selected,1], ".pdf")
      },
      content = function(file) {
        pdf(file)
        print(selected_plot_input())
        dev.off()
      }
    )

    output$downloadPlot1 <- downloadHandler(
      filename = 'plot1.pdf',
      content = function(file) {
        pdf(file, height = 10, paper = "a4")
        print(plot1_input())
        dev.off()
      }
    )
    output$downloadPlot2 <- downloadHandler(
      filename = 'plot2.pdf',
      content = function(file) {
        pdf(file, height = 10, paper = "a4")
        print(plot2_input())
        dev.off()
      }
    )

    output$downloadPlot3 <- downloadHandler(
      filename = 'plot3.pdf',
      content = function(file) {
        pdf(file, height = 10, paper = "a4")
        print(plot3_input())
        dev.off()
      }
    )



  })
})


shinyApp(ui = ui, server = server)
