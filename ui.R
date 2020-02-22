#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)


# Define UI for application that draws a histogram
shinyUI(
    fluidPage(
        theme = shinytheme("flatly"),
        tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "mystyle.css")),   
        navbarPage(
            
            title = strong("Hydroxychloroquine suflate (HCQ suflate)"),
            
            tabPanel(
                strong("HCQ concentration prediction"),
                fluidRow(
                    column(4, 
                           # h3("Dose regimen"),
                           fluidRow(
                               column(6),
                               column(6,
                                      sliderInput("cycle_pk", label = "cycles", min = 1, max = 28, step = 1, value = c(1))
                               )
                           ),
                           fluidRow(
                               column(6,
                                      sliderInput("amt_pk1", label = "HCQ sulfate dose (mg)", min = 0, max = 1000, step = 50, value = c(100))
                               ),
                               column(6,
                                      sliderInput("amt_pk2", label = "dose interruption (mg)", min = 0, max = 1000, step = 50, value = c(0))
                               )
                           ),
                           fluidRow(
                               column(6,
                                      sliderInput("interval_pk1", label = "dose interval (h)", min = 0, max = 96, step = 2, value = c(12))
                               ),
                               column(6,
                                      sliderInput("interval_pk2", label = "interrupted interval (h)", min = 0, max = 96, step = 2, value = c(0))
                               )
                           ),
                           fluidRow(
                               column(6,
                                      sliderInput("n_pk1", label = "dose times", min = 1, max = 28, step = 1, value = c(14))
                               ),
                               column(6,
                                      sliderInput("n_pk2", label = "interruption times", min = 1, max = 28, step = 1, value = c(1))
                               )
                           )
                    ),
                column(8,
                       plotOutput("distPlotPK"))
                ),
                
                fluidRow(
                    column(2,
                           sliderInput("cl_pk", label = "CL/F (L/h)", min = 1, max = 100, step = 0.5, value = c(10.9))),
                    column(2,
                           sliderInput("vc_pk", label = "Vc (L)", min = 10, max = 2000, step = 10, value = c(437))),
                    column(2,
                           sliderInput("vp_pk", label = "Vp (L)", min = 10, max = 2000, step = 10, value = c(1390))),
                    column(2,
                           sliderInput("q_pk", label = "Q/F (L)", min = 10, max = 100, step = 10, value = c(45.1))),
                    column(2,
                           sliderInput("ka_pk", label = "Ka (1/h)", min = 0.1, max = 2, step = 0.05, value = c(1.15))),
                    column(2,
                           sliderInput("p_pk", label = "tissue penetration", min = 1, max = 1000, step = 10, value = c(1)))
                ),
                # absolutePanel(
                #     top = 130, right = 30, width = 200, height = 10, draggable = TRUE,
                #     HTML(
                #         paste0("<strong>C<sub>", "max", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                #                textOutput(outputId = "pkcmax1", inline = T), "</code> ng/mL</strong>")
                #     )
                # ),
                # absolutePanel(
                #     top = 190, right = 20, width = 250, height = 10, draggable = TRUE,
                #     HTML(
                #         paste0("<strong>AUC<sub>", "tau", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                #                textOutput(outputId = "pkauc1", inline = T), "</code> hr*ng/mL</strong>")
                #     )
                # ),
                absolutePanel(
                    top = 110, left = 55, width = 100, height = 10, draggable = TRUE,
                    img(src="LOGOdMed.png", height = 50)
                )
            ),
            footer = h5(HTML("dMed Copyright 2019 : 
                       <strong style='color:#ec4c3c;background-color:#F8F9F9'> E </strong>
                       arly <strong style='color:#ec4c3c;background-color:#F8F9F9'> D </strong> evelepment and 
                       <strong style='color:#ec4c3c;background-color:#F8F9F9'> C </strong> linical 
                       <strong style='color:#ec4c3c;background-color:#F8F9F9'> P </strong>harmacology"), align = "right")
        )
    )
)
