ui = fluidPage(title = "HX Modelling", shinyjs::useShinyjs(),
  fluidPage(id = "tabset",
    #---------------------------------------------- 
    # UI code
    #----------------------------------------------
    tabsetPanel(id = "tabs",
      tabPanel("Input",setBackgroundColor("#E2E2E2"),
                fluidRow(
                  column(12, align = "center",
                         br()
                  )
               ),
               fluidRow(
                 column(3,offset = 0.5, align = "center", style = "height:675px;background-color:#f5f5f5;border-color:#E0E0E0;border-style:inset;",
                        titlePanel(strong("Fluid Properties")),
                        numericInput("Fs", HTML("Shell Flowrate (kg/s)"),"", value = 0.0012708),
                        numericInput("Ft", HTML("Tube Flowrate (kg/s)"),"", value = 0.014833333),
                        numericInput("Tsi", HTML("Inlet Shell Temperature (K)"),"", value = 373.15),
                        numericInput("Tti", HTML("Inlet Tube Temperature (K)"),"", value = 281.15),
                        numericInput("P", HTML("Shell Pressure (Pa)"),"", value = 101325),
                        numericInput("RH", HTML("Inlet Relative Humidity (%)"),"", value = 75),
                        selectInput("fluid", HTML("Inlet Shell Fluid"),choices = c("Liquid Water" = "lw", "Vapour Water" = "vw",
                                                                                   "Vapour Water + Ethanol" = "vwe"),selected = "vw"),
                        numericInput("fracw", HTML("Inlet Water Vapour Mass Fraction"),"", value = NA)
                 ),
                 column(1),
                 column(4, align = "center", style = "height:675px;background-color:#f5f5f5;border-color:#E0E0E0;border-style:inset;",
                        titlePanel(strong("Dimensions")),
                        numericInput("tid", HTML("Tube Inner Diameter (m)"),"", value = 0.014453),
                        numericInput("tt", HTML("Tube Thickness (m)"),"", value = 0.001422),
                        numericInput("sid", HTML("Shell Inner Diameter (m)"),"", value = 0.1587),
                        numericInput("st", HTML("Shell Thickness (m)"),"", value = 0.0096),
                        numericInput("L", HTML("Pipe Length (m)"),"", value = 1.2),
                        radioButtons("config", label = HTML("Tube configuration"), inline = T,
                                     choices = list("In-line" = "square", "Staggered" = "triangle"), 
                                     selected = "triangle"),
                        textInput("ntubesrow", HTML("Number of Tubes Per Row"),"2,5,6,7,6,5,2", placeholder = "1,3,5,3,1"),
                        numericInput("baffles", HTML("# of Baffles"),"",value = 11)
                 ),
                 column(1),
                 column(3, align = "center", style = "height:675px;background-color:#f5f5f5;border-color:#E0E0E0;border-style:inset;",
                        titlePanel(strong("Model Parameters")),
                        numericInput("modelsteps", HTML("# of Model Steps"), value = 24),
                        numericInput("error", HTML("Maximum Acceptable Error (%)"), value = 0.01),
                        numericInput("iterations", HTML("Maximum Acceptable # of Iterations"), value = 50),
                        
                        actionButton("start", "Calculate")
                        
                 )
               )
      ),
      tabPanel("Output",setBackgroundColor("#E2E2E2"),
               fluidRow(
                 column(4,align = "center",
                        plotlyOutput("plot2",width = "500px", height = "500px")
                 ),
                 column(8,align = "center",
                        plotlyOutput("plot1",height = "500px")
                 )
               ),
               fluidRow(
                 column(4,align = "center", 
                        sliderInput("slider", min = NA, max = NA, value = NA, label = "", width = "550px", ticks = F, round = -3)
                        
                 ),
                 column(8, align = "center",
                        tableOutput('resulttable')
                       
                 )
               )
      )
    )
  )
)
