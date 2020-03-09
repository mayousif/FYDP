ui = fluidPage(title = "HX Modelling", shinyjs::useShinyjs(),
  fluidPage(id = "tabset",
    #---------------------------------------------- 
    # UI code
    #----------------------------------------------
    tabPanel("Start",setBackgroundColor("#E2E2E2"),
             fluidRow(
               column(12, align = "center",
                      titlePanel(strong("Heat Exchanger Simulator"))
               )
             ),
             fluidRow(
               column(4,
                      numericInput("Fs", HTML("Shell Flowrate (kg/s)"),"", value = 0.5),
                      numericInput("Ft", HTML("Tube Flowrate (kg/s)"),"", value = 1),
                      numericInput("Tsi", HTML("Inlet Shell Temperature (K)"),"", value = 373.5),
                      numericInput("P", HTML("Shell Pressure (Pa)"),"", value = 101325),
                      numericInput("Tti", HTML("Inlet Tube Temperature (K)"),"", value = 283.4),
                      numericInput("RH", HTML("Inlet Relative Humidity (%)"),"", value = 75),
                      selectInput("fluid", HTML("Inlet Shell Fluid"),choices = c("Liquid Water" = "lw", "Vapour Water" = "vw"),selected = "vw")
                      #numericInput("Cps", HTML("Shell Heat Capacity (J/(kg K)"),"", value = 4184),
                      #numericInput("Cpt", HTML("Tube Heat Capacity (J/(kg K)"),"",value = 4184),
                      #numericInput("hs", HTML("Shell Heat Transfer Coefficient (W/(m<sup>2</sup> K)"),"", value = 1000),
                      #numericInput("ht", HTML("Tube Heat Transfer Coefficient (W/(m<sup>2</sup> K)"),"", value = 1000)

               ),
               column(4, align = "center",
                      numericInput("tid", HTML("Tube Inner Diameter (m)"),"", value = 0.02523),
                      numericInput("tt", HTML("Tube Thickness (m)"),"", value = 0.001),
                      numericInput("sid", HTML("Shell Inner Diameter (m)"),"", value = 0.1),
                      numericInput("st", HTML("Shell Thickness (m)"),"", value = 1),
                      numericInput("L", HTML("Pipe Length (m)"),"", value = 1),
                      radioButtons("config", label = HTML("Tube configuration"),
                                   choices = list("In-line" = "square", "Staggered" = "triangle"), 
                                   selected = "square"),
                      textInput("ntubesrow", HTML("Number of Tubes Per Row"),"1,1", placeholder = "1,3,5,3,1")
               ),
               column(4, align = "right",
                      numericInput("modelsteps", HTML("# of Model Steps"),"", value = 4),
                      numericInput("baffles", HTML("# of Baffles"),"",value = 1)
               )
             ),
             fluidRow(
               column(4,align = "center",
                      plotlyOutput("plot2",width = "500px", height = "500px")
               ),
               column(8,align = "center", 
                      actionButton("start", "Start"),
                      plotlyOutput("plot1",height = "440px")
               )
             ),
             fluidRow(
               column(4,align = "center", 
                      sliderInput("slider", min = NA, max = NA, value = NA, label = "",round = -3, width = "550px")
                      
               )
             )
    )
  )
)
