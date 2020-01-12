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
               column(4, align = "right",
                      numericInput("Fs", HTML("Shell Flowrate (kg/s)"),"", value = 0.0498),
                      numericInput("Ft", HTML("Tube Flowrate (kg/s)"),"", value = 0.0265),
                      numericInput("Tsi", HTML("Inlet Shell Temperature (K)"),"", value = 322.3),
                      numericInput("Tti", HTML("Inlet Tube Temperature (K)"),"", value = 283.4),
                      numericInput("Cps", HTML("Shell Heat Capacity (J/(kg K)"),"", value = 4184),
                      numericInput("Cpt", HTML("Tube Heat Capacity (J/(kg K)"),"",value = 4184),
                      numericInput("hs", HTML("Shell Heat Transfer Coefficient (W/(m<sup>2</sup> K)"),"", value = 1000),
                      numericInput("ht", HTML("Tube Heat Transfer Coefficient (W/(m<sup>2</sup> K)"),"", value = 1000)

               ),
               column(4, align = "center",
                      numericInput("tid", HTML("Tube Inner Diameter (m)"),"", value = 0.01),
                      numericInput("tt", HTML("Tube Thickness (m)"),"", value = 0.001),
                      numericInput("sid", HTML("Shell Inner Diameter (m)"),"", value = 1),
                      numericInput("st", HTML("Shell Thickness (m)"),"", value = 1),
                      numericInput("L", HTML("Pipe Length (m)"),"", value = 10),
                      numericInput("ntubes", HTML("Number of Tubes"),"", value = 1),
                      textInput("ntubesrow", HTML("Number of Tubes Per Row"),"", placeholder = "1,3,5,3,1")
               ),
               column(4, align = "left",
                      numericInput("modelsteps", HTML("# of Model Steps"),"", value = 100),
                      numericInput("baffles", HTML("# of Baffles"),"",value = 4)
               )
             ),
             fluidRow(
               column(12,align = "center",
                      actionButton("start", "Start"),
                      plotlyOutput("plot1")
               )
             )
    )
  )
)
