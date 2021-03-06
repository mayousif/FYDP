# Location that program was installed to
MainLocation = getwd()

# Location of run.R, ui.R, and server.R (must also be changed in batch file)
AppLocation = paste(MainLocation,"/shiny", sep = "")

# Location of R packages library
.libPaths("C:/Users/megue/Documents/R/win-library/3.5")
LibraryLocation = paste(MainLocation,"/library", sep = "")

# Install/load required R packages
usePackage = function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE,repos = c("http://cran.us.r-project.org","http://r-forge.r-project.org"))
  require(p, character.only = TRUE)
}
usePackage("shiny")
usePackage("shinyWidgets")
usePackage("shinyjs")
usePackage("ggplot2")
usePackage("plotly")
usePackage("deSolve")
usePackage("bvpSolve")
usePackage("pracma")
usePackage("nleqslv")

# Launch App
runApp(appDir = AppLocation,launch.browser = T)