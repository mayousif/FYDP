# Location that program was installed to
MainLocation = getwd()

# Location of run.R, ui.R, and server.R (must also be changed in batch file)
AppLocation = paste(MainLocation,"/shiny", sep = "")

# Location of R packages library
.libPaths(paste(MainLocation,"/library", sep = ""))
LibraryLocation = paste(MainLocation,"/library", sep = "")

# Install/load required R packages
usePackage = function(p) {
  if (!is.element(p, installed.packages(lib.loc = LibraryLocation)[,1]))
    install.packages(p, dep = TRUE,lib = LibraryLocation,repos = c("http://cran.us.r-project.org","http://r-forge.r-project.org"))
  require(p, character.only = TRUE,lib.loc = LibraryLocation)
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