library(shiny)
library(shinythemes)
library(SuscTesting)
library(BayesianMonoErrorModels)
library(shinycssloaders)


# library(devtools)
# install_github('gdepalma/SuscTesting')


source('ui.R')
source('server.R')

shinyApp(ui, server)
