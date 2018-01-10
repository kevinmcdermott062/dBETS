library(shiny)
library(shinythemes)
library(SuscTesting)
library(BayesianMonoErrorModels)
library(shinycssloaders)


# library(devtools)
# install_github('gdepalma/SuscTesting')
# devtools::install_github('andrewsali/shinycssloaders')


source('ui.R')
source('server.R')

shinyApp(ui, server)
