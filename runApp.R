library(shiny)
library(shinythemes)
library(SuscTesting)
library(BayesianMonoErrorModels)

source('FunctionsERB.R')
source('FunctionsModel.R')
source('FunctionsLogistic.R')
source('FunctionsSpline.R')


# library(devtools)
# install_github('gdepalma/SuscTesting')
# 

source('ui.R')
source('server.R')

shinyApp(ui, server)
