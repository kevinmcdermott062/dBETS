library(shiny)


fluidPage(theme=shinythemes::shinytheme("cosmo"),
  
    
  # Application title
  h3("dBETS - diffusion Breakpoint Estimation Testing Software"),
  
  sidebarPanel(
    selectInput("Page", strong("Navigate"),list('Data Entry'=1,'Error-rate Bounded Method'=2,
        'Model-Based Approaches'=3),selected=1),
    # tags$style(type="text/css", '#Page { width: 70%;}'),
    # tags$style(".span2 {width: 100px; float: left; margin: 5px;}"),
                     
    ### Load Data
    conditionalPanel(condition="input.Page==1",
      HTML("<a class='btn' href='https://dbets.shinyapps.io/dBETS/'><b>Refresh / Clear</b></a>"),br(),
      wellPanel(                
      checkboxInput(inputId= "dataSource", label="Use a file stored on my local machine.", value=FALSE),
      conditionalPanel(condition = "input.dataSource == false",
        textInput(inputId="url", label="File URL:", value="https://dbets.shinyapps.io/dBETS/data1.csv")
        # tags$style(type='text/css', "#url { width: 95%; }"),
      ),
      conditionalPanel(condition = "input.dataSource == true",            
        fileInput(inputId = "file", label="")),
      selectInput(inputId="convertToLog", label="Are MIC values on log2 scale?",choices=list('Yes'=1,'No'=0),selected=1),
      conditionalPanel(condition = "input.oneBrkpt == false",
        p('MIC Breakpoints (log2 scale)'),
        fluidRow(
          column(6,selectInput("MICBrkptL", "Lower (<=)", choices =seq(-8,8,by=1),selected=-1)),
          column(6,selectInput("MICBrkptU", "Upper (>=)", choices =seq(-8,8,by=1),selected=1))
        )),
      conditionalPanel(condition = "input.oneBrkpt == true",
        br(),
        selectInput("MICBrkpt", "MIC Breakpoint Log2 Scale (<=)", choices =seq(-8,8,by=1),selected=0)),
      checkboxInput(inputId= "oneBrkpt", label="One MIC Breakpoint", value=FALSE),
      br(),
      actionButton('downloadData',strong('Load Data and Get Started')),
      tags$style(type='text/css', "#downloadData { width: 100%; }"),
      br(),
      conditionalPanel(condition = "input.downloadData!=0",
        br(),div("Graph Options: ", style="font-weight: bold"),br(),
        fluidRow(
          column(6,selectInput("FlipS", "X-Axis:",list('MIC'=TRUE,'DIA'=FALSE),selected=TRUE)),
          column(6,selectInput(inputId='miclogS',label='Log2 Labels',choices=list('Yes'=TRUE,'No'=FALSE),selected=TRUE))
        ),
      conditionalPanel(condition="input.downloadData!=0",
        br(),actionButton('startGraph',strong('Plot Graph'))),
      conditionalPanel(condition="input.startGraph!=0",
        br(),downloadButton('downloadStart',strong('Download Graph'))))
    )),
      
    ### ERB
    conditionalPanel(condition="input.Page==2 & input.panel !=1",
      br(),
      wellPanel(
        tags$style(".span2 {width: 100px; float: left; margin: 5px;}"),
        
        ### first panel
        conditionalPanel(condition="input.panel==2",
          conditionalPanel(condition = "input.oneBrkpt == false",
             fluidRow(
               column(5,selectInput("minWidth1", "Min Width:",choices =seq(2,8,by=1),selected=3)),
               column(5,selectInput("maxWidth1", "Max Width:",choices =seq(3,20,by=1),selected=10))
            ),
            br(),
            div("Error Percentages within one of Intermediate Range"),
            fluidRow(
                column(4,selectInput("VM11", "VM:",choices =seq(0,30,by=.5),selected=10)),
                column(4,selectInput("M11", "Major:",choices =seq(0,30,by=.5),selected=10)),
                column(4,selectInput("m11", "Minor:",choices =seq(0,60,by=.5),selected=40))
            ),
            br(),
            div("Error Percentages outside one of Intermediate Range"),
            fluidRow(
              column(4,selectInput("VM21", "VM:",choices =seq(0,10,by=.5),selected=2)),
              column(4,selectInput("M21", "Major:",choices =seq(0,10,by=.5),selected=2)),
              column(4,selectInput("m21", "Minor:",choices =seq(0,20,by=.5),selected=5))
            ),
            br()),
          conditionalPanel(condition = "input.oneBrkpt == true",
            fluidRow(
              column(6,selectInput("VMOneBrkpt", "VM:",choices =seq(0,30,by=.5),selected=1)),
              column(6,selectInput("MOneBrkpt", "Major:",choices =seq(0,30,by=.5),selected=5))
            ),br()),
          actionButton('actionERB',strong('Run')),br(),br(),
          tags$style(type="text/css", '#actionERB { width: 100%;}'),
          conditionalPanel(condition="input.actionERB!=0",
             div("Graph Options: ", style="font-weight: bold"),br(),          
             fluidRow(
               column(5,selectInput("FlipERB1", "X-Axis:",list('MIC'=TRUE,'DIA'=FALSE),selected=TRUE)),
               column(5,selectInput('miclogERB1', 'Log2 Labels',list('Yes'=TRUE,'No'=FALSE),selected=TRUE))
             ),br(),
             actionButton('plotERB1',strong('Plot ERB Graph')), br(),
             conditionalPanel(condition="input.plotERB1!=0",
                br(),downloadButton('downloadERB',strong('Download ERB Graph'))))),     
        
        ### second panel
        conditionalPanel(condition="input.panel==3",
          conditionalPanel(condition = "input.oneBrkpt == false",
             fluidRow(
                 column(5,uiOutput("minWidth2")),
                 column(5,uiOutput("maxWidth2"))
             ),
             br(),
             div("Error Percentages within one of Intermediate Range"),
            fluidRow(
              column(4,selectInput("VM12", "VM:",choices =seq(0,30,by=.5),selected=10)),
              column(4,selectInput("M12", "Major:",choices =seq(0,30,by=.5),selected=10)),
              column(4,selectInput("m12", "Minor:",choices =seq(0,60,by=.5),selected=40))
            ),
            br(),
            div("Error Percentages outside one of Intermediate Range"),
            fluidRow(
              column(4,selectInput("VM22", "VM:",choices =seq(0,30,by=.5),selected=10)),
              column(4,selectInput("M22", "Major:",choices =seq(0,30,by=.5),selected=10)),
              column(4,selectInput("m22", "Minor:",choices =seq(0,60,by=.5),selected=40))
            )),
            conditionalPanel(condition = "input.oneBrkpt == true",
              fluidRow(
                column(5,selectInput("VMOneBrkpt2", "VM:",choices =seq(0,30,by=.5),selected=1)),
                column(5,selectInput("MOneBrkpt2", "Major:",choices =seq(0,30,by=.5),selected=5))
               )),
            br(),actionButton('actionBoot',strong('Run')),
            # tags$style(type="text/css", '#actionBoot { width: 100%;}'),
            conditionalPanel(condition="input.actionBoot!=0  & input.oneBrkpt == false",
              br(),actionButton('plotBootERB',strong('Plot Distribution DIA Breakpoints')))),
        
        ###third panel
        conditionalPanel(condition="input.panel==4",
          conditionalPanel(condition = "input.oneBrkpt == false",
            fluidRow(
                column(5,uiOutput("D1ERB")),
                column(5,uiOutput("D2ERB"))
            ),
            br(),
            div("Error Percentages within one of Intermediate Range"),
            fluidRow(
              column(4,selectInput("VM13", "VM:",choices =seq(0,30,by=.5),selected=10)),
              column(4,selectInput("M13", "Major:",choices =seq(0,30,by=.5),selected=10)),
              column(4,selectInput("m13", "Minor:",choices =seq(0,60,by=.5),selected=40))
            ),
            br(),
            div("Error Percentages outside one of Intermediate Range"),
            fluidRow(
              column(4,selectInput("VM23", "VM:",choices =seq(0,30,by=.5),selected=10)),
              column(4,selectInput("M23", "Major:",choices =seq(0,30,by=.5),selected=10)),
              column(4,selectInput("m23", "Minor:",choices =seq(0,60,by=.5),selected=40))
            )),
          conditionalPanel(condition = "input.oneBrkpt == true",
            HTML('<table border=0 width="100%"><tr bgcolor="#f5f5f5"><td>'),
            uiOutput("D11ERB"),
            HTML('</td></tr></table>'),
            tags$style(type="text/css", '#D11ERB { width: 50%;}'),             
            br(),
            div("Error Percentages", style="font-weight: bold"),
            fluidRow(
              column(5,selectInput("VMOneBrkpt3", "VM:",choices =seq(0,30,by=.5),selected=1)),
              column(5,selectInput("MOneBrkpt3", "Major:",choices =seq(0,30,by=.5),selected=5))
            )),
          br(),actionButton('actionERBSelected',strong('Run')),
          # tags$style(type="text/css", '#actionERBSelected { width: 100%;}'),
          conditionalPanel(condition="input.actionERBSelected!=0",
            br(),div("Graph Options: ", style="font-weight: bold"),br(),
            fluidRow(
              column(5,selectInput("FlipERB3", "X-Axis:",list('MIC'=TRUE,'DIA'=FALSE),selected=TRUE)),
              column(5,selectInput('miclogERB3', 'Log2 Labels',list('Yes'=TRUE,'No'=FALSE),selected=TRUE))
            ),br(),
            actionButton('plotERBSelected',strong('Plot ERB Selected Graph')), br(),
            conditionalPanel(condition="input.plotERBSelected!=0",
              br(),downloadButton('downloadERBSelected',strong('Download ERB Selected Graph')))))                   
    )),
    
    ### Model
    conditionalPanel(condition="input.Page==3 & input.panel2 !=1",
       br(),
       wellPanel(
         tags$style(".span2 {width: 100px; float: left; margin: 5px;}"),
         
         #### first panel
          conditionalPanel(condition="input.panel2==2",
            conditionalPanel(condition = "input.oneBrkpt == false",             
               fluidRow(
                 column(5,selectInput("minWidthM1", "Min Width:",choices =seq(2,8,by=1),selected=3)),
                 column(5,selectInput("maxWidthM1", "Max Width:",choices =seq(3,20,by=1),selected=7))
               ),br()),
            selectInput("lossType1", "Loss Function:",list('Squared'=1),selected=1),br(),
            actionButton('actionLog',strong('Run')),br(),br(),
            # tags$style(type="text/css", '#actionLog { width: 100%;}'),
          conditionalPanel(condition="input.actionLog!=0",
            # br(),div("Graph Options: ", style="font-weight: bold"),br(),
            # fluidRow(
            #   column(5,selectInput("FlipMod1", "X-Axis:",list('MIC'='Yes','DIA'='No'),selected='No')),
            #   column(5,selectInput('miclogM1', 'Log2 Labels',list('Yes'=1,'No'=0),selected=1))
            # ),
            # br(),
            actionButton('actionLogPlot',strong('Plot Logistic Fit Graph')))),
          # conditionalPanel(condition="input.actionLogPlot!=0",
            # br(),downloadButton('downloadLog',strong('Download Logistic Graph')))),
          # tags$style(type="text/css", '#actionLog { width: 90%;}'),
          # tags$style(type="text/css", '#actionLogPlot { width: 90%;}'),
             
         ### Second panel
          conditionalPanel(condition="input.panel2==3",
            conditionalPanel(condition = "input.oneBrkpt == false",
               fluidRow(
                 column(5,selectInput("minWidthM2", "Min Width:",choices =seq(2,8,by=1),selected=3)),
                column(5,selectInput("maxWidthM2", "Max Width:",choices =seq(3,20,by=1),selected=7))
               ),br()),
            
            selectInput("lossType2", "Loss Function:",list('Squared'=1),selected=1),br(),
            actionButton('actionSpline',strong('Run')),br(),br(),
            # tags$style(type="text/css", '#actionSpline { width: 100%;}'),
          conditionalPanel(condition="input.actionSpline!=0",
            # br(),
            # fluidRow(
            #   column(5,selectInput("FlipMod2", "X-Axis:",list('MIC'='Yes','DIA'='No'),selected='No')),
            #   column(5,selectInput('miclogM2', 'Log2 Labels',list('Yes'=1,'No'=0),selected=1))
            # ),
            # br(),
            actionButton('actionSplinePlot',strong('Plot Spline Fit Graph')))),
          # conditionalPanel(condition="input.actionSplinePlot!=0",
            # br(),downloadButton('downloadSpline',strong('Download Spline Graph')))),
          # tags$style(type="text/css", '#actionSpline { width: 90%;}'),
          # tags$style(type="text/css", '#actionSplinePlot { width: 90%;}'),
         
         ### Third Panel
         conditionalPanel(condition="input.panel2==4",
              # div(class="row-fluid",
              #     div(class="span2",selectInput("FlipMod3", "X-Axis:",list('MIC'='Yes','DIA'='No'),selected='No')),
              #     div(class="span2",selectInput('miclogM3', 'Log2 Labels',list('Yes'=1,'No'=0),selected=1))
              # ),
              # br(),
              actionButton('actionCompare',strong('Plot Compare Fit Graph'))),
              # conditionalPanel(condition="input.actionCompare!=0",
                # br(),downloadButton('downloadCompare',strong('Download Graph')))),
              # tags$style(type="text/css", '#actionCompare { width: 90%;}'),
         
         ### Fourth Panel
         conditionalPanel(condition="input.panel2==5",
          selectInput("modelSelect", "Model:",list('Spline'='spline','Logistic'='logistic'),selected='spline'),
          conditionalPanel(condition = "input.oneBrkpt == false",
             fluidRow(
               column(5,uiOutput("D1Set1")),
               column(5,uiOutput("D2Set1"))
             ),
            checkboxInput(inputId='includeSecondDIA', label='Include Second Set of DIA Breakpoints?',value=FALSE),
            conditionalPanel(condition="input.includeSecondDIA!=0",
              fluidRow(
                column(5,uiOutput("D1Set2")),
                column(5,uiOutput("D2Set2")),br()
               ))),
          conditionalPanel(condition = "input.oneBrkpt == true",
            uiOutput("D1One"),
            checkboxInput(inputId='includeSecondDIAOne', label='Include Second DIA Breakpoint?',value=FALSE),
            conditionalPanel(condition="input.includeSecondDIAOne!=0",
              uiOutput("D2One"))),
          br(),
          div("Graph Options: ", style="font-weight: bold"),
          br(),selectInput('miclogM4', 'Log2 Labels',list('Yes'=1,'No'=0),selected=1),
          tags$style(type="text/css", '#miclogM4 { width: 40%;}'),br(),br(),
          actionButton('actionProbDIA',strong('Run')),
          tags$style(type="text/css", '#actionProbDIA { width: 100%;}'))
          # conditionalPanel(condition="input.actionProbDIA!=0",
          # br(),downloadButton('downloadProbDIA',strong('Download Graph')))),
          # tags$style(type="text/css", '#actionProbDIA { width: 90%;}')
    ))
  ),
      
  mainPanel(
  
  ###Load Data
  conditionalPanel(condition="input.Page==1",
    list(
      conditionalPanel("input.downloadData==0",
        br(),
        h5(strong('Welcome to dBETS software! (Version 1.5)')),
        br(),
        div(HTML("Software package created by Glen DePalma (<a href=https://github.com/gdepalma>https://github.com/gdepalma</a>) 
                 and Bruce A. Craig (<a href=http://www.stat.purdue.edu/people/faculty/bacraig>http://www.stat.purdue.edu/people/faculty/bacraig</a>), with advice from John Turnidge.")),br(),
        div(HTML("The determination of diffusion test breakpoints has become a challenging issue due to the increasing resistance of microorganisms to antibiotics.   
            dBETS (<b>d</b>iffusion <b>B</b>reakpoint <b>E</b>stimation <b>T</b>esting <b>S</b>oftware) helps clinicians  easily analyze data from susceptibility experiments through visualization, error-rate bounded, and model-based approaches.")),br(),
        div(HTML("There is a navigational panel near the top left corner to jump between the three sections: 1) Data Entry, 
            2) Error-Rate Bounded Method, and 3) Model-Based Approaches.  On each page is a left-side panel of user-specified inputs
          and options.  Once these inputs are set, a click of a run button  (located on the panel below the inputs) activates the procedure.
          The output will appear in the main section on the right once the procedure completes.  Data can be uploaded via the local hard drive or from a URL. For information on how to format data for use with dBETS 
            <a href=https://dbets.shinyapps.io/dBETS/Formatting_Data.pdf>click here.</a>")),br(),
        div(HTML("Please contact Glen DePalma, <a href=mailto:glen.depalma@gmail.com>glen.depalma@gmail.com</a>, to report problems or suggestions.  The code and helper packages are located 
                 on Glen's GitHub page. All code is under GPL-3 license (feel free to use, modify, and distribute; but source code must be made available).")),br(),
        div(HTML("<hr>")),
        div(HTML("<u>Paper describing the use of dBETS: </u>")),br(),
        div(HTML("DePalma, Glen and Turnidge, John and Craig, Bruce A.  <b>Determination of disk diffusion susceptibility testing interpretive criteria using model-based analysis: development and implementation.</b>  <i>Diagnostic Microbiology & Infectious Disease.</i> 2017. 87[2] 143-149.")),
        div(HTML("<a href=http://www.dmidjournal.com/article/S0732-8893(16)30044-X/pdf>http://www.dmidjournal.com/article/S0732-8893(16)30044-X/pdf</a>")),br(),
        div(HTML("<u>Paper describing the models behind dBETS: </u>")),br(),
        div(HTML("DePalma, Glen and Craig, Bruce A.  <b>Bayesian monotonic errors-in-variables models with applications to pathogen susceptibility testing.</b> <i>Statistics in Medicine.</i> 2018. 37[3] 487-502. DOI: 10.1002/sim.7533.")),
        div(HTML("<a href=https://www.ncbi.nlm.nih.gov/pubmed/29156492>https://www.ncbi.nlm.nih.gov/pubmed/29156492</a>")),
        div(HTML("<a href=https://arxiv.org/abs/1806.06974>https://arxiv.org/abs/1806.06974</a>")),br(),
        div(HTML("<u>Papers citing dBETS: </u>")),br(),
        div(HTML("Badger, S., Abraham, S., Saputra, S., Trott, D.J., Turnidge, J., Mitchell, T., Caraguel, C.G.B., Jordan, D. 
          <b>Relative performance of antimicrobial susceptibility assays on clinical Escherichia coli isolates from animals.</b>
          <i>Veterinary Microbiology.</i> 2018. 214 [56-64]")),
        div(HTML("<a href=https://www.ncbi.nlm.nih.gov/pubmed/29408033>https://www.ncbi.nlm.nih.gov/pubmed/29408033</a>"))
        
        
        ),
      conditionalPanel("input.downloadData>0",br(),
        verbatimTextOutput("loadData"),
        plotOutput("plotData",width = "800px", height = "600px")))),
  
  
  #\textbf{DePalma, Glen} and Craig, Bruce A.  \emph{Bayesian monotonic errors-in-variables models with applications to pathogen susceptibility testing.}  Statistics in Medicine. 2018. 37[3] 487-502. DOI: 10.1002/sim.7533.
  # \textbf{DePalma, Glen} and Turnidge, John and Craig, Bruce A.  \emph{Determination of disk diffusion susceptibility testing interpretive criteria using model-based analysis: development and implementation.}  Diagnostic Microbiology \& Infectious Disease. 2017. 87[2] 143-149.

  ###ERB
  conditionalPanel(condition="input.Page==2",
    tabsetPanel(
       tabPanel("Find optimal DIA breakpoints", list(br(),
                verbatimTextOutput("ERBCalc"),
                plotOutput("ERBPlot",width = "800px", height = "600px")),value=2),
       tabPanel("Assess uncertainty in estimate", list(br(),
                 verbatimTextOutput("bootERB"),plotOutput("bootPlot",width="700px",height="500px")),value=3),
       tabPanel("Evaluate index for selected breakpoints", list(br(),
                verbatimTextOutput("ERBCalc2"),
                plotOutput("ERBPlot2",width = "800px", height = "600px")),value=4),
       id='panel'
    )),
  
  ###Model
  conditionalPanel(condition="input.Page==3",
       tabsetPanel(
         tabPanel("Logistic Model", list(p(strong("This function will take several minutes to run.")),
                  br(),verbatimTextOutput("logistic"),br(),br(),
                  plotOutput("logisticPlot",width = "800px", height = "600px")),value=2),
         tabPanel("Spline Model", list(p(strong("This function will take several minutes to run.")),
                   br(),
                   verbatimTextOutput("spline"),
                   plotOutput("splinePlot",width = "800px", height = "600px"),
                   plotOutput("splineBrkptPlot",width = "700px", height = "500px")),value=3),
         tabPanel("Compare Model Fits", list(br(),
                   verbatimTextOutput("compareFits"),
                   fluidRow(
                     splitLayout(cellWidths = c("50%", "50%"), verbatimTextOutput("CompareDescLog"),
                                 verbatimTextOutput("CompareDescSpline"))
                   ),
                   plotOutput("compareFitsPlot",width = "800px", height = "600px")),value=4),
         tabPanel("DIA Classification Probability", list(
                 p(strong("MIC breakpoint(s) shifted left 0.5 units to account for rounding.")),
                 br(),                                
                 verbatimTextOutput("probDIA"),
                 plotOutput("probDIAPlot",width = "800px", height = "600px")),value=5),
         id='panel2'
    ))
  )
)



