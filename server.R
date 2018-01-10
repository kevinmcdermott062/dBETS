
server <- function(input, output, session){
  
  ### Reactive UI functions
  #ERB Panel 2
  output$minWidth2 <- renderUI({
    selectInput("minWidth2", "Min Width:",choices =seq(2,8,by=1),selected=as.numeric(input$minWidth1))
  })
  output$maxWidth2 <- renderUI({
    selectInput("maxWidth2", "Max Width:",choices =seq(3,20,by=1),selected=as.numeric(input$maxWidth1))
  })
  output$VM12 <- renderUI({
    selectInput("VM12", "Very Major:",choices =seq(0,30,by=.5),selected=as.numeric(input$VM11))
  })
  output$M12 <- renderUI({
    selectInput("M12", "Major:",choices =seq(0,30,by=.5),selected=as.numeric(input$M11))
  })
  output$m12 <- renderUI({
    selectInput("m12", "Minor:",choices =seq(0,60,by=.5),selected=as.numeric(input$m11))
  })
  output$VM22 <- renderUI({
    selectInput("VM22", "Very Major:",choices =seq(0,30,by=.5),selected=as.numeric(input$VM21))
  })
  output$M22 <- renderUI({
    selectInput("M22", "Major:",choices =seq(0,30,by=.5),selected=as.numeric(input$M21))
  })
  output$m22 <- renderUI({
    selectInput("m22", "Minor:",choices =seq(0,60,by=.5),selected=as.numeric(input$m21))
  })
  ###ERB Panel 3
  output$VM13 <- renderUI({
    selectInput("VM13", "Very Major:",choices =seq(0,30,by=.5),selected=as.numeric(input$VM11))
  })
  output$M13 <- renderUI({
    selectInput("M13", "Major:",choices =seq(0,30,by=.5),selected=as.numeric(input$M11))
  })
  output$m13 <- renderUI({
    selectInput("m13", "Minor:",choices =seq(0,60,by=.5),selected=as.numeric(input$m11))
  })
  output$VM23 <- renderUI({
    selectInput("VM23", "Very Major:",choices =seq(0,30,by=.5),selected=as.numeric(input$VM21))
  })
  output$M23 <- renderUI({
    selectInput("M23", "Major:",choices =seq(0,30,by=.5),selected=as.numeric(input$M21))
  })
  output$m23 <- renderUI({
    selectInput("m23", "Minor:",choices =seq(0,60,by=.5),selected=as.numeric(input$m21))
  })
  
  
  output$loadData <- renderPrint({
    if(input$downloadData!=0){
      return(isolate({
        
        MICLog2 <<- as.numeric(input$convertToLog)

        if (input$dataSource == FALSE){            
          url=input$url
          if(url=='https://dbets.shinyapps.io/dBETS/data1.csv')
            url='data1.csv'
          if(grepl('https',url)==TRUE) stop ('https:// not supported, please use http://')
          a1=read_csv(url)
        } else{      
          df <- input$file
          path <- df$datapath  
          a1 <- read_csv(path)
        }
        parms=parse_data(a1)
          
        xcens=rep(0,nrow(a1))
        ycens=rep(0,nrow(a1))
        MIC=rep(NA,nrow(a1))
        DIA=rep(NA,nrow(a1))
        
        #check for censoring
        for(i in 1:nrow(a1)){
          if(grepl('>',a1[i,1])==TRUE){
            MIC[i]=substr(a1[i,1],2,nchar(a1[i,1]))
            xcens[i]=1
          }else if(grepl('<',a1[i,1])==TRUE){
            MIC[i]=substr(a1[i,1],2,nchar(a1[i,1]))
            xcens[i]=-1
          }else
            MIC[i]=a1[i,1]
          
          if(grepl('>',a1[i,2])==TRUE){
            DIA[i]=substr(a1[i,2],2,nchar(a1[i,2]))
            ycens[i]=1
          }else if(grepl('<',a1[i,2])==TRUE){
            DIA[i]=substr(a1[i,2],2,nchar(a1[i,2]))
            ycens[i]=-1
          }else
            DIA[i]=a1[i,2]
        }
        
        MIC=parms$MIC
        
        if(MICLog2==0)
          MIC<-log(MIC,2)
        
        if(sum(is.na(MIC))>1) stop('Problem converting MIC to log2 Scale!')
        
        MIC<<- MIC
        DIA<<-parms$DIA
        xcens<<-parms$xcens
        ycens<<-parms$ycens
        
        if(input$oneBrkpt==FALSE){
          MICBrkptL <<- as.numeric(input$MICBrkptL)
          MICBrkptU <<- as.numeric(input$MICBrkptU)
          if(MICBrkptL>=MICBrkptU) stop('Lower MIC Breakpoint must be smaller then Upper MIC Breakpoint')
          descriptiveStat(MIC,DIA,xcens,ycens,MICBrkptL,MICBrkptU)
        }else{
          MICBrkpt <<- as.numeric(input$MICBrkpt)
          descriptiveStatOne(MIC,DIA,xcens,ycens,MICBrkpt)
        }
                
        
        
        xgrid <<- seq(min(MIC)-1,max(MIC)+1,length=1000)
        
        ### ERB Panel 3
        output$D1ERB <- renderUI({
          selectInput("D1ERB", "Lower DIA Breakpoint:",choices =seq(6,60,by=1),selected=round(quantile(DIA,probs=.4)))
        })
        output$D2ERB <- renderUI({
          selectInput("D2ERB", "Upper DIA Breakpoint:",choices =seq(6,60,by=1),selected=round(quantile(DIA,probs=.6)))
        })
        output$D11ERB <- renderUI({
          selectInput("D11ERB", "DIA Breakpoint:",choices =seq(6,60,by=1),selected=round(quantile(DIA,probs=.4)))
        })
        
        ### Model Panel 4
        output$D1Set1 <- renderUI({
          selectInput("D1Set1", "Set 1 Lower DIA Brkpt:",choices =seq(6,50,by=1),selected=round(quantile(DIA,probs=.4)))
        })
        output$D2Set1 <- renderUI({
          selectInput('D2Set1', 'Set 1 Upper DIA Brkpt:',choices =seq(6,50,by=1),selected=round(quantile(DIA,probs=.6)))
        })
        
        output$D1Set2 <- renderUI({
          selectInput("D1Set2", "Set 2 Lower DIA Brkpt:",choices =seq(6,50,by=1),selected=round(quantile(DIA,probs=.45)))
        })
        output$D2Set2 <- renderUI({
          selectInput("D2Set2", "Set 2 Lower DIA Brkpt:",choices =seq(6,50,by=1),selected=round(quantile(DIA,probs=.65)))
        })
        output$D1One <- renderUI({
          selectInput("D1One", "DIABrkpt 1:",choices =seq(6,50,by=1),selected=round(quantile(DIA,probs=.4)))
        })
        output$D2One <- renderUI({
          selectInput('D2One', 'DIABrkpt 2:',choices =seq(6,50,by=1),selected=round(quantile(DIA,probs=.5)))
        })
              
        invisible()
      }
    ))
    }else{return(invisible())}
  })  
  
  output$plotData <- renderPlot({
    if(input$startGraph!=0){
      return(isolate({
        print(input$miclogS)
        print(input$FlipS)
        if(input$oneBrkpt==FALSE)
          fit=basicPlot(MIC,DIA,xcens,ycens,MICBrkptL,MICBrkptU,MICXaxis=input$FlipS,log2MIC=input$miclogS)
        else
          fit=basicPlotOne(MIC,DIA,xcens,ycens,MICBrkpt,MICXaxis=input$FlipS,log2MIC=input$miclogS)
        plot1 <<- fit
        plot(fit)
        invisible()
      }
    ))
    }else{return(NULL)}
  })  
  
  output$downloadStart <- downloadHandler(
    filename=function(){paste('BasicPlot', Sys.Date(), '.png', sep='')},
    content=function(file){
      png(file,width=1000,height=800,units='px')
      print(plot1)
      dev.off()
    }
  )  
      
  ### ERB Panel 1
  output$ERBCalc <- renderPrint({
    if(input$actionERB!=0){
      return(isolate({
        if(input$panel==2){
          if(exists('MIC')==FALSE) stop('Error!  Was data loaded?')
          if(input$oneBrkpt==FALSE){
            VM1=as.numeric(input$VM11)
            M1=as.numeric(input$M11)
            m1=as.numeric(input$m11)
            VM2=as.numeric(input$VM21)
            M2=as.numeric(input$M21)
            m2=as.numeric(input$m21)
            minWidth=as.numeric(input$minWidth1)
            maxWidth=as.numeric(input$maxWidth1)
            parms=findBrkptsERB(MIC,DIA,MICBrkptL=MICBrkptL,MICBrkptU=MICBrkptU,
                          VM1=VM1,M1=M1,m1=m1,VM2=VM2,M2=M2,m2=m2,
                          minWidth=minWidth,maxWidth=maxWidth)
            output$D1ERB <- renderUI({
              selectInput("D1ERB", "Lower DIA Breakpoint:",choices =seq(6,60,by=1),selected=parms$D1)}
            )
            output$D2ERB <- renderUI({
              selectInput("D2ERB", "Upper DIA Breakpoint:",choices =seq(6,60,by=1),selected=parms$D2)}
            )
          }else{
            VM=as.numeric(input$VMOneBrkpt)
            M=as.numeric(input$MOneBrkpt)
            parms=findBrkptsERBOne(MIC,DIA,VM,M,MICBrkpt)
            output$D11ERB <- renderUI({
              selectInput("D11ERB", "DIA Breakpoint:",choices =seq(6,60,by=1),selected=parms$DIABrkpt)}
            )
          }
          invisible()
        }
      }))
    }else{return(invisible())}
  })
  
  output$ERBPlot <- renderPlot({
    if(input$plotERB1!=0){
      return(isolate({
        if(input$panel==2){
          if(input$oneBrkpt==FALSE){
            VM1=as.numeric(input$VM11)
            M1=as.numeric(input$M11)
            m1=as.numeric(input$m11)
            VM2=as.numeric(input$VM21)
            M2=as.numeric(input$M21)
            m2=as.numeric(input$m21)
            minWidth=as.numeric(input$minWidth1)
            maxWidth=as.numeric(input$maxWidth1)
            convert=FALSE
            if(input$miclogE1=='1') convert=TRUE
            fit=PlotBrkptsERB2(MIC=MIC,DIA=DIA,xcens=xcens,ycens=ycens,MICBrkptL=MICBrkptL,MICBrkptU=MICBrkptU,
                           VM1=VM1,M1=M1,m1=m1,VM2=VM2,M2=M2,m2=m2,
                           minWidth=minWidth,maxWidth=maxWidth,flipGraph=input$FlipERB1,logConvert=convert)
          }else{
            VM=as.numeric(input$VMOneBrkpt)
            M=as.numeric(input$MOneBrkpt)
            convert=FALSE
            if(input$miclogE1=='1') convert=TRUE
            fit=PlotBrkptsERB2One(MIC,DIA,xcens,ycens,VM,M,MICBrkpt,flipGraph=input$FlipERB1,logConvert=convert)
          }
          plot2 <<- fit
          plot(fit)
        }
      }))
    }else{return(invisible())}
  })
  
  output$downloadERB <- downloadHandler(
    filename=function(){paste('ERB', Sys.Date(), '.png', sep='')},
    content=function(file){
      png(file,width=1000,height=800,units='px')
      print(plot2)
      dev.off()
    }
  )  
  
  ### ERB Panel 2
  output$bootERB <- renderPrint({
    if(input$actionBoot!=0){
      return(isolate({
        if(input$panel==3){
          if(exists('MIC')==FALSE) stop('Error!  Was data loaded?')
          if(input$oneBrkpt==FALSE){
            VM1=as.numeric(input$VM12)
            M1=as.numeric(input$M12)
            m1=as.numeric(input$m12)
            VM2=as.numeric(input$VM22)
            M2=as.numeric(input$M22)
            m2=as.numeric(input$m22)
            minWidth=as.numeric(input$minWidth2)
            maxWidth=as.numeric(input$maxWidth2)
            a1=bootStrapERB(MIC,DIA,MICBrkptL=MICBrkptL,MICBrkptU=MICBrkptU,
                         VM1=VM1,M1=M1,m1=m1,VM2=VM2,M2=M2,m2=m2,
                         minWidth=minWidth,maxWidth=maxWidth,session)
            bootData<<-a1
          }else{
            VM=as.numeric(input$VMOneBrkpt2)
            M=as.numeric(input$MOneBrkpt2)
            a1=bootStrapERBOne(MIC,DIA,MICBrkpt,VM,M,session)
          }
          invisible()
        }
      }))
    }else{return(invisible())}
  })
  
    
  output$bootPlot <- renderPlot({
    if(input$plotBootERB!=0){
      return(isolate({
        if(input$panel==3){
          plotBootDataERB(bootData)
          invisible()
        }
      }))
    }else{return(invisible())}
  })
  
  ### ERB Panel 3
  output$ERBCalc2 <- renderPrint({
    if(input$actionERBSelected!=0){
      return(isolate({
        if(input$panel==4){
          if(exists('MIC')==FALSE) stop('Error!  Was data loaded?')
          if(input$oneBrkpt==FALSE){
            D1=as.numeric(input$D1ERB)
            D2=as.numeric(input$D2ERB)
            VM1=as.numeric(input$VM13)
            M1=as.numeric(input$M13)
            m1=as.numeric(input$m13)
            VM2=as.numeric(input$VM23)
            M2=as.numeric(input$M23)
            m2=as.numeric(input$m23)
            ERBGivenDIA(MIC,DIA,MICBrkptL=MICBrkptL,MICBrkptU=MICBrkptU,DIABrkptL=D1,DIABrkptU=D2,
                        VM1=VM1,M1=M1,m1=m1,VM2=VM2,M2=M2,m2=m2)
          }else{
            VM=as.numeric(input$VMOneBrkpt3)
            M=as.numeric(input$MOneBrkpt3)
            ERBGivenDIAOne(MIC,DIA,xcens,ycens,MICBrkpt,DIABrkpt=as.numeric(input$D11ERB)-.5,VM,M)
          }
          invisible()
        }
      }))
    }else{return(invisible())}
  })
  
  output$ERBPlot2 <- renderPlot({
    if(input$plotERBSelected!=0){
      return(isolate({
        if(input$panel==4){
          if(input$oneBrkpt==FALSE){
            D1=as.numeric(input$D1ERB)
            D2=as.numeric(input$D2ERB)
            VM1=as.numeric(input$VM13)
            M1=as.numeric(input$M13)
            m1=as.numeric(input$m13)
            VM2=as.numeric(input$VM23)
            M2=as.numeric(input$M23)
            m2=as.numeric(input$m23)
            convert=FALSE
            if(input$miclogE3=='1') convert=TRUE
            fit=PlotBrkptsERBGiven(MIC=MIC,DIA=DIA,xcens=xcens,ycens=ycens,MICBrkptL=MICBrkptL,
                 MICBrkptU=MICBrkptU,D1=D1,D2=D2, VM1=VM1,M1=M1,m1=m1,VM2=VM2,M2=M2,m2=m2,
                 flipGraph=input$FlipERB3,logConvert=convert)
          }else{
            VM=as.numeric(input$VMOneBrkpt3)
            M=as.numeric(input$MOneBrkpt3)
            convert=FALSE
            if(input$miclogE3=='1') convert=TRUE
            fit=PlotBrkptsERBGivenOne(MIC,DIA,xcens,ycens,VM,M,MICBrkpt,DIABrkpt=as.numeric(input$D11ERB)-.5,flipGraph=input$FlipERB3,logConvert=convert)
          }
          plot3 <<- fit
          plot(fit)
        }
      }))
    }else{return(invisible())}
  })
  
  output$downloadERBSelected <- downloadHandler(
    filename=function(){paste('ERBGivenDIA', Sys.Date(), '.png', sep='')},
    content=function(file){
      png(file,width=1000,height=800,units='px')
      print(plot3)
      dev.off()
    }
  )  

  #### Model Panel 1
  output$logistic <- renderPrint({
    if(input$actionLog!=0){
      return(isolate({
        if(input$panel2==2){
          if(exists('MIC')==FALSE) stop('Error!  Was data loaded?')
          if(input$oneBrkpt==FALSE){
            minWidth=as.numeric(input$minWidthM1)
            maxWidth=as.numeric(input$maxWidthM1)
            parms=runLogistic(MIC,DIA,xcens,ycens,MICBrkptL,MICBrkptU,minWidth,maxWidth,session)
            dataLog <<- data.frame(dens=parms$medianDensity,fit=parms$medFit)
            medianDensityLog<<-parms$medianDensity
            medFitLog<<-parms$medFit
            lowerFitLog<<-parms$lowerFit
            upperFitLog<<-parms$upperFit
            modelBrkptLog<<-parms$a1
            D1Log <<- parms$D1
            D2Log <<- parms$D2
            ### Model Panel 4
            output$D1Set1 <- renderUI({
              selectInput("D1Set1", "Set 1 Lower DIA Brkpt:",choices =seq(6,50,by=1),selected=parms$D1)
            })
            output$D2Set1 <- renderUI({
              selectInput('D2Set1', 'Set 1 Upper DIA Brkpt:',choices =seq(6,50,by=1),selected=parms$D2)
            })
            output$D1Set2 <- renderUI({
              selectInput("D1Set2", "Set 2 Lower DIA Brkpt:",choices =seq(6,50,by=1),selected=parms$D1)
            })
            output$D2Set2 <- renderUI({
              selectInput('D2Set2', 'Set 2 Upper DIA Brkpt:',choices =seq(6,50,by=1),selected=parms$D2)
            })
          }else{
            parms=runLogisticOne(MIC,DIA,xcens,ycens,MICBrkpt,session)
            dataLog <<- data.frame(dens=parms$medianDensity,fit=parms$medFit)
            medianDensityLog<<-parms$medianDensity
            medFitLog<<-parms$medFit
            lowerFitLog<<-parms$lowerFit
            upperFitLog<<-parms$upperFit
            DIABrkpt <<- parms$DIABrkpt
            output$D1One <- renderUI({
              selectInput("D1One", "DIABrkpt 1:",choices =seq(6,50,by=1),selected=parms$DIABrkpt)
            })
            output$D2One <- renderUI({
              selectInput('D2One', 'DIABrkpt 2:',choices =seq(6,50,by=1),selected=parms$DIABrkpt+1)
            })
          }
          invisible()
        }
      }))
    }else{return(invisible())}
  })
  
  output$logisticPlot <- renderPlot({
    if(input$actionLogPlot!=0){
      return(isolate({
        if(input$panel2==2){
          convert=FALSE
          if(input$miclogM1=='1') convert=TRUE
          if(input$oneBrkpt==FALSE){
            fit=plotLogistic(MIC,DIA,xcens,ycens,MICBrkptL,MICBrkptU,D1Log,D2Log,
               medianDensityLog,medFitLog,lowerFitLog,upperFitLog,
               flipGraph=input$FlipMod1,logConvert=convert)

          }else{
            fit=plotLogisticOne(MIC,DIA,xcens,ycens,MICBrkpt,DIABrkpt-.5,
                 medianDensityLog,medFitLog,lowerFitLog,upperFitLog,
                 flipGraph=input$FlipMod1,logConvert=convert)
          }
          plot4 <<- fit
          options(warn=-1)
          plot(fit)
          options(warn=1)
        }
      }))  
    }else{return(invisible())}
  })
  
  output$logisticBrkptPlot <- renderPlot({
    if(input$actionLogPlot!=0 & input$oneBrkpt==FALSE){
      return(isolate({
        if(input$panel2==2){
          plotBrkptModel(modelBrkptLog)
          invisible()
        }
      }))
    }else{return(NULL)}
  })
  
  output$downloadLog <- downloadHandler(
    filename=function(){paste('Logistic', Sys.Date(), '.png', sep='')},
    content=function(file){
      png(file,width=1000,height=800,units='px')
      print(plot4)
      dev.off()
    }
  )  
  
  ### Model Panel 2
  output$spline <- renderPrint({
    if(input$actionSpline!=0){
      return(isolate({
        if(input$panel2==3){
          if(exists('MIC')==FALSE) stop('Error!  Was data loaded?')
          if(input$oneBrkpt==FALSE){
            minWidth=as.numeric(input$minWidthM2)
            maxWidth=as.numeric(input$maxWidthM2)
            parms=runSpline(MIC,DIA,xcens,ycens,MICBrkptL,MICBrkptU,minWidth,maxWidth,session)
            dataSpline <<- data.frame(dens=parms$medianDensity,fit=parms$medFit)
            medianDensitySpline<<-parms$medianDensity
            medFitSpline<<-parms$medFit
            lowerFitSpline<<-parms$lowerFit
            upperFitSpline<<-parms$upperFit
            modelBrkptSpline<<-parms$a1
            D1Spline <<- parms$D1
            D2Spline <<- parms$D2
            ### Model Panel 4
            output$D1Set1 <- renderUI({
              selectInput("D1Set1", "Set 1 Lower DIA Brkpt:",choices =seq(6,50,by=1),selected=parms$D1)
            })
            output$D2Set1 <- renderUI({
              selectInput('D2Set1', 'Set 1 Upper DIA Brkpt:',choices =seq(6,50,by=1),selected=parms$D2)
            })
            output$D1Set2 <- renderUI({
              selectInput("D1Set2", "Set 2 Lower DIA Brkpt:",choices =seq(6,50,by=1),selected=parms$D1)
            })
            output$D2Set2 <- renderUI({
              selectInput('D2Set2', 'Set 2 Upper DIA Brkpt:',choices =seq(6,50,by=1),selected=parms$D2)
            })
          }else{
            parms=runSplineOne(MIC,DIA,xcens,ycens,MICBrkpt,session)
            dataSpline <<- data.frame(dens=parms$medianDensity,fit=parms$medFit)
            medianDensitySpline<<-parms$medianDensity
            medFitSpline<<-parms$medFit
            lowerFitSpline<<-parms$lowerFit
            upperFitSpline<<-parms$upperFit
            DIABrkpt <<- parms$DIABrkp
            output$D1One <- renderUI({
              selectInput("D1One", "DIABrkpt 1:",choices =seq(6,50,by=1),selected=parms$DIABrkp)
            })
            output$D2One <- renderUI({
              selectInput('D2One', 'DIABrkpt 2:',choices =seq(6,50,by=1),selected=parms$DIABrkp+1)
            })
          }
          invisible()
        }
      }))
    }else{return(invisible())}
  })
  
  output$splinePlot <- renderPlot({
    if(input$actionSplinePlot!=0){
      return(isolate({
        if(input$panel2==3){
          convert=FALSE
          if(input$miclogM2=='1') convert=TRUE
          if(input$oneBrkpt==FALSE){
            fit=plotSpline(MIC,DIA,xcens,ycens,MICBrkptL,MICBrkptU,D1Spline,D2Spline,
               medianDensitySpline,medFitSpline,lowerFitSpline,upperFitSpline,
               flipGraph=input$FlipMod2,logConvert=convert)
          }else{
            fit=plotSplineOne(MIC,DIA,xcens,ycens,MICBrkpt,DIABrkpt-.5,
                     medianDensitySpline,medFitSpline,lowerFitSpline,upperFitSpline,
                     flipGraph=input$FlipMod2,logConvert=convert)
          }
          plot5 <<- fit
          options(warn=-1)
          plot(fit)
          options(warn=1)
        }
      }))
    }else{return(invisible())}
  })
  
  output$splineBrkptPlot <- renderPlot({
    if(input$actionSplinePlot!=0 & input$oneBrkpt==FALSE){
      return(isolate({
        if(input$panel2==3){
          plotBrkptModel(modelBrkptSpline)
          invisible()
        }
      }))
    }else{return(invisible())}
  })
  
  output$downloadSpline <- downloadHandler(
    filename=function(){paste('Spline', Sys.Date(), '.png', sep='')},
    content=function(file){
      png(file,width=1000,height=800,units='px')
      print(plot5)
      dev.off()
    }
  )
  
  ### Model Panel 3
  output$compareFits <- renderPrint({
    if(input$actionCompare!=0){
      return(isolate({
        if(input$panel2==4){
          if(exists('dataLog')==FALSE | exists('dataSpline')==FALSE)
            stop('Please run both the spline and logistic models first.')
          if(exists('dataSpline')==TRUE & exists('dataLog')==TRUE)
            cat('Models found. Comparing fits....\n')  
          return(invisible())
        }
      }))
    }else{return(invisible())}
  })
  
  output$compareFitsPlot <- renderPlot({
    if(input$actionCompare!=0){
      return(isolate({
        if(input$panel2==4){
          if(exists('dataSpline')==TRUE & exists('dataLog')==TRUE){
            convert=FALSE
            if(input$miclogM3=='1') convert=TRUE
            if(input$oneBrkpt==FALSE){
              fit=compareFitsPlot(dataSpline,dataLog,MIC,DIA,xcens,ycens,MICBrkptL,MICBrkptU,
                flipGraph=input$FlipMod3,logConvert=convert)
            }else{
              fit=compareFitsPlotOne(dataSpline,dataLog,MIC,DIA,xcens,ycens,MICBrkpt,flipGraph=input$FlipMod3,logConvert=convert)
            }
            plotC <<- fit
            options(warn=-1)
            plot(fit)
            options(warn=1)
          }
          return(invisible())
        }
      }))
    }else{return(NULL)}
  })
  
  output$downloadCompare <- downloadHandler(
    filename=function(){paste('CompareFits', Sys.Date(), '.tiff', sep='')},
    content=function(file){
      tiff(file,width=900,height=900,units='px')
      print(plotC)
      dev.off()
    }
  )
  
  
  ### Model Panel 4
  calcProb <- reactive({
    paste(capture.output({
      if(input$actionProbDIA!=0){
        return(isolate({
          if(input$panel2==5){
            if(input$oneBrkpt==FALSE){
              DIASet2=c(0,0)
              
              if(input$modelSelect=='spline')
                if(exists('dataSpline')==FALSE)
                  stop('Run spline model first (second tab).\n')
              if(input$modelSelect=='logistic')
                if(exists('dataLog')==FALSE)
                  stop('Run logistic model first (first tab).\n')
              
              if(input$modelSelect=='spline') a1=dataSpline
              if(input$modelSelect=='logistic') a1=dataLog
              
              DIASet1=c(as.numeric(input$D1Set1),as.numeric(input$D2Set1))
              if(input$includeSecondDIA!=0){
                DIASet2=c(as.numeric(input$D1Set2),as.numeric(input$D2Set2))
              }else{DIASet2=c(NA,NA)}
              
              if(DIASet1[2]<DIASet1[1])
                stop('Lower DIA breakpoint Set 1 must be less than upper DIA breakpoint. \n')
              if(is.na(DIASet2[1])==FALSE)
                if(DIASet2[2]<DIASet2[1])
                  stop('Lower DIA breakpoint Set 2 must be less than upper DIA breakpoint. \n')
              probDIAData <<- probDIAClass(a1,MICBrkptL,MICBrkptU,DIASet1,DIASet2)
            }else{
              DIA2=0
              if(input$modelSelect=='spline')
                if(exists('dataSpline')==FALSE)
                  stop('Run spline model first (second tab).\n')
              if(input$modelSelect=='logistic')
                if(exists('dataLog')==FALSE)
                  stop('Run logistic model first (first tab).\n')
              if(input$modelSelect=='spline') a1=dataSpline
              if(input$modelSelect=='logistic') a1=dataLog
              
              DIA1=as.numeric(input$D1One)
              if(input$includeSecondDIAOne!=0){
                DIA2=as.numeric(input$D2One)
              }else{DIA2=NA}
              probDIAData <<-probDIAClassOne(a1,MICBrkpt,DIA1,DIA2)
            }
            return(invisible())
          }
        }))
      }else{invisible()}
    }), collapse="\n")
  })
  
  output$probDIA <- renderText({
   return(calcProb())
  })
  
  output$probDIAPlot <- renderPlot({
    calcProb()
    if(input$actionProbDIA!=0){
      return(isolate({
        if(input$panel2==5){
          
          convert=FALSE
          if(input$miclogM4=='1') convert=TRUE
          
          if(input$oneBrkpt==FALSE){
            
            DIASet1=c(input$D1Set1,input$D2Set1)
            if(input$includeSecondDIA!=0){
              DIASet2=c(input$D1Set2,input$D2Set2)
            }else{DIASet2=NA}
          
            plotProb <<- plotProbDIAClass(probDIAData,MICBrkptL,MICBrkptU,DIASet1,DIASet2,logConvert=convert)
            
          }else{
            
            DIA1=as.numeric(input$D1One)
            if(input$includeSecondDIAOne!=0){
              DIA2=as.numeric(input$D2One)
            }else{DIA2=NA}
            
            plotProb <<- plotProbDIAClassOne(probDIAData,MICBrkpt,DIA1,DIA2,logConvert=convert)
            
          }
          
          options(warn=-1)
          plot(plotProb)
          options(warn=1)
          return(invisible())
        }
      }))
    }else{return(invisible())}
  })
  
  output$downloadProbDIA <- downloadHandler(
    filename=function(){paste('Probability Classification', Sys.Date(), '.tiff', sep='')},
    content=function(file){
      tiff(file,width=900,height=900,units='px')
      print(plotProb)
      dev.off()
    }
  )
  
  output$plotDist <- renderPlot({
    if(input$actionProbDIA!=0){
      return(isolate({
        if(input$panel2==5){
          
          convert=FALSE
          if(input$miclogM4=='1') convert=TRUE
          
          if(input$oneBrkpt==FALSE){
            plotDist <<- plotUnderlyingDistibution(probDIAData,convert,MICBrkptL,MICBrkptU)
          }else{
            plotDist <<- plotUnderlyingDistibutionOne(probDIAData,convert,MICBrkpt)
          }
          
          options(warn=-1)
          plot(plotDist)
          options(warn=1)
          return(invisible())
        }
      }))
    }else{return(invisible())}
  })
  
  
}


