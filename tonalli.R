library(shiny)
library(shinyWidgets)
library(tidyr)
library(dplyr)
library(ggplot2)
tonalli <- read.csv('data/tonalli.csv')

# FALTA cosas de ciencia

# Define UI ----
ui <- fluidPage(
    fluidRow(
     column(4, div( style = "margin-top: 15px; color: #006400; font-size: 40px;", tags$b("BrachyTonalli 1.0"))),
     column(4, div( style = "margin-top: 0px;", img(src = "brachy.jpg", height = 140, width = 330))),
     column(4, tags$blockquote("BrachyTonalli is named after the Mexica concept of 'Tonalli', meaning 'of the Sun', 
                               that explains the vital force living beings use to associate and palpitate with the Sun rhythms")
            )),
    fluidRow(
     column(4, ),
     column(4, h5("BrachyTonalli is a companion ShinyApp to explore the transcriptome of Brachypodium distachyon (leaves) under submergence stress, 
                  low light stress and normal growth conditions reported by: "), 
            tags$a(href="https://doi.org/10.3390/ijms24108555", "Medina-Chávez et al (2022)"), style = "background-color: #FDF5E6", ),
            ),
    fluidRow(
      column(4, style = "background-color: #FFFFF0", textInput(inputId = "caption", label = "Enter your Brachypodium gene: e.g. Bradi3g16515", value = "Bradi3g16515"), submitButton("Run")),
    ),
              fluidRow(
                column(12, style = "background-color: #FFF8DC;", tableOutput(outputId = "TuGen"), downloadButton('downloadFullYourGene', 'Save .CSV'))),
              fluidRow(
                column(6, style = "background-color: #F4A460;", plotOutput(outputId = "sensitive")), 
                column(6, style = "background-color: #F4A460;", plotOutput(outputId = "tolerant"))
                ),
                  fluidRow(
                    column(6, style = "background-color: #F4A460;", downloadButton('downloadBd21cpm', 'Save .PNG')), 
                    column(6, style = "background-color: #F4A460;", downloadButton('downloadBd213cpm', 'Save .PNG'))
                ),
              fluidRow(
                column(4, style = "background-color: white", plotOutput(outputId = "LOGS"), tableOutput(outputId = "LOGS_table")),
                column(4, style = "background-color: white", plotOutput(outputId = "LOGLL"), tableOutput(outputId = "LOGLL_table")),
                column(4, style = "background-color: white", plotOutput(outputId = "LOGNG"), tableOutput(outputId = "LOGNG_table"))
              ),
                  fluidRow(
                    column(4, style = "background-color: #D2691E;", downloadButton('downloadLOG_SvsLL', 'Save .PNG')),
                    column(4, style = "background-color: #D2691E;", downloadButton('downloadLOG_LLvsNG', 'Save .PNG')),
                    column(4, style = "background-color: #D2691E;", downloadButton('downloadLOG_NGvs8NG', 'Save .PNG'))
              ),
    
              fluidRow(
                column(4, style = "background-color: white;", plotOutput(outputId = "LOGSvs"), tableOutput(outputId = "LOGSvs_table")),
                column(4, style = "background-color: white;", plotOutput(outputId = "LOGLLvs"), tableOutput(outputId = "LOGLLvs_table")),
                column(4, style = "background-color: white;", plotOutput(outputId = "LOGNGvs"), tableOutput(outputId = "LOGNGvs_table"))
              ),
                  fluidRow(
                    column(4, style = "background-color: #FFFFE0;", downloadButton('downloadLOGSvs', 'Save .PNG')),
                    column(4, style = "background-color: #FFFFE0;", downloadButton('downloadLOGLLvs', 'Save .PNG')),
                    column(4, style = "background-color: #FFFFE0;", downloadButton('downloadLOGNGvs', 'Save .PNG'))
                  ),
              fluidRow(
               column(4, style = "background-color: #FFFFE0;", h5("
                                                                  In our original report, we considered a Log2FC value of 1.5/-1.5 as threshold for
                                                                  up/downregulation, counts per million (CPM) of at least 15 for expression, and a 
                                                                  FDR of maximum 0.000005 as significant. However, these are stringent parameters that 
                                                                  can be changed. For example, for the discovery of low light responsive genes (a less shocking stress), 
                                                                  we have observed that they are too strict.
                                                                  ")),
               column(4, style = "background-color: #FFFFE0;", h5("
                                                                  The raw dataset can be downloaded at NCBI-GEO as project GSE215159. In the following link you
                                                                  can download the analyzed dataset as CSV:
                                                                  "), downloadButton('downloadTonalliCSV', 'Save tonalli.CSV')
                                                                  ),
               column(4, style = "background-color: #FFFFE0;", h5("
                                                                  These experiments, sequencing and datasets were funded by grant
                                                                  287137 (CONHCYT-México). Many thanks to R developers for making all these
                                                                  tools avaialble.
                                                                  "))
            )
)

# Define server logic ---- 
server <- function(input, output, session) {

       # Preparing Tonalli for download
  
       tonalliCSVdownload <- reactive({
       tonalliCSV <- tonalli
       })
  
       output$tonalliCSV2 <- renderTable({tonalliCSVdownload()})
  
       output$downloadTonalliCSV <- downloadHandler(
       filename = "tonalli_full.csv",
       content = function(file){
       write.csv(tonalliCSVdownload(), file)
       })
  
  
      # Displaying gene identity
      
      output$TuGen <- renderTable({
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        TuGen <- FoundGene[c(1),1:4]
      })
      
      # Downloading full data for your gene
      
      FullYourGeneTable <- reactive({
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        FullYourGene <- FoundGene[c(1),1:239]
      })
      
      output$FullYourGene <- renderTable({FullYourGeneTable()})
      
      output$downloadFullYourGene <- downloadHandler(
        filename = "your_Gene.csv",
        content = function(file){
          write.csv(FullYourGeneTable(), file)
        })
      
      # Graphing, displaying and downloading Bd21 cpm expression
      
      Bd21cpm <- reactive({UserGene <-  input$caption
                                FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
                                ExpSub21 <- FoundGene[c(1),200:204]
                                ExpLL21 <- FoundGene[c(1),190:194]
                                ExpNG21 <- FoundGene[c(1),180:184]
                                SDSub21 <- FoundGene[c(1),205:209]
                                SDLL21 <- FoundGene[c(1),195:199]
                                SDNG21 <- FoundGene[c(1),185:189]
                                
                                ExpSub213 <- FoundGene[c(1),230:234]
                                ExpLL213 <- FoundGene[c(1),220:224]
                                ExpNG213 <- FoundGene[c(1),210:214]
                                SDSub213 <- FoundGene[c(1),235:239]
                                SDLL213 <- FoundGene[c(1),225:229]
                                SDNG213 <- FoundGene[c(1),215:219]
                                
                                ExpSub21Col <- as.data.frame(t(ExpSub21))
                                ExpLL21Col <- as.data.frame(t(ExpLL21))
                                ExpNG21Col <- as.data.frame(t(ExpNG21))
                                SDSub21Col <- as.data.frame(t(SDSub21))
                                SDLL21Col <- as.data.frame(t(SDLL21))
                                SDNG21Col <- as.data.frame(t(SDNG21))
                                
                                ExpSub213Col <- as.data.frame(t(ExpSub213))
                                ExpLL213Col <- as.data.frame(t(ExpLL213))
                                ExpNG213Col <- as.data.frame(t(ExpNG213))
                                SDSub213Col <- as.data.frame(t(SDSub213))
                                SDLL213Col <- as.data.frame(t(SDLL213))
                                SDNG213Col <- as.data.frame(t(SDNG213))
                                
                                SUB21 <- bind_cols(ExpSub21Col, SDSub21Col)
                                LL21 <- bind_cols(ExpLL21Col, SDLL21Col)
                                NG21 <- bind_cols(ExpNG21Col, SDNG21Col)
                                
                                SUB213 <- bind_cols(ExpSub213Col, SDSub213Col)
                                LL213 <- bind_cols(ExpLL213Col, SDLL213Col)
                                NG213 <- bind_cols(ExpNG213Col, SDNG213Col)
                                
                                colnames(SUB21)[1] <- "Expression"
                                colnames(SUB21)[2] <- "DS"
                                colnames(LL21)[1] <- "Expression"
                                colnames(LL21)[2] <- "DS"
                                colnames(NG21)[1] <- "Expression"
                                colnames(NG21)[2] <- "DS"
                                colnames(SUB213)[1] <- "Expression"
                                colnames(SUB213)[2] <- "DS"
                                colnames(LL213)[1] <- "Expression"
                                colnames(LL213)[2] <- "DS"
                                colnames(NG213)[1] <- "Expression"
                                colnames(NG213)[2] <- "DS"
                                
                                SUB21data <- data.frame(z=SUB21, Condition=c("Submergence"))
                                LL21data <- data.frame(z=LL21, Condition=c("Low Light"))
                                NG21data <- data.frame(z=NG21, Condition=c("Normal Growth"))
                                SUB213data <- data.frame(z=SUB213, Condition=c("Submergence"))
                                LL213data <- data.frame(z=LL213, Condition=c("Low Light"))
                                NG213data <- data.frame(z=NG213, Condition=c("Normal Growth"))
                                
                                colnames(SUB21data)[1] <- "Expression"
                                colnames(SUB21data)[2] <- "DS"
                                colnames(LL21data)[1] <- "Expression"
                                colnames(LL21data)[2] <- "DS"
                                colnames(NG21data)[1] <- "Expression"
                                colnames(NG21data)[2] <- "DS"
                                colnames(SUB213data)[1] <- "Expression"
                                colnames(SUB213data)[2] <- "DS"
                                colnames(LL213data)[1] <- "Expression"
                                colnames(LL213data)[2] <- "DS"
                                colnames(NG213data)[1] <- "Expression"
                                colnames(NG213data)[2] <- "DS"
                                SUB21dataxgraph<-bind_rows(SUB21data, LL21data, NG21data)
                                SUB21dataxgraph2 <- data.frame(z=SUB21dataxgraph, Time=c("0","8","16","20","24"))
                                SUB213dataxgraph<-bind_rows(SUB213data, LL213data, NG213data)
                                SUB213dataxgraph2 <- data.frame(z=SUB213dataxgraph, Time=c("0","8","16","20","24"))
                                
                                colnames(SUB21dataxgraph2)[1] <- "Expression"
                                colnames(SUB21dataxgraph2)[2] <- "DS"
                                colnames(SUB21dataxgraph2)[3] <- "Condition"
                                colnames(SUB213dataxgraph2)[1] <- "Expression"
                                colnames(SUB213dataxgraph2)[2] <- "DS"
                                colnames(SUB213dataxgraph2)[3] <- "Condition"
                                
                                SUB21dataxgraph2$Time <- as.numeric(as.character(SUB21dataxgraph2$Time))
                                SUB213dataxgraph2$Time <- as.numeric(as.character(SUB213dataxgraph2$Time))
                                
                                YourGene <- UserGene
                                
                                ggplot(SUB21dataxgraph2, aes(x=Time, y=Expression))+
                                  geom_point(aes(color=Condition, shape=Condition),size=4)+
                                  geom_line(aes(group=Condition, color=Condition), linewidth=1.5)+
                                  geom_errorbar(aes(ymin=Expression-DS,ymax=Expression+DS), size=0.4, width=1)+
                                  theme(aspect.ratio = 1)+
                                  ggtitle(YourGene, "Brachypodium distachyon Bd21 CPM expression")
                                })
      
      output$sensitive <- renderPlot({Bd21cpm()})
      
      output$downloadBd21cpm <- downloadHandler(
        filename = "Bd21cpm.png",
        content = function(file){
          ggsave(file, plot = Bd21cpm(), width = 7, height = 5, device = "png")
        })
      
      
      # Graphing, displaying and downloading Bd21-3 cpm expression
      
      Bd213cpm <- reactive({
        
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        ExpSub21 <- FoundGene[c(1),200:204]
        ExpLL21 <- FoundGene[c(1),190:194]
        ExpNG21 <- FoundGene[c(1),180:184]
        SDSub21 <- FoundGene[c(1),205:209]
        SDLL21 <- FoundGene[c(1),195:199]
        SDNG21 <- FoundGene[c(1),185:189]
        
        ExpSub213 <- FoundGene[c(1),230:234]
        ExpLL213 <- FoundGene[c(1),220:224]
        ExpNG213 <- FoundGene[c(1),210:214]
        SDSub213 <- FoundGene[c(1),235:239]
        SDLL213 <- FoundGene[c(1),225:229]
        SDNG213 <- FoundGene[c(1),215:219]
        
        ExpSub21Col <- as.data.frame(t(ExpSub21))
        ExpLL21Col <- as.data.frame(t(ExpLL21))
        ExpNG21Col <- as.data.frame(t(ExpNG21))
        SDSub21Col <- as.data.frame(t(SDSub21))
        SDLL21Col <- as.data.frame(t(SDLL21))
        SDNG21Col <- as.data.frame(t(SDNG21))
        
        ExpSub213Col <- as.data.frame(t(ExpSub213))
        ExpLL213Col <- as.data.frame(t(ExpLL213))
        ExpNG213Col <- as.data.frame(t(ExpNG213))
        SDSub213Col <- as.data.frame(t(SDSub213))
        SDLL213Col <- as.data.frame(t(SDLL213))
        SDNG213Col <- as.data.frame(t(SDNG213))
        
        SUB21 <- bind_cols(ExpSub21Col, SDSub21Col)
        LL21 <- bind_cols(ExpLL21Col, SDLL21Col)
        NG21 <- bind_cols(ExpNG21Col, SDNG21Col)
        
        SUB213 <- bind_cols(ExpSub213Col, SDSub213Col)
        LL213 <- bind_cols(ExpLL213Col, SDLL213Col)
        NG213 <- bind_cols(ExpNG213Col, SDNG213Col)
        
        colnames(SUB21)[1] <- "Expression"
        colnames(SUB21)[2] <- "DS"
        colnames(LL21)[1] <- "Expression"
        colnames(LL21)[2] <- "DS"
        colnames(NG21)[1] <- "Expression"
        colnames(NG21)[2] <- "DS"
        colnames(SUB213)[1] <- "Expression"
        colnames(SUB213)[2] <- "DS"
        colnames(LL213)[1] <- "Expression"
        colnames(LL213)[2] <- "DS"
        colnames(NG213)[1] <- "Expression"
        colnames(NG213)[2] <- "DS"
        
        SUB21data <- data.frame(z=SUB21, Condition=c("Submergence"))
        LL21data <- data.frame(z=LL21, Condition=c("Low Light"))
        NG21data <- data.frame(z=NG21, Condition=c("Normal Growth"))
        SUB213data <- data.frame(z=SUB213, Condition=c("Submergence"))
        LL213data <- data.frame(z=LL213, Condition=c("Low Light"))
        NG213data <- data.frame(z=NG213, Condition=c("Normal Growth"))
        
        colnames(SUB21data)[1] <- "Expression"
        colnames(SUB21data)[2] <- "DS"
        colnames(LL21data)[1] <- "Expression"
        colnames(LL21data)[2] <- "DS"
        colnames(NG21data)[1] <- "Expression"
        colnames(NG21data)[2] <- "DS"
        colnames(SUB213data)[1] <- "Expression"
        colnames(SUB213data)[2] <- "DS"
        colnames(LL213data)[1] <- "Expression"
        colnames(LL213data)[2] <- "DS"
        colnames(NG213data)[1] <- "Expression"
        colnames(NG213data)[2] <- "DS"
        SUB21dataxgraph<-bind_rows(SUB21data, LL21data, NG21data)
        SUB21dataxgraph2 <- data.frame(z=SUB21dataxgraph, Time=c("0","8","16","20","24"))
        SUB213dataxgraph<-bind_rows(SUB213data, LL213data, NG213data)
        SUB213dataxgraph2 <- data.frame(z=SUB213dataxgraph, Time=c("0","8","16","20","24"))
        
        colnames(SUB21dataxgraph2)[1] <- "Expression"
        colnames(SUB21dataxgraph2)[2] <- "DS"
        colnames(SUB21dataxgraph2)[3] <- "Condition"
        colnames(SUB213dataxgraph2)[1] <- "Expression"
        colnames(SUB213dataxgraph2)[2] <- "DS"
        colnames(SUB213dataxgraph2)[3] <- "Condition"
        
        SUB21dataxgraph2$Time <- as.numeric(as.character(SUB21dataxgraph2$Time))
        SUB213dataxgraph2$Time <- as.numeric(as.character(SUB213dataxgraph2$Time))
        
        YourGene <- UserGene
        
        ggplot(SUB213dataxgraph2, aes(x=Time, y=Expression))+
          geom_point(aes(color=Condition, shape=Condition),size=4)+
          geom_line(aes(group=Condition, color=Condition), linewidth=1.5)+
          geom_errorbar(aes(ymin=Expression-DS,ymax=Expression+DS), size=0.4, width=1)+
          theme(aspect.ratio = 1)+
          ggtitle(YourGene, "Brachypodium distachyon Bd21-3 CPM expression") 
      })
      
      output$tolerant <- renderPlot({Bd213cpm()})
      
      output$downloadBd213cpm <- downloadHandler(
        filename = "Bd21-3cpm.png",
        content = function(file){
          ggsave(file, plot = Bd213cpm(), width = 7, height = 5, device = "png")
        })
      
      # Graphing, displaying and downloading Bd21 and Bd21-3 individual Submergence/Low Light Log2FC comparison
  
      LOG_SvsLL <- reactive({
        
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        LogSub21 <- FoundGene[c(1),35:39]
        LogSub213 <- FoundGene[c(1),45:49]
        
        
        LogSub21Col <- as.data.frame(t(LogSub21))
        LogSub213Col <- as.data.frame(t(LogSub213))
        
        LogSub21data <- data.frame(z=LogSub21Col, Ecotype=c("Bd21"))
        LogSub213data <- data.frame(z=LogSub213Col, Ecotype=c("Bd21-3"))
        
        colnames(LogSub21data)[1] <- "Log"
        colnames(LogSub213data)[1] <- "Log"
        
        LogSUBdataxgraph<-bind_rows(LogSub21data, LogSub213data)
        LogSUBdataxgraph2 <- data.frame(z=LogSUBdataxgraph, Time=c("0","8","16","20","24"))
        
        
        colnames(LogSUBdataxgraph2)[1] <- "Log2FC"
        colnames(LogSUBdataxgraph2)[2] <- "Ecotype"
        
        LogSUBdataxgraph2$Time <- as.numeric(as.character(LogSUBdataxgraph2$Time))
        
        YourGene <- UserGene
        
        ggplot(LogSUBdataxgraph2, aes(x=Time, y=Log2FC))+
          geom_point(aes(color=Ecotype, shape=Ecotype),size=4)+
          geom_line(aes(group=Ecotype, color=Ecotype), linewidth=1.5)+
          theme(aspect.ratio = 1)+
          ggtitle(YourGene, "Brachypodium distachyon Log2FC Submergence/Low Light") 
      })
      
      output$LOGS <- renderPlot({LOG_SvsLL()})
      
      output$downloadLOG_SvsLL <- downloadHandler(
        filename = "LOG_SvsLL.png",
        content = function(file){
          ggsave(file, plot = LOG_SvsLL(), width = 7, height = 5, device = "png")
        })
     
      LOG_SvsLL_table <- reactive({
        
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        LogSub21fdr <- FoundGene[c(1),40:44]
        LogSub213fdr <- FoundGene[c(1),50:54]
        
        
        LogSub21Colfdr <- as.data.frame(t(LogSub21fdr))
        LogSub213Colfdr <- as.data.frame(t(LogSub213fdr))
        
        LogSub21datafdr <- data.frame(z=LogSub21Colfdr, Ecotype=c("Bd21"))
        LogSub213datafdr <- data.frame(z=LogSub213Colfdr, Ecotype=c("Bd21-3"))
        
        colnames(LogSub21datafdr)[1] <- "FDR"
        colnames(LogSub213datafdr)[1] <- "FDR"
        
        LogSUBdataxgraphfdr<-bind_rows(LogSub21datafdr, LogSub213datafdr)
        LogSUBdataxgraph2fdr <- data.frame(z=LogSUBdataxgraphfdr, Time=c("0","8","16","20","24"))
        
        
        colnames(LogSUBdataxgraph2fdr)[1] <- "FDR"
        colnames(LogSUBdataxgraph2fdr)[2] <- "Ecotype"
        
        LogSUBdataxgraph2fdr$Time <- as.numeric(as.character(LogSUBdataxgraph2fdr$Time))
        
        LOG_SvsLL_table <- LogSUBdataxgraph2fdr
      })
      
      output$LOGS_table <- renderTable({LOG_SvsLL_table()}, digits = 6)
      
       # Graphing, displaying and downloading Bd21 and Bd21-3 individual Low Light/Normal Growth Log2FC comparison
      
      LOG_LLvsNG <- reactive({
        
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        LogLL21 <- FoundGene[c(1),5:9]
        LogLL213 <- FoundGene[c(1),15:19]
        
        
        LogLL21Col <- as.data.frame(t(LogLL21))
        LogLL213Col <- as.data.frame(t(LogLL213))
        
        LogLL21data <- data.frame(z=LogLL21Col, Ecotype=c("Bd21"))
        LogLL213data <- data.frame(z=LogLL213Col, Ecotype=c("Bd21-3"))
        
        colnames(LogLL21data)[1] <- "Log"
        colnames(LogLL213data)[1] <- "Log"
        
        LogLLdataxgraph<-bind_rows(LogLL21data, LogLL213data)
        LogLLdataxgraph2 <- data.frame(z=LogLLdataxgraph, Time=c("0","8","16","20","24"))
        
        
        colnames(LogLLdataxgraph2)[1] <- "Log2FC"
        colnames(LogLLdataxgraph2)[2] <- "Ecotype"
        
        LogLLdataxgraph2$Time <- as.numeric(as.character(LogLLdataxgraph2$Time))
        
        YourGene <- UserGene
        
        ggplot(LogLLdataxgraph2, aes(x=Time, y=Log2FC))+
          geom_point(aes(color=Ecotype, shape=Ecotype),size=4)+
          geom_line(aes(group=Ecotype, color=Ecotype), linewidth=1.5)+
          theme(aspect.ratio = 1)+
          ggtitle(YourGene, "Brachypodium distachyon Log2FC Low Light / Normal Growth") 
      })
      
      output$LOGLL <- renderPlot({LOG_LLvsNG()})
      
      output$downloadLOG_LLvsNG <- downloadHandler(
        filename = "LOG_LLvsNG.png",
        content = function(file){
          ggsave(file, plot = LOG_LLvsNG(), width = 7, height = 5, device = "png")
        })
      
      LOG_LLvsNG_table <- reactive({
        
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        LogLL21fdr <- FoundGene[c(1),10:14]
        LogLL213fdr <- FoundGene[c(1),20:24]
        
        LogLL21Colfdr <- as.data.frame(t(LogLL21fdr))
        LogLL213Colfdr <- as.data.frame(t(LogLL213fdr))
        
        LogLL21datafdr <- data.frame(z=LogLL21Colfdr, Ecotype=c("Bd21"))
        LogLL213datafdr <- data.frame(z=LogLL213Colfdr, Ecotype=c("Bd21-3"))
        
        colnames(LogLL21datafdr)[1] <- "FDR"
        colnames(LogLL213datafdr)[1] <- "FDR"
        
        LogLLdataxgraphfdr <- bind_rows(LogLL21datafdr, LogLL213datafdr)
        LogLLdataxgraph2fdr <- data.frame(z=LogLLdataxgraphfdr, Time=c("0","8","16","20","24"))
        
        
        colnames(LogLLdataxgraph2fdr)[1] <- "FDR"
        colnames(LogLLdataxgraph2fdr)[2] <- "Ecotype"
        
        LogLLdataxgraph2fdr$Time <- as.numeric(as.character(LogLLdataxgraph2fdr$Time))
        
        LOG_LLvsNG_table <- LogLLdataxgraph2fdr
        

      })
      
      output$LOGLL_table <- renderTable({LOG_LLvsNG_table()}, digits = 6)
      
      # Graphing, displaying and downloading Bd21 and Bd21-3 individual X Normal Growth / 8h Normal Growth Log2FC comparison
      
      LOG_NGvs8NG <- reactive({
        
        
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        LogNG21 <- FoundGene[c(1),65:68]
        LogNG213 <- FoundGene[c(1),73:76]
        
        
        LogNG21Col <- as.data.frame(t(LogNG21))
        LogNG213Col <- as.data.frame(t(LogNG213))
        
        LogNG21data <- data.frame(z=LogNG21Col, Ecotype=c("Bd21"))
        LogNG213data <- data.frame(z=LogNG213Col, Ecotype=c("Bd21-3"))
        
        colnames(LogNG21data)[1] <- "Log"
        colnames(LogNG213data)[1] <- "Log"
        
        LogNGdataxgraph<-bind_rows(LogNG21data, LogNG213data)
        LogNGdataxgraph2 <- data.frame(z=LogNGdataxgraph, Time=c("0","16","20","24"))
        
        
        colnames(LogNGdataxgraph2)[1] <- "Log2FC"
        colnames(LogNGdataxgraph2)[2] <- "Ecotype"
        
        LogNGdataxgraph2$Time <- as.numeric(as.character(LogNGdataxgraph2$Time))
        
        YourGene <- UserGene
        
        ggplot(LogNGdataxgraph2, aes(x=Time, y=Log2FC))+
          geom_point(aes(color=Ecotype, shape=Ecotype),size=4)+
          geom_line(aes(group=Ecotype, color=Ecotype), linewidth=1.5)+
          theme(aspect.ratio = 1)+
          ggtitle(YourGene, "Log2FC=(Xh Normal Growth / 8h Normal Growth)") 
      })
      
      output$LOGNG <- renderPlot({LOG_NGvs8NG()})
      
      output$downloadLOG_NGvs8NG <- downloadHandler(
        filename = "LOG_NGvs8NG.png",
        content = function(file){
          ggsave(file, plot = LOG_NGvs8NG(), width = 7, height = 5, device = "png")
        })
      
      LOG_NGvs8NG_table <- reactive({UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        LogNG21fdr <- FoundGene[c(1),69:72]
        LogNG213fdr <- FoundGene[c(1),77:80]
      
      
        LogNG21Colfdr <- as.data.frame(t(LogNG21fdr))
        LogNG213Colfdr <- as.data.frame(t(LogNG213fdr))
      
        LogNG21datafdr <- data.frame(z=LogNG21Colfdr, Ecotype=c("Bd21"))
        LogNG213datafdr <- data.frame(z=LogNG213Colfdr, Ecotype=c("Bd21-3"))
      
        colnames(LogNG21datafdr)[1] <- "FDR"
        colnames(LogNG213datafdr)[1] <- "FDR"
      
        LogNGdataxgraphfdr <- bind_rows(LogNG21datafdr, LogNG213datafdr)
        LogNGdataxgraph2fdr <- data.frame(z=LogNGdataxgraphfdr, Time=c("0","16","20","24"))
      
      
        colnames(LogNGdataxgraph2fdr)[1] <- "FDR"
        colnames(LogNGdataxgraph2fdr)[2] <- "Ecotype"
      
        LogNGdataxgraph2fdr$Time <- as.numeric(as.character(LogNGdataxgraph2fdr$Time))
        
        LOG_NGvs8NG_table <- LogNGdataxgraph2fdr
        })
      
      output$LOGNG_table <- renderTable({LOG_NGvs8NG_table()}, digits = 6)
      
      
      # Graphing, displaying and downloading Bd21 vs Bd21-3 Log2FC comparison under Submergence
      
      SubVS <- reactive({
        
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        Log21vs213Sub <- FoundGene[c(1),55:59]
        Log21vs213SubCol <- as.data.frame(t(Log21vs213Sub))
        Log21vs213Subxgraph2 <- data.frame(z=Log21vs213SubCol, Time=c("0","8","16","20","24"))
        Log21vs213Subxgraph2$Time <- as.numeric(as.character(Log21vs213Subxgraph2$Time))
        colnames(Log21vs213Subxgraph2)[1] <- "Log2FC"
        
        YourGene <- UserGene
        
        ggplot(Log21vs213Subxgraph2, aes(x=Time, y=Log2FC))+
          geom_point(aes(), color="#619cff", size=4)+
          geom_line(aes(group=1), color ="#619cff", linewidth=1.5)+
          theme(aspect.ratio = 1)+
          ggtitle(YourGene, "Brachypodium distachyon Bd21/Bd21-3 Log2FC under submergence")
      })
      
      output$LOGSvs <- renderPlot({SubVS()})
      
      output$downloadLOGSvs <- downloadHandler(
        filename = "Bd21vsBd21-3_logSub.png",
        content = function(file){
          ggsave(file, plot = SubVS(), width = 7, height = 5, device = "png")
        })
      
      SubVStable <- reactive({
        
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        Log21vs213Subfdr <- FoundGene[c(1),60:64]
        Log21vs213SubColfdr <- as.data.frame(t(Log21vs213Subfdr))
        Log21vs213Subxgraph2fdr <- data.frame(z=Log21vs213SubColfdr, Time=c("0","8","16","20","24"))
        Log21vs213Subxgraph2fdr$Time <- as.numeric(as.character(Log21vs213Subxgraph2fdr$Time))
        colnames(Log21vs213Subxgraph2fdr)[1] <- "FDR"
        SubVStable <- Log21vs213Subxgraph2fdr
      })
      
      output$LOGSvs_table <- renderTable({SubVStable()}, digits = 6)
      
      # Graphing, displaying and downloading Bd21 vs Bd21-3 Log2FC comparison under Low Light
      
      LLVS <- reactive({
        
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        Log21vs213LL <- FoundGene[c(1),25:29]
        Log21vs213LLCol <- as.data.frame(t(Log21vs213LL))
        Log21vs213LLxgraph2 <- data.frame(z=Log21vs213LLCol, Time=c("0","8","16","20","24"))
        Log21vs213LLxgraph2$Time <- as.numeric(as.character(Log21vs213LLxgraph2$Time))
        colnames(Log21vs213LLxgraph2)[1] <- "Log2FC"
        
        YourGene <- UserGene
        
        ggplot(Log21vs213LLxgraph2, aes(x=Time, y=Log2FC))+
          geom_point(aes(), color="#f8766d", size=4)+
          geom_line(aes(group=1), color ="#f8766d", linewidth=1.5)+
          theme(aspect.ratio = 1)+
          ggtitle(YourGene, "Brachypodium distachyon Bd21/Bd21-3 Log2FC under Low light")
      })
      
      output$LOGLLvs <- renderPlot({LLVS()})
      
      output$downloadLOGLLvs <- downloadHandler(
        filename = "Bd21vsBd21-3_logLL.png",
        content = function(file){
          ggsave(file, plot = LLVS(), width = 7, height = 5, device = "png")
        })
      
      LLVStable <- reactive({
        
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        Log21vs213LLfdr <- FoundGene[c(1),30:34]
        Log21vs213LLColfdr <- as.data.frame(t(Log21vs213LLfdr))
        Log21vs213LLxgraph2fdr <- data.frame(z=Log21vs213LLColfdr, Time=c("0","8","16","20","24"))
        Log21vs213LLxgraph2fdr$Time <- as.numeric(as.character(Log21vs213LLxgraph2fdr$Time))
        colnames(Log21vs213LLxgraph2fdr)[1] <- "FDR"
        LLVStable <- Log21vs213LLxgraph2fdr
        
      })
      
      output$LOGLLvs_table <- renderTable({LLVStable()}, digits = 6)
      
      # Graphing, displaying and downloading Bd21 and Bd21-3 versus Log2FC comparison under Normal Growth
      
      NGVS <- reactive({
        
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        Log21vs213NG <- FoundGene[c(1),81:85]
        Log21vs213NGCol <- as.data.frame(t(Log21vs213NG))
        Log21vs213NGxgraph2 <- data.frame(z=Log21vs213NGCol, Time=c("0","8","16","20","24"))
        Log21vs213NGxgraph2$Time <- as.numeric(as.character(Log21vs213NGxgraph2$Time))
        colnames(Log21vs213NGxgraph2)[1] <- "Log2FC"
        
        YourGene <- UserGene
        
        ggplot(Log21vs213NGxgraph2, aes(x=Time, y=Log2FC))+
          geom_point(aes(), color="#00ba2e", size=4)+
          geom_line(aes(group=1), color ="#00ba2e", linewidth=1.5)+
          theme(aspect.ratio = 1)+
          ggtitle(YourGene, "Brachypodium distachyon Bd21/Bd21-3 Log2FC under Normal Growth")
      })
      
      output$LOGNGvs <- renderPlot({NGVS()})
      
      output$downloadLOGNGvs <- downloadHandler(
        filename = "Bd21vsBd21-3_logNG.png",
        content = function(file){
          ggsave(file, plot = NGVS(), width = 7, height = 5, device = "png")
        })
      
      NGVS_table <- reactive({
        
        UserGene <- input$caption
        FoundGene <- tonalli[which(tonalli$Gene==UserGene, arr.ind=TRUE)[1],]
        Log21vs213NGfdr <- FoundGene[c(1),86:90]
        Log21vs213NGColfdr <- as.data.frame(t(Log21vs213NGfdr))
        Log21vs213NGxgraph2fdr <- data.frame(z=Log21vs213NGColfdr, Time=c("0","8","16","20","24"))
        Log21vs213NGxgraph2fdr$Time <- as.numeric(as.character(Log21vs213NGxgraph2fdr$Time))
        colnames(Log21vs213NGxgraph2fdr)[1] <- "FDR"
        NGVS_table <- Log21vs213NGxgraph2fdr
      })
      
      output$LOGNGvs_table <- renderTable({NGVS_table()}, digits = 6)
      
      
}
  
  # Run the app ----
  shinyApp(ui = ui, server = server)