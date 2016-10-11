#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

###################################################################################################################
#
# CONSTANTS
#--------------------------------------------------------
#
# Calculate concentrations based on spiked of ISTD (CUSTOMIZE to your data)...
# ToDo: Transfer these info to seperate input files
ISTD_VOL = 50 # uL
SAMPLE_VOL = 5 # uL

# Constants used to split the data and perform statistical analysis (t.test so far)
# ToDo : Make these more flexible (different numbers of parameters as needed)
expGrp = c("ParameterA", "ParameterB", "ParameterC")
filterParameterA = c("ACTH") #Vector include any value
#
####################################################################################################################

library(shiny)
library(artyfarty)
library(highcharter)
library(broom)
library(tidyr)
library(dplyr)
library(dtplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(xlsx)

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
   
   # Application title
   titlePanel("masshunteR"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         fileInput("mainInput", "mainINPUT"),
         fileInput("ISTDMapping", "ISTDMapping"),
         fileInput("ISTDConc", "ISTDConc"),
         actionButton("Submit", "Submit"),
         verbatimTextOutput("Status"),
         downloadButton("DownloadData", "DownloadData"),
         downloadButton("DownloadISTD", "DownloadISTD"),
         downloadButton("DownloadQC", "DownloadQC")),
      
      # Show a plot of the generated distribution
      mainPanel(
        selectInput("QCorISTD", "QC or ISTD", choices=list("QC" = 1, "ISTD" = 2)),
        uiOutput("selectCompound"),
        plotOutput("compoundPlot"),
        plotOutput("plot")
        
      )
   )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output, session) {
  session$onSessionEnded(stopApp)
   observeEvent(input$Submit,
                {
                  output$Status <- renderText({"processing"})
                  #output$Status <- renderText("Processing")
                print("Processing")
                #datWide <- read.csv("20160615_Pred_ACTH_Original_Data.csv", header = FALSE, sep = ",", na.strings=c("#N/A", "NULL"), check.names=FALSE, as.is=TRUE, strip.white=TRUE )
                #mapISTD <- read.csv("CompoundISTDList_SLING-PL-Panel_V1.csv", header = TRUE, sep = ",", check.names=TRUE, as.is=TRUE, strip.white = TRUE)
                #ISTDDetails <- read.xlsx("ISTD-map-conc_SLING-PL-Panel_V1.xlsx", sheetIndex = 2)
                datWide <- read.csv(input$mainInput$datapath, header=FALSE, sep=",", na.strings=c("#N/A", "NULL"), check.names=FALSE, as.is=TRUE, strip.white=TRUE)
                mapISTD <- read.csv(input$ISTDMapping$datapath, header = TRUE, sep = ",", check.names=TRUE, as.is=TRUE, strip.white = TRUE)  
                ISTDDetails <- read.xlsx(input$ISTDConc$datapath, sheetIndex = 2)
                #ISTDDetails <- read.csv(input$ISTDConc$datapath)
                ISTDDetails$ISTD <- trimws(ISTDDetails$ISTD)
                
                datWide[1,] <- lapply(datWide[1,],function(y) gsub(" Results","",y))
                if(datWide[2,2]=="" & !is.na(datWide[2,2])){
                  datWide[2,2] <- "a"
                  count = 6
                } else {
                  count = 5
                }
                
                # Fill in compound name in empty columns (different parameters of the same compound)
                for(c in 1:ncol(datWide)){
                  val = datWide[1,c]
                  if((!is.na(val)) && nchar(val,keepNA = TRUE) > 0 ){
                    colname=val
                  } else {
                    datWide[1,c]=colname
                  }
                }
                
                # Concatenate rows containing parameters + compounds to the form parameter.compound for conversion with reshape() 
                datWide[1,] <- paste(datWide[2,], datWide[1,],sep = ".")
                colnames(datWide) <- datWide[1, ]
                datWide = datWide[-1:-2,]
                
                # Replace column header names and remove columns that are not needed or not informative 
                setnames(datWide, old=c(".Sample","Data File.Sample", "Name.Sample", "Acq. Date-Time.Sample", "Type.Sample"), new=c("QuantWarning","SampleFileName","SampleName", "AcqTime", "SampleTypeMethod"))
                datWide = datWide[, !names(datWide) %in% c("NA.Sample","Level.Sample")]
                
                # Transform wide to long (tidy) table format
                datLong=reshape(datWide,idvar = "SampleFileName", varying = colnames(datWide[,-1:-count]), direction = "long",sep = "." )
                row.names(datLong) <- NULL
                setnames(datLong, old=c("time"), new=c("Compound"))
                
                # Covert to data.table object and change column types
                dat <- dplyr::tbl_df(datLong)
                datWide <- dplyr::tbl_df(datWide)
                numCol <- c("RT","Area","Height", "FWHM")
                dat[, names(dat) %in% numCol] <- lapply(dat[,names(dat) %in% numCol], as.numeric)
                factCol <- c("QuantWarning","SampleName","SampleFileName","SampleTypeMethod")
                dat[, names(dat) %in% factCol] <- lapply(dat[,names(dat) %in% factCol], as.factor)
                dat$Compound <- trimws(dat$Compound)
                ###################################################################################################################
                # ISTD normalization and calculation of absolute concentrations
                ###################################################################################################################
                
                # Try to guess sample type based on sample file name
                dat <- dat %>% 
                  mutate(SampleType=factor(ifelse(grepl("PQC", SampleName), "PQC", ifelse(grepl("TQC", SampleName), "TQC", ifelse(grepl("BLANK", SampleName), "BLANK", "Sample"))))) 
                
                # Normalises the data and adds the result to a new column
                Rprof(line.profiling = TRUE)
                dat <- as.data.table(dat)
                print("x")
                dat <- dat  %>% group_by(SampleFileName) %>% left_join(mapISTD[,c("Compound","ISTD")], by="Compound", copy=TRUE)# %>% 
                dat <- dat %>% mutate(isISTD = (Compound %in% ISTD)) %>%  group_by(ISTD) %>% mutate(NormArea = Area/Area[isISTD])
                View(dat)
                print(paste(dat$NormArea,dat$ISTD))
                Rprof(NULL)
                print("y")
                # Groups the data for later processing
                dat <- dat %>% group_by(SampleFileName)
                #debug()
                print("a")
                # Guess sample type of all runs
                dat[,SampleType:=ifelse(grepl("QC", SampleName), "QC", ifelse(grepl("BLK", SampleName), "BLANK", "Sample"))]
                print("z")
                # Functions to calculate the concentrations and then add them to dat
                # Each function takes an entire row from dat as input
                uMValue <- function(row){
                  istd <- row[["ISTD"]]
                  ISTD_CONC <- ISTDDetails[ISTDDetails$ISTD==istd,"ISTDconcNGML"]
                  ISTD_MW <- ISTDDetails[ISTDDetails$ISTD==istd,"ISTD_MW"]
                  normalisedArea <- row[["NormArea"]]
                  umVal <- (as.numeric(normalisedArea)   * (ISTD_VOL/1000 * ISTD_CONC/ISTD_MW*1000) / SAMPLE_VOL * 1000)/1000
                  umVal
                }
                
                ngmlValue <- function(row){
                  istd <- row[["ISTD"]]
                  ISTD_CONC <- ISTDDetails[ISTDDetails$ISTD==istd,"ISTDconcNGML"]
                  ISTD_MW <- ISTDDetails[ISTDDetails$ISTD==istd,"ISTD_MW"]
                  normalisedArea <- row[["NormArea"]]
                  ngmlVal <- as.numeric(normalisedArea)   * (ISTD_VOL/1000 * ISTD_CONC) / SAMPLE_VOL * 1000
                  ngmlVal
                }
                
                # <- mutate(dat, uM = uMValue(ISTD,NormArea))
                dat[,uM := apply(dat,1,uMValue)]
                dat[,ngml := apply(dat,1,ngmlValue)]
                
                
                ###############################
                # Basic Statistics and Plots
                ###############################
                
                # Estimate experimental groups, factors etc
                # ------------------------------------------
                # Extract groups/factors from sample names to new fields: 
                # e.g. Plasma_Control_Female, Plasma_TreatmentA_Female...
                # Alternatively: yet another input table containing sample information
                
                # Assuming 3 factors... (can this be made flexible, assuming all samples names are consistent?)
                #expGrp = c("FactorA", "FactorB", "FactorC")
                
                datSamples <- dat %>% filter(SampleType =="Sample")  
                datSamples <- separate(datSamples,col = SampleName, into = expGrp, convert=TRUE, remove=FALSE, sep ="-")
                datSamples$ParameterA <- gsub("^.*?_","",datSamples[[4]])
                
                # necessary? (can this be made flexible, in case there are more factors, e.g. using col indices?)   
                datSamples$ParameterA <- as.factor(datSamples$ParameterA)
                datSamples$ParameterB <- as.factor(datSamples$ParameterB)
                datSamples$ParameterC <- as.factor(datSamples$ParameterC)
                
                # Basic statistics: mean +/- SD, t Test...
                # ------------------------------------------
                
                #### Work in progress....
                
                # Wrapper function for t.test p-value which returns NA instead of an error if the data is invalid
                # e.g. insufficient data points now return NA instead of throwing an error
                # Function by Tony Plate at https://stat.ethz.ch/pipermail/r-help/2008-February/154167.html
                my.t.test.p.value <- function(...) {
                  obj<-try(t.test(...,paired=TRUE), silent=TRUE)
                  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
                }
                
                
                #function to calculate p-value given dataframe from a single group
                pValueFromGroup <- function(data){
                  bValues <- unique(data$ParameterB)
                  if(!(length(bValues==2))){
                    stop("length(bValues)!=2")
                  }
                  dat1 <- data[data$ParameterB==bValues[1],]
                  dat2 <- data[data$ParameterB==bValues[2],]
                  pValue <- my.t.test.p.value(dat1$NormArea,dat2$NormArea)
                  #pValue <- t.test(dat1$NormArea,dat2$NormArea)$p.value
                  pValue
                }
                
                meanNormArea <- datSamples %>% filter(ParameterA %in% filterParameterA) %>%
                  filter(NormArea!=1)
                #group_by(Compound,ParameterB) #%>%
                #pVal <- by(meanNormArea, as.factor(meanNormArea$Compound),pValueFromGroup, simplify = TRUE)
                
                #datFiltered <- datSamples %>% 
                #  filter(ParameterA %in% filterParameterA) %>%
                #  filter(grep("LPC 20:1",Compound)) %>%
                #  filter(NormArea!=1) %>%
                #  droplevels() %>%
                #  group_by(Compound) %>%
                #  do(tidy(t.test(uM~ParameterB,data=., paired=TRUE)))
                
                # Calculate average and SD of all replicates
                datSelected <- datSamples %>% group_by(Compound, ParameterA, ParameterB) %>% 
                  summarise(meanNormArea=mean(NormArea), SDNormarea = sd(NormArea), meanuM=mean(uM), SDuM = sd(uM), nArea = n()) %>%
                  filter(ParameterA %in% filterParameterA) %>%
                  filter(meanNormArea!=1)# %>%
                #  mutate(pValue = pVal[[Compound]])
                print("done")
                print("test")
                output$DownloadData <- downloadHandler(
                  filename='data.csv',
                  content=function(file) {
                    write.csv(datSelected, file)
                  }
                )
                # Plot peak areas of compounds in all QC samples
                # --------------------------------------------------     
                
                datQC <- dat[SampleType=="QC"]
                
                QCplot <- ggplot(data=datQC, mapping=aes(x=AcqTime,y=NormArea, group=1, ymin=0)) +
                  ggtitle("Peak Areas of QC samples") +
                  geom_point(size=0.8) +
                  geom_line(size=1) +
                  #scale_y_log10() +
                  facet_wrap(~Compound, scales="free") +
                  xlab("AcqTime") +
                  ylab("Peak Areas") +
                  theme(axis.text.x=element_blank())# +
                  #ggsave("QCplot.png",width=30,height=30) 
                output$DownloadQC <- downloadHandler(
                  filename='QCplot.png',
                  content=function(file){
                    ggsave(file,plot=QCplot,width=30,height=30)
                  }
                )
                
                
                # Plot peak areas of ISTDs in all samples, colored by sampleType
                # --------------------------------------------------------------     
                
                
                datISTD <- dat[grepl("(IS)",Compound),]         
                
                ISTDplot <- ggplot(data=datISTD, mapping=aes(x=AcqTime,y=NormArea,color=SampleType, group=1, ymin=0))+
                  ggtitle("Peak ares of ISTDs in all samples") +
                  geom_point(size=0.8) +
                  geom_line(size=1) +
                  #scale_y_log10() +
                  facet_wrap(~Compound, scales="free") +
                  xlab("AcqTime") +
                  ylab("Peak Areas") +
                  theme(axis.text.x=element_blank()) #+
                  #ggsave("ISTDplot.png",width=30,height=30)
                output$DownloadISTD <- downloadHandler(
                  filename='ISTDplot.png',
                  content=function(file){
                    ggsave(file, plot=ISTDplot, width=30,height=30)
                  }
                )
                
                
                print("finished")
                output$Status <- renderText("finished")
                output$plot <- renderPlot({
                  if(input$QCorISTD==1){
                    plotData <- QCplot
                    data1 <<- datQC
                  } else {
                    plotData <- ISTDplot
                    data1 <<- datISTD
                  }
                  output$selectCompound <- renderUI({
                    selectInput("CompoundList", "Select Compound", unique(data1$Compound))
                  })
                  print(plotData)
                })
                
                output$compoundPlot <- renderPlot({
                  data1 <- data1[data1$Compound==input$CompoundList,]

                g1 <- ggplot(data1, mapping=aes(x=AcqTime,y=NormArea, group=1, ymin=0))+
                  ggtitle("Peak ares of ISTDs in all samples") +
                  geom_point(size=0.8) +
                  geom_line(size=1) +
                  theme_scientific() +
                  #scale_y_log10() +
                  #facet_wrap(~Compound, scales="free") +
                  xlab("AcqTime") +
                  ylab("Peak Areas") +
                  theme(axis.text.x=element_blank())
                if(input$QCorISTD!=1){
                  g1 <- g1 + aes(color=SampleType)
                }
                print(g1)
                })
                }
   )
#  output$compoundPlot <- renderPlot({
#    input$plotGraph
#    isolate(
#      data1 <- data1[data1$Compound==input$CompoundList,]
#    )
#    
#    g1 <- ggplot(data1, mapping=aes(x=AcqTime,y=NormArea, group=1, ymin=0))+
#      ggtitle("Peak ares of ISTDs in all samples") +
#      geom_point(size=0.8) +
#      geom_line(size=1) +
#      theme_scientific() +
#      #scale_y_log10() +
#      #facet_wrap(~Compound, scales="free") +
#      xlab("AcqTime") +
#      ylab("Peak Areas") +
#      theme(axis.text.x=element_blank())
#    if(input$QCorISTD!=1){
#      g1 <- g1 + aes(color=SampleType)
#    }
#    print(g1)
#  })
  

})

# Run the application 
shinyApp(ui = ui, server = server)

