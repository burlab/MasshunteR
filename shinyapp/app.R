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

library(plotly)
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
        checkboxGroupInput("QCorSample", "QC or Sample", c("All", "QC","Sample"), selected = "All"),
        radioButtons("ISTDyesorno", "ISTD?",c("All", "OnlyISTD", "NoISTD")),
        checkboxGroupInput("Variable", "Variable", c("Normalised Area", "Area", "Retention Time"), selected="Normalised Area"),
        uiOutput("selectCompound"),
        plotlyOutput("compoundPlot")
        
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
                #ISTDDetails <- read.xlsx(input$ISTDConc$datapath, sheetIndex = 2)
                ISTDDetails <- read.csv(input$ISTDConc$datapath)
                ISTDDetails$ISTD <- trimws(ISTDDetails$ISTD)
                
                dat <- tidyData(datWide)
                
                dat <- normaliseData(dat, mapISTD, ISTDDetails)
                
                
                
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
                datISTD
                
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
                
                
                output$selectCompound <- renderUI({
                  
                  if(input$ISTDyesorno=="OnlyISTD"){
                    dat <- dat[dat$isISTD,]
                    View(dat)
                  } else if(input$ISTDyesorno=="NoISTD"){
                    dat <- dat[!(dat$isISTD),]
                  } else if(input$ISTDyesorno=="All"){
                  }
                                    selectInput("CompoundList", "Select Compound", unique(dat$Compound))
                                    })
                print("finished")
                output$Status <- renderText("finished")

## NOT NEEDED ANYMORE?                
#               output$plot <- renderPlot({
#                  if(input$QCorISTD==1){
#                    plotData <- QCplot
#                    data1 <<- datQC
#                  } else {
#                    plotData <- ISTDplot
#                    data1 <<- datISTD
#                  }
#                  output$selectCompound <- renderUI({
#                    selectInput("CompoundList", "Select Compound", unique(data1$Compound))
#                  })
#                  print(plotData)
#                })
                
                output$compoundPlot <- renderPlotly({
                  data1 <- dat[dat$Compound==input$CompoundList,]
                  if("All" %in% input$QCorSample){
                  }else if("QC" %in% input$QCorSample & "Sample" %in% input$QCorSample){
                    data1 <- data1[data1$SampleType=="QC"|data1$SampleType=="Sample",]
                  } else if(input$QCorSample=="QC"){
                    data1 <- data1[data1$SampleType=="QC",]
                  } else if(input$QCorSample=="Sample"){
                    data1 <- data1[data1$SampleType=="Sample",]
                  }
                  View(data1)
                  updateSelectInput(session, "selectCompound", "Select Compound", unique(data1$Compound))

                g1 <- ggplot(data1, mapping=aes(x=AcqTime))+
                  #ggtitle("Peak ares of ISTDs in all samples") +
                  #geom_point(size=5) +
                  #geom_line(size=1) +
                  aes(ymin=0)+ 
                  #theme_scientific() +
                  #scale_y_log10() +
                  #facet_wrap(~Compound, scales="free") +
                  xlab("AcqTime") +
                  ylab("Peak Areas") +
                  theme(axis.text.x=element_blank())
                if("Area" %in% input$Variable){
                  g1 <- g1 + geom_point(aes(y = Area), size=3)
                } 
                if("Normalised Area" %in% input$Variable){
                  g1 <- g1 + geom_point(aes(y = NormArea), size=3)
                }
                if("Retention Time" %in% input$Variable){
                  g1 <- g1 + geom_point(aes(y = RT), size=3)
                }
                ggplotly(g1)
                })
                }
   )
})

# Function which takes the rawdata (wide data) and returns a tidy dataset(long data)
# This function is just a wrapper for the first part of the original logic 
tidyData <- function(rawDat){
  rawDat[1,] <- lapply(rawDat[1,],function(y) gsub(" Results","",y))
  if(rawDat[2,2]=="" & !is.na(rawDat[2,2])){
    rawDat[2,2] <- "a"
    count = 6
  } else {
    count = 5
  }
  
  # Fill in compound name in empty columns (different parameters of the same compound)
  for(c in 1:ncol(rawDat)){
    val = rawDat[1,c]
    if((!is.na(val)) && nchar(val,keepNA = TRUE) > 0 ){
      colname=val
    } else {
      rawDat[1,c]=colname
    }
  }
  
  # Concatenate rows containing parameters + compounds to the form parameter.compound for conversion with reshape() 
  rawDat[1,] <- paste(rawDat[2,], rawDat[1,],sep = ".")
  colnames(rawDat) <- rawDat[1, ]
  rawDat = rawDat[-1:-2,]
  
  # Replace column header names and remove columns that are not needed or not informative 
  setnames(rawDat, old=c(".Sample","Data File.Sample", "Name.Sample", "Acq. Date-Time.Sample", "Type.Sample"), new=c("QuantWarning","SampleFileName","SampleName", "AcqTime", "SampleTypeMethod"))
  rawDat = rawDat[, !names(rawDat) %in% c("NA.Sample","Level.Sample")]
  
  # Transform wide to long (tidy) table format
  datLong=reshape(rawDat,idvar = "SampleFileName", varying = colnames(rawDat[,-1:-count]), direction = "long",sep = "." )
  row.names(datLong) <- NULL
  setnames(datLong, old=c("time"), new=c("Compound"))
  
  # Covert to data.table object and change column types
  dat <- dplyr::tbl_df(datLong)
  rawDat <- dplyr::tbl_df(rawDat)
  numCol <- c("RT","Area","Height", "FWHM")
  dat[, names(dat) %in% numCol] <- lapply(dat[,names(dat) %in% numCol], as.numeric)
  factCol <- c("QuantWarning","SampleName","SampleFileName","SampleTypeMethod")
  dat[, names(dat) %in% factCol] <- lapply(dat[,names(dat) %in% factCol], as.factor)
  dat$Compound <- trimws(dat$Compound)
  
  return(dat)
}

# Function which takes the tidy data and calculates the Normalised area, uM and ngML
# This function is just a wrapper for the relevant code in the original logic
normaliseData <- function(longData, ISTDMapping, ISTDDetails){
  ###################################################################################################################
  # ISTD normalization and calculation of absolute concentrations
  ###################################################################################################################
  
  # Try to guess sample type based on sample file name
  longData <- longData %>% 
    mutate(SampleType=factor(ifelse(grepl("PQC", SampleName), "PQC", ifelse(grepl("TQC", SampleName), "TQC", ifelse(grepl("BLANK", SampleName), "BLANK", "Sample"))))) 
  
  # Normalises the data and adds the result to a new column
  Rprof(line.profiling = TRUE)
  longData <- as.data.table(longData)
  print("x")
  longData <- longData  %>% group_by(SampleFileName) %>% left_join(ISTDMapping[,c("Compound","ISTD")], by="Compound", copy=TRUE)# %>% 
  longData <- longData %>% mutate(isISTD = (Compound %in% ISTD)) %>%  group_by(ISTD) %>% mutate(NormArea = Area/Area[isISTD])
  #View(longData)
  #print(paste(longData$NormArea,longData$ISTD))
  Rprof(NULL)
  print("y")
  # Groups the data for later processing
  longData <- longData %>% group_by(SampleFileName)
  #debug()
  print("a")
  # Guess sample type of all runs
  longData[,SampleType:=ifelse(grepl("QC", SampleName), "QC", ifelse(grepl("BLK", SampleName), "BLANK", "Sample"))]
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
  longData[,uM := apply(longData,1,uMValue)]
  longData[,ngml := apply(longData,1,ngmlValue)]
}

# Run the application 
shinyApp(ui = ui, server = server)

