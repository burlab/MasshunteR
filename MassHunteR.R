###################################################################################################################
#
# Processing of Agilent MassHunter Quant results
# --------------------------------------------------------
#
# This script processes Agilent MassHunter Quant Results:
#    - Imports the output table (csv) of the Agilent MassHunter Quant Software  (containing peak areas, RT, FWHM etc)
#    - Imports a text file (csv) with compound - ISTD mappings
#    - Normalizes peak areas with ISTD  
#    - Calculates concentrations based spiked ISTD concentration/amount
#    - Predicts sample type (sample, QC, Blank) based on sample file name   
#    - Plots some QC charts
# 
#
# Bo Burla / Singapore Lipidomics Incubator (SLING)
# National University of Singapore
#
# 13.06.2016 -
#
###################################################################################################################

#setwd("D:/Bo/Data/RawData/LCMS/ExperimentA")
#setwd("D:/Bo/Data/RawData/LCMS/GL08D_StabilityTests")
setwd("D://Adithya//Sample data")

library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(xlsx)

###################################################################################################################
# Import Agilent MassHunter Quant Export file (CSV) and convert to long (tidy) data format
###################################################################################################################

# Read Agilent MassHunter Quant Export file (CSV)
datWide <- read.csv("20160615_Pred_ACTH_Original_Data.csv", header = FALSE, sep = ",", na.strings=c("#N/A", "NULL"), check.names=FALSE, as.is=TRUE, strip.white=TRUE )
mapISTD <- read.csv("CompoundISTDList_SLING-PL-Panel_V1.csv", header = TRUE, sep = ",", check.names=TRUE, as.is=TRUE, strip.white = TRUE)
ISTDDetails <- read.xlsx("ISTD-map-conc_SLING-PL-Panel_V1.xlsx", sheetIndex = 2)
ISTDDetails$ISTD <- trimws(ISTDDetails$ISTD)
#datWide <- read.csv("Results.csv", header = FALSE, sep = ",", na.strings=c("#N/A", "NULL"), check.names=FALSE, as.is=TRUE)
#mapISTD <- read.csv("CompoundISTDList.csv", header = TRUE, sep = ",", check.names=TRUE, as.is=TRUE, strip.white = TRUE)

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

# Normalize with corresponding ISTD, according to external data file (compound-ISTD mapping file)

# Write sample type (guessed based on sample name) of all samples into an additional column

# add the ISTD data to the dataset
dat <- dat %>% mutate(ISTD = sapply(Compound,function(y) mapISTD[which(y == mapISTD$Compound),2])) 
print("x")

# Function which takes a compound and returns it's normalised area
# input is a complete row from the data frame

normalise <- function(row){
    compo <- trimws(row[["Compound"]])
    fileName <- row[["SampleFileName"]]
    istd <- row[["ISTD"]]
    compArea <- as.numeric(row[["Area"]])
    istdArea <- as.numeric(dat[dat$SampleFileName == fileName & dat$Compound == istd,][["Area"]])
    normalisedArea <- compArea / istdArea
    normalisedArea
}

# Normalises the data and adds the result to a new column
Rprof(line.profiling = TRUE)
#x <- apply(dat,1, normalise)
dat <- as.data.table(dat)
dat[,NormArea := apply(dat,1,normalise)]
#dat <- mutate(dat, NormArea = apply(dat,1,normalise))
Rprof(NULL)

# Groups the data for later processing
#dat <- dat %>% group_by(SampleFileName)

# Guess sample type of all runs
dat[,SampleType:=ifelse(grepl("QC", SampleName), "QC", ifelse(grepl("BLK", SampleName), "BLANK", "Sample"))]
#dat <- dat %>% 
#  mutate(SampleType=ifelse(grepl("QC", SampleName), "QC", ifelse(grepl("BLK", SampleName), "BLANK", "Sample")))


# Calculate concentrations based on spiked of ISTD (CUSTOMIZE to your data)...
# ToDo: Transfer these info to seperate input files
#ISTD_CONC = 20 # ng/mL
#ISTD_MW = 383.47 # g mol-1
ISTD_VOL = 50 # uL
SAMPLE_VOL = 5 # uL

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
print("y")
expGrp = c("FactorA", "FactorB", "FactorC")
print("z")

datSamples <- dat %>% filter(SampleType =="Sample")  
print("a")
datSamples <- separate(datSamples,col = SampleName, into = expGrp, convert=TRUE, remove=FALSE, sep ="-")
print("x")

# necessary? (can this be made flexible, in case there are more factors, e.g. using col indices?)   
datSamples$FactorA = as.factor(datSamples$FactorA)
datSamples$FactorB <- as.factor(datSamples$FactorB)
datSamples$FactorC <- as.factor(datSamples$FactorC)

# Basic statistics: mean +/- SD, t Test...
# ------------------------------------------

#### Work in progress....
    # Calculate average and SD of all replicates
    datSelected <- datSamples %>% group_by(Compound, FactorA, FactorB) %>% 
      summarise(meanArea=mean(Area), SDarea = sd(Area), meanuM=mean(uM), SDuM = sd(uM))
    
    # Calculate t tests...  
    # ! Don't think this works yet...  
    
    # Filter for specific FactorA values
    filterFactorA = c("Control", "TreatmentA")
    datSelected = datSamples %>%filter(FactorA %in% filterFactorA)

    meanNormArea <- datSelected %>% group_by(Compound, FactorB) %>% 
      summarise_each(funs(t.test(.[vs == 0], .[vs == 1])$p.value), vars = disp:qsec)
#### ......

    
# --------------------------------
# Plots
# --------------------------------
# Plot concentrations vs FactorA for each compound, different line and colored according to FactorC (one compound per panel) 

ggplot(data=pdat, mapping=aes(x = FactorB, y = uM, group=FactorC, color = FactorC)) +
  ggtitle("Treatment") +
  geom_point(size = 3) +
  geom_line(size=0.8)  +
  scale_colour_brewer(palette = "Set1") +
  theme_grey(base_size = 10) +
  facet_wrap(~LipidName, scales="free") +
  aes(ymin=0) +
  #geom_errorbar(aes(ymax = meanConc + SDfmol, ymin=meanConc - SDfmol), width=1)  +
  #geom_smooth(method='lm', se = FALSE, level=0.95)
  xlab("Days under  treatment") +
  ylab("uM in plasma") + 
  theme(axis.text=element_text(size=9), axis.title=element_text(size=12,face="bold"), 
        strip.text = element_text(size=10, face="bold"),
        legend.title=element_text(size=10, face="bold"),
        #legend.position=c(0.89,0.1),
        plot.title = element_text(size=16, lineheight=2, face="bold", margin=margin(b = 20, unit = "pt"))) +
        annotate("text", x = 1.5, y = 1, label = "Some text")

      
###############################
# QC Plots
###############################      
      
# Plot retention time of all compounds in all samples
# --------------------------------------------------
 
datSelected = datSamples %>%filter(SampleType %in% c("BLANK"))    
       
ggplot(data=datSelected, mapping=aes(x = SampleName, y = RT, color = SampleType)) +
  ggtitle("Retention Time") +
  geom_point(size = 3) +
  geom_line(size=0.8)  +
  scale_colour_brewer(palette = "Set1") +
  theme_grey(base_size = 10) +
  facet_wrap(~LipidName, scales="free") +
  aes(ymin=0) +
  xlab("Sample") +
  ylab("Retention time [min]") + 
  theme(axis.text=element_text(size=9), axis.title=element_text(size=12,face="bold"), 
        strip.text = element_text(size=10, face="bold"),
        legend.title=element_text(size=10, face="bold"),
        #legend.position=c(0.89,0.1),
        plot.title = element_text(size=16, lineheight=2, face="bold", margin=margin(b = 20, unit = "pt"))) +
        annotate("text", x = 1.5, y = 1, label = "Some text")


        
# Plot peak areas of compounds in all QC samples
# --------------------------------------------------     


# Plot peak areas of ISTDs in all samples, colored by sampleType
# --------------------------------------------------------------     
            

