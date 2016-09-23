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

###################################################################################################################
# Import Agilent MassHunter Quant Export file (CSV) and convert to long (tidy) data format
###################################################################################################################

# Read Agilent MassHunter Quant Export file (CSV)
datWide <- read.csv("Results.csv", header = FALSE, sep = ",", na.strings=c("#N/A", "NULL"), check.names=FALSE, as.is=TRUE)
mapISTD <- read.csv("CompoundISTDList.csv", header = TRUE, sep = ",", check.names=TRUE, as.is=TRUE)
datWide[1,] <- lapply(datWide[1,],function(y) gsub(" Results","",y))

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
datLong=reshape(datWide,idvar = "SampleFileName", varying = colnames(datWide[,-1:-5]), direction = "long",sep = "." )
row.names(datLong) <- NULL
setnames(datLong, old=c("time"), new=c("Compound"))

# Covert to data.table object and change column types
dat <- dplyr::tbl_df(datLong)
numCol <- c("RT","Area","Height", "FWHM")
dat[, names(dat) %in% numCol] <- lapply(dat[,names(dat) %in% numCol], as.numeric)
factCol <- c("QuantWarning","SampleName","SampleFileName","SampleTypeMethod", "Compound")
dat[, names(dat) %in% factCol] <- lapply(dat[,names(dat) %in% factCol], as.factor)

###################################################################################################################
# ISTD normalization and calculation of absolute concentrations
###################################################################################################################

# Try to guess sample type based on sample file name
dat <- dat %>% 
  mutate(SampleType=factor(ifelse(grepl("PQC", SampleName), "PQC", ifelse(grepl("TQC", SampleName), "TQC", ifelse(grepl("BLANK", SampleName), "BLANK", "Sample"))))) 

# Normalize with corresponding ISTD, according to external data file (compound-ISTD mapping file)

#dat1 <- dat %>% group_by(SampleFileName) %>% mutate(NormArea = Compound)
dat <- dat %>% mutate(ISTD = sapply(Compound,function(y) mapISTD[which(y == mapISTD$Compound),2])) # DOES NOT WORK YET


# Write sample type (guessed based on sample name) of all samples into an additional column

# add the ISTD data to the dataset
dat <- dat %>% mutate(ISTD = sapply(Compound,function(y) mapISTD[which(y == mapISTD$Compound),2])) # DOES NOT WORK YET

# Function which takes a compound and returns it's normalised area
# com is the column name (use Area.X) and row is each row(timestamp)
normalise <- function(row){
  compo <- row[["Compound"]]
  fileName <- row[["SampleFileName"]]
  compo <- gsub("^.*\\.","",compo) # Extracts the name of the compound to reference from setnames
  istd <- mapISTD[which(mapISTD$Compound==compo),2] # Finds the relevant ISTD
  row <- filter(datWide,SampleFileName==fileName)
  # Finds the areas of the compound and the istd before calculating the normalised area
  compArea <- as.numeric(select(row,contains(paste("Area",compo,sep=".")))[,1])
  istdArea <- as.numeric(select(row,contains(paste("Area",istd,sep=".")))[,1])
  normalisedArea <- compArea/istdArea
  normalisedArea
}

# Vectors containing all the compounds and times
compoundList <- unique(dat$Compound)
timeList <- unique(dat$AcqTime)
# Normalises the data and adds the result to a new column
dat <- mutate(dat, NormArea = apply(dat,1,normalise))

# Groups the data for later processing
dat <- dat %>% group_by(SampleFileName)


# Guess sample type of all runs
dat <- dat %>% 
  mutate(SampleType=ifelse(grepl("QC", SampleName), "QC", ifelse(grepl("BLK", SampleName), "BLANK", "Sample")))


# Calculate concentrations based on spiked of ISTD (CUSTOMIZE to your data)...
# ToDo: Transfer these info to seperate input files
ISTD_CONC = 20 # ng/mL
ISTD_MW = 383.47 # g mol-1
ISTD_VOL = 50 # uL
SAMPLE_VOL = 5 # uL
dat <- dat %>% 
  mutate(uM = (NormArea   * (ISTD_VOL/1000 * ISTD_CONC/ISTD_MW*1000) / SAMPLE_VOL * 1000)/1000) %>%  
  mutate(ngml = (NormArea   * (ISTD_VOL/1000 * ISTD_CONC) / SAMPLE_VOL * 1000))



###############################
# Basic Statistics and Plots
###############################

# Estimate experimental groups, factors etc
# ------------------------------------------
# Extract groups/factors from sample names to new fields: 
# e.g. Plasma_Control_Female, Plasma_TreatmentA_Female...
# Alternatively: yet another input table containing sample information

# Assuming 3 factors... (can this be made flexible, assuming all samples names are consistent?)
expGrp = c("FactorA", "FactorB", "FactorC")

datSamples <- dat %>% filter(SampleType =="Sample") %>% 
  separate(col = SampleName, into = expGrp, convert=TRUE, remove=FALSE, sep ="-")

# necessary? (can this be made flexible, in case there are more factors, e.g. using col indices?)   
datSamples$FactorA = as.factor(pdat$FactorA)
datSamples$FactorB <- as.factor(pdat$FactorB)
datSamples$FactorC <- as.factor(pdat$FactorC)

# Basic statistics: mean +/- SD, t Test...
# ------------------------------------------

#### Work in progress....
    # Calculate average and SD of all replicates
    datSelected <- datSelected %>% group_by(LipidName, FactorA, FactorB) %>% 
      summarise(meanArea=mean(Area), SDarea = sd(Area), meanuM=mean(uM), SDuM = sd(uM))
    
    # Calculate t tests...  
    # ! Don't think this works yet...  
    
    # Filter for specific FactorA values
    filterFactorA = c("Control", "TreatmentA")
    datSelected = datSamples %>%filter(FactorA %in% filterFactorA)

    meanNormArea <- datSelected %>% group_by(LipidName, FactorB) %>% 
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
            

