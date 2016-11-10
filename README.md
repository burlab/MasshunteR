# MasshunteR

Processing of Agilent MassHunter Quant results
--------------------------------------------------------

This script aims to process Agilent MassHunter Quant Results:
- Imports the output table (csv) of the Agilent MassHunter Quant Software  (containing peak areas, RT, FWHM etc)
- Imports a text file (csv) with compound - ISTD mappings
- Normalizes peak areas with ISTD  
- Calculates concentrations based spiked ISTD concentration/amount
- Predicts sample type (sample, QC, Blank) based on sample file name   
- Plots some QC charts


Bo Burla, Singapore Lipidomics Incubator (SLING), National University of Singapore

##Shiny
All of the Above functionality is implemented through a shiny application (found in shinyapp\). This shiny application features interactive (and colourblind friendly) plots through plotly as well as some additional features. Specifically:
- csv files containing the normalised areas and concentrations can be downloaded seperately
- a csv file with the summarised data (mean and sd) can also be downloaded

Adithya Diddapur, Singapore Lipidomics Incubator (SLING), National University of Singapore
