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
 
Bo Burla & Adithya Diddapur

Singapore Lipidomics Incubator (SLING), National University of Singapore
