#this part of the program is the importation of the datasets
#from the sanger format into a format that works for our datasets
#It will not import all data but rather it will use
library(dplyr)
library(tidyverse)
#the function for importing the data, descriptions of what everything means will be below
tttt<-list(C2BBe1_c903R1.read_count, C2BBe1_c903R2.read_count, C2BBe1_c903R3.read_count)
patr<-"/Users/user/Documents/CRISPR_thesis/gRNA Sanger"
grr<-list.files(path=patr, full.names = TRUE)
grr
reduce(full_join, by = c('sgRNA', 'gene', 'ERS717283.plasmid'))

grr <- list.files(path = patr, pattern = "\\.tsv$", full.names = TRUE)

file.info("/Users/user/Documents/CRISPR_thesis/gRNA Sanger")

ggggg<-basename(grr, read.tsv)
ggggg

prefixes <- str_extract(ggggg, "R\\d+(?=[._]|$)")
prefixes
# Step 1: Read all files into a list
file_list <- lapply(grr, read.delim)

# Step 2: Get column names of each file
all_colnames <- lapply(file_list, colnames)
all_colnames

xyz

xyz<-Reduce(function(x, y)
  merge(x, y, by = c('sgRNA', 'gene') ,all = TRUE),
  lapply(grr, read.delim))
colnames(xyz)


data<-(C2BBe1_c903R1.read_count)
#Here I am writing the normalisation
importfunct<- function(file_path, dataset) {
  #path<-
}


Imports <- function(types,  dataset) {
  
  #The type will define if the data will be imported from Achilles or from Sanger and thus what format it will use
  #originally I planned to use if else if and else however I found this function which seemed a lot easier to read and use
  switch(types,
          "sanger" = sangerfunc(),
          "achilles" =achillesfunc(),
          "Unknown data type please use sanger or achilles as formats"
  )
  
  #Sanger data importation you can see above the code, I will move it here later, I needed it outside of the function
  #to experiment, it is meant to call upon this function. I might change it to make it modular but it is intended for
  #our specific dataset by defining it by the -plasmid in the name rather than the full name.
  sangerfunc<-(dataset) {
    
  }
  
  #Achilles data importation once we are both done we can merge the function this more or less exists as a place holder
  achillesfunc
  # OBSOLETE storage, depending on the storage folder or file it will know if it will just extract the data from the file names
  #or a complete folder, name will be the the name of this file or folder
  
}
##FDR P values
           