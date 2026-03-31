# this part of the program is the importation of the datasets
# The original format is sanger but has individual tsv files, this part will do the following:
# 1.merge all of the into tables based on their shared Plasmid as it is a control
# 2.Merge results that come from the same cell line
# variable elements sep_ means it is extracting and grouping the data
# grouped are part of the grouping of the replicates of the cell lines
# av_ are the variables that are used as part of merging the replicant cell lines into a single variable
# as they they are the same cell line we cannot use them as separate variables without swaying the dataset
# this part of the repository is made by Vasileios Theocharis


#These are libraries that are used the two with # are ones that were tested but not used
#library(dplyr)
#library(tidyverse)
library(edgeR)

#This is the path to the data on my computer, in the fuction it will be made as a variable
patr<-"/Users/user/Documents/CRISPR_thesis/gRNA Sanger"
grr<-list.files(path=patr, full.names = TRUE)


# here I extract the column names and add the names of the files back
sep_extract <- lapply(grr, function(f) {
  colnames(read.delim(f, nrows = 1, header = TRUE, sep = "\t"))
})
#separated on the 4th column
names(sep_extract) <- basename(grr)
sep_final<- split(sep_extract, sapply(sep_extract, `[`, 4))

#this part is the is now taking the heads from the extraction and matching them 
#match them based on the R which is the runs and then based on the _ which connects the names
grouped_R_ <- lapply(sep_final, function(group) {
  groupedR <- gsub("R\\d+\\.read_count\\.tsv$", "", names(group))
  grouped_ <- sapply(strsplit(groupedR, "_"), `[`, 1)
  split(group, grouped_)
})
#after this runs we now have all the data that is part of the same cell line and use the same plasmid together


#This part here uses these lists to finally merge them together It will go through each of these list

av_merge <- lapply(grouped_R_, function(av_pla) {
  av_funct<-lapply(av_pla, function(av_line) {
    replicates <- lapply(names(av_line), function(name) {
      read.delim(file.path(patr, name), header = TRUE)
    })
    #mean rows to create
    mean_sample <- rowMeans(sapply(replicates, function(x) x[, 3]))
    plasmid <- replicates[[1]][, 4]
    #creating a new dataframe with each of the cell lines that have been merged
    #by renaming them to plasmid it makes it so we no longer need to account for the names
    # up to now names have been an issue as they have been different between the two groups
    data.frame(sgRNA = replicates[[1]][, 1], gene = replicates[[1]][, 2], mean_sample = mean_sample, plasmid = plasmid)
  })
  #Last part, we use reduce and map reduce to iterate within the the lists and map to iterate the tables in the lists
  Reduce(function(x,y) merge(x,y, by = c("sgRNA", "gene", "plasmid")), 
  #map takes the column and gives it the name of the dataframe to preserve the structure
    Map(function(f_merge, name) {
      name->colnames(f_merge)[3] #Abnormal syntax, personal preference because name is the new variable
      f_merge
    }, av_funct, names(av_funct)))
})

#For normalisation here we will use TMM from EdgeR we have to extract the data into the environment
list2env(av_merge, envir = .GlobalEnv)
#for TMM we need groups to define the plasmid

#Making the count groups for TMM in order to set plasmid as control
group <-factor(ifelse(grepl("plasmid", 
#Sapply to make it so the columns have to be numberic, this removes the non count columns
  colnames(ERS717283.plasmid[sapply(ERS717283.plasmid, is.numeric)])), "plasmid", "sample"))


con_DGE<-DGEList(counts = ERS717283.plasmid, group = factor(ifelse(grepl("plasmid", 
#Sapply to make it so the columns have to be numberic, this removes the non count columns
  colnames(ERS717283.plasmid[sapply(ERS717283.plasmid, is.numeric)])), "plasmid", "sample")))

con_norm<- calcNormFactors(con_DGE, method = "TMM", refColumn = 1)

Imports <- function(types,  datapath) {
  
  #The type will define if the data will be imported from Achilles or from Sanger and thus what format it will use
  #originally I planned to use if else if and else however I found this function which seemed a lot easier to read and use
  switch(types,
          "sanger" = sangerfunc(),
          "achilles" =achillesfunc(),
          "Unknown data type or file path please use sanger or achilles data or correct your datapath"
  )
  
  #Sanger data importation you can see above the code, I will move it here later, I needed it outside of the function
  #to experiment, it is meant to call upon this function. I might change it to make it modular but it is intended for
  #our specific dataset by defining it by the -plasmid in the name rather than the full name.
  sangerfunc<-(datapath) {
    patr<-datapath
    grr<-list.files(path=patr, full.names = TRUE)
    
    
    # here I extract the column names and add the names of the files back
    sep_extract <- lapply(grr, function(f) {
      colnames(read.delim(f, nrows = 1, header = TRUE, sep = "\t"))
    })
    #separated on the 4th column
    names(sep_extract) <- basename(grr)
    sep_final<- split(sep_extract, sapply(sep_extract, `[`, 4))
    
    #this part is the is now taking the heads from the extraction and matching them 
    #match them based on the R which is the runs and then based on the _ which connects the names
    grouped_R_ <- lapply(sep_final, function(group) {
      groupedR <- gsub("R\\d+\\.read_count\\.tsv$", "", names(group))
      grouped_ <- sapply(strsplit(groupedR, "_"), `[`, 1)
      split(group, grouped_)
    })
    #after this runs we now have all the data that is part of the same cell line and use the same plasmid together
    
    
    #This part here uses these lists to finally merge them together It will go through each of these list
    
    av_merge <- lapply(grouped_R_, function(av_pla) {
      av_funct<-lapply(av_pla, function(av_line) {
        replicates <- lapply(names(av_line), function(name) {
          read.delim(file.path(patr, name), header = TRUE)
        })
        #mean rows to create
        mean_sample <- rowMeans(sapply(replicates, function(x) x[, 3]))
        plasmid <- replicates[[1]][, 4]
        #creating a new dataframe with each of the cell lines that have been merged
        #by renaming them to plasmid it makes it so we no longer need to account for the names
        # up to now names have been an issue as they have been different between the two groups
        data.frame(sgRNA = replicates[[1]][, 1], gene = replicates[[1]][, 2], mean_sample = mean_sample, plasmid = plasmid)
      })
      #Last part, we use reduce and map reduce to iterate within the the lists and map to iterate the tables in the lists
      Reduce(function(x,y) merge(x,y, by = c("sgRNA", "gene", "plasmid")), 
             #map takes the column and gives it the name of the dataframe to preserve the structure
             Map(function(f_merge, name) {
               name->colnames(f_merge)[3] #Abnormal syntax, personal preference because name is the new variable
               f_merge
             }, av_funct, names(av_funct)))
    })
    list2env(av_merge, envir = .GlobalEnv)
    
    
  }
  
  #Achilles data importation once we are both done we can merge the function this more or less exists as a place holder
  achillesfunc
  # OBSOLETE storage, depending on the storage folder or file it will know if it will just extract the data from the file names
  #or a complete folder, name will be the the name of this file or folder
  
}
##FDR P values
           