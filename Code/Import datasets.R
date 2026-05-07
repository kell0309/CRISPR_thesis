# this part of the pipeline is the importation of the datasets
# The original format is sanger but has individual tsv files, this part will do the following:
# 1.Merge all of the into tables based on their shared Plasmid as it is a control
# 2.Merge results that come from the same cell line
# 3.Additional variable was found time, it took priority over cell lines
# this part of the repository is made by Vasileios Theocharis

#These are libraries that are used the two with # are ones that were tested but not used
library(dplyr)

#This is the path to the data on my computers
#patr<-"/Users/user/Documents/CRISPR_thesis/gRNA Sanger"
patr<-"/home/vasilis/SchoolFolder/bioinformatics/R scripts/CRISPR_thesis/gRNA Sanger"
grr<-list.files(path=patr, full.names = TRUE)

#here I extract the column names and add the names of the files back
#this is meant for the grouping
sep_extract <- lapply(grr, function(f) {
  colnames(read.delim(f, nrows = 1, header = TRUE, sep = "\t"))
})
#separated on the 4th column meaning the sample
#we are only extracting the column heads as this is about grouping said names
names(sep_extract) <- basename(grr)
sep_final<- split(sep_extract, sapply(sep_extract, `[`, 4))
#this part is the is now taking the heads from the extraction and matching them 
#match them based on the R which is the runs and then based on the _ which connects the names
grouped_R_ <- lapply(sep_final, function(group) {
  #Most of our data is only technical replicates however some are also time based and with unknown conditions
  #this is a regex made to define both the normal format of our data for recognition and
  #the two groups the regex specifically only rearches for the simple format
  grouped_time <- ifelse(grepl("R\\d+\\.read_count\\.tsv$", names(group)), "simple", "time")
  #here we go through the data and take up only the names rather than
  #keeping the whole file name
  lapply(split(group, grouped_time), function(time) {
    groupedR <- gsub("R\\d+(_[A-Z0-9]+)?\\.read_count\\.tsv$", "\\1", names(time))
  #Here is the simpler group only defined by the replicates
  #groupedR <- gsub("R\\d+\\.read_count\\.tsv$", "", names(group))
  
  #the processing depends on there being multiples we before it is returned ones 
  #with only a single replicate are removed
    grouped_ <- gsub("_$", "", groupedR)
    split(time, grouped_)
  })
})
#extracting out the data out of the resutls into the two groups
time_non<- lapply(grouped_R_, `[[`, "simple")
#time_pos<- lapply(grouped_R_,  `[[`, "time")
#after this runs we now have all the data that is part of the same cell line
#and use the same plasmid together

#This part here uses these lists to  merge them together It will go through each of these list
call_merge<-function(time, patr) {
  
  lapply(time, function(av_pla) {
    av_funct<-lapply(av_pla, function(av_line) {
      replicates <- lapply(names(av_line), function(name) {
        #accessing the file based on path
        read.delim(file.path(patr, name), header = TRUE)
      })
      #mean rows to create
      mean_sample <- rowMeans(sapply(replicates, function(x) x[, 3]))
      plasmid <- replicates[[1]][, 4]
      #creating a new dataframe with each of the cell lines that have been merged
      #by renaming them to plasmid it makes it so we no longer need to account for the names
      #up to now names have been an issue as they have been different between the two groups
      data.frame(sgRNA = replicates[[1]][, 1], gene = replicates[[1]][, 2], mean_sample = mean_sample, plasmid = plasmid)
    })
    #Last part, we use reduce and map. reduce to iterate within the lists and map to iterate the tables in the lists
    #by= is the three common columns among our dataset thus they are used as consistent variables
    red_merge<-Reduce(function(x,y) merge(x,y, by = c("sgRNA", "gene", "plasmid")), 
                      #map takes the column and gives it the name of the dataframe to preserve the structure
                      Map(function(f_merge, name) {
                        name->colnames(f_merge)[3] #Abnormal syntax, personal preference
                        #this is where everything is iterated over during the function
                        f_merge
                        #close the function, return the variables
                      }, av_funct, names(av_funct)))
    # here the low already essential genes are removed if they already die in the control
    #var<-red_merge[red_merge$plasmid > 10, ] 
  })
}

#running the function in order to process the data and merge all the means
#if identical sample cell lines
non_merge<-call_merge(time_non, patr)
#pos_merge<-call_merge(time_pos, patr)

#removing from the list the two different controls, in order to be processed further.
list2env(setNames(non_merge, paste0("non_", names(non_merge))), envir = .GlobalEnv)
#list2env(setNames(pos_merge, paste0("pos_", names(pos_merge))), envir = .GlobalEnv)




           