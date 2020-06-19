## by Simon Schwab 2019, 2020

library(testit)

# Reads the data from single studies
pb.readStudies = function(path, file) {
  
  data = read.csv(file.path(path, file),
                  colClasses=c("factor",                             # id
                               "integer", "character", "character",  # comparison :  nr, name, id
                               "integer", "character", "character",  # outcome    :  nr, name, measure
                               "character", "factor",                #               id, flag
                               "integer",  "character", "character", # subgroup   :  nr, name, id
                               "character","character","character",  # study      :  id, name, year
                               "factor",                             #               data_source  
                               "numeric", "numeric",                 # effect, std.err.
                               "numeric", "numeric",                 # group 1    :  events, total
                               "numeric", "numeric",                 # group 1    :  mean, sd
                               "numeric", "numeric",                 # group 2    :  events, total
                               "numeric", "numeric"                  # group 2    :  mean, sd
                  ))
  
  # Fix invalid years from study.year
  idx_invalid = (!grepl("^(19|20)\\d{2}$", data$study.year, perl = TRUE))
  idx_hasyear = grepl(".*(19|20)(\\d{2}).*", data$study.year, perl = TRUE)
  fixme = data$study.year[idx_invalid & idx_hasyear]
  fixed = gsub(".*(19|20)(\\d{2}).*", "\\1\\2", fixme, perl = TRUE)
  data$study.year[idx_invalid & idx_hasyear] = fixed
  
  # Fix remaining invalid years from study.id
  idx_invalid = (!grepl("^(19|20)\\d{2}$", data$study.year, perl = TRUE))
  idx_hasyear = grepl(".*(19|20)(\\d{2}).*", data$study.id, perl = TRUE)
  fixme = data$study.id[idx_invalid & idx_hasyear]
  fixed = gsub(".*(19|20)(\\d{2}).*", "\\1\\2", fixme, perl = TRUE)
  data$study.year[idx_invalid & idx_hasyear] = fixed
  
  # set remaining to NA
  idx_invalid = (!grepl("^(19|20)\\d{2}$", data$study.year, perl = TRUE))
  data$study.year[idx_invalid] = NA
  data$study.year = as.integer(data$study.year)
  
  return(data)
}

# Reads the data from meta-analysis
pb.readMA = function(path, file) {
  
  data = read.csv(file.path(path, file),
                  colClasses=c("factor",                             # id
                               "character", "character", "factor",   # outcome id, name, flag
                               "logical",  "logical",                # hasSubgroups, isSubgroup
                               "numeric", "numeric", "numeric",      # effect, ci_start, ci_end
                               "numeric", "numeric", "factor",       # z, p_z, estimable
                               "numeric", "numeric", "numeric",      # studies, hetero.chi2, hetero.df
                               "numeric", "numeric", "numeric",      # hetero.p_chi2, I2, tau2
                               "numeric", "numeric"                  # total1, total2
                  ))
  
  return(data)
}

# Cleans the raw data with some reg expressions
pb.merge.outcome.measures = function(data) {
  
  names = unique(data$outcome.measure)
  data$outcome.measure.merged = data$outcome.measure
  
  # Define aliases, first element will be used for all aliases
  # https://handbook-5-1.cochrane.org/chapter_9/9_2_2_5_what_is_the_event.htm
  aliases=list()
  aliases[[1]] = grep("risk ratio|^RR$|RR Ratio[s]*|IRR|relat[ive]* risk", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[2]] = grep("^mean dif[ference]*$|^mean dif[ference]*\\s+.[^%]|^MD", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[3]] = grep("change in", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[4]] = grep("Std. mean|SMD", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[5]] = grep("hedges", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[6]] = sort(grep("% change|% increase|% rate|risk difference \\(%\\)|mean difference \\[%\\]|changes in MVC \\[%]|changes \\[%]",
                           names, perl = TRUE, ignore.case = TRUE, value = TRUE))
  aliases[[7]] = grep("^OR$|^odds ratio[s]*|[^peto] odds ratio", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[8]] = grep("peto", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[9]] = grep("hazard|^HR$|Survival HR|HR and variance", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[10]] = grep("rate ratio", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[11]] = grep("risk dif[ference]*$|risk dif[ference]*\\s+.[^%]", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[12]] = grep("^rate dif", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[13]] = grep("prevented fraction", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  
  for (i in 1:length(aliases)) {
    idx =  data$outcome.measure.merged %in% aliases[[i]] # find all aliases
    data$outcome.measure.merged[idx] = aliases[[i]][1] # replace consistently with 1st element
  }
  
  return(list(data, aliases))
}

# Search keyword in database and returns all rows that match.
pb.search = function(keyword, data) {
  if ("rev.title" %in% names(data)) {
    idx1 = grepl(keyword, data$rev.title, ignore.case = T)
    return(data[idx1,c("file.nr", "doi", "rev.title", 
                       "rev.year","nr.studies","pool.count1","pool.count2")])
    
  } else {
    idx1 = grepl(keyword, data$outcome.name, ignore.case = T)
    idx2 = grepl(keyword, data$comparison.name, ignore.case = T)
    idx3 = grepl(keyword, data$study.name, ignore.case = T)
    idx4 = grepl(keyword, data$subgroup.name, ignore.case = T)
    return(data[idx1|idx2|idx3|idx4,c("file.nr","comparison.name","outcome.name","outcome.measure", 
                                      "subgroup.name","study.name","total1","total2")])
  }
  
}

# cleans journals names by removing some special characters and caring about
# common inconsistencies for a better matching of journal titles.
pb.cleanJournalNames = function(titles) {
  titles = tolower(titles)
  titles = gsub(",|-|‐|–|:|\\.|\\sand of\\s|\\sof the\\s|\\sand the\\s|^the\\s|\\sthe\\s|\\sof\\s|&\\s|\\sand\\s|\\(.*\\)|\\[.*\\]|(19|20)\\d{2}|\\/", " ", titles) # remove stuff
  titles = gsub("jounal|jounral", "journal", titles) # spelling
  titles = gsub("medicial", "medical", titles) # spelling
  titles = gsub("\\s{2,}", " ", titles) # replace double spaces
  titles = gsub("^\\s|\\s$", "", titles) # replace leading/ending spaces
  
  return(titles)
}