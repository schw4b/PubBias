# by Simon Schwab, 2019

library(testit)

# Reads the data
pb.readData = function(path, file) {
  
  data = read.csv(file.path(path, file),
                  colClasses=c("integer", "character", "character",  # file       :  nr, name, doi
                               "integer", "integer",                 # file       :  index, version
                               "integer", "character",               # comparison :  nr, name
                               "integer", "character", "character",  # outcome    :  nr, name, measure
                               "integer", "character",               # subgroup   :  nr, name
                               "character", "integer",               # study      :  name, year
                               "numeric", "numeric",                 # effect, std.err.
                               "numeric", "numeric",                 # group 1    :  events, total
                               "numeric", "numeric",                 # group 1    :  mean, se
                               "numeric", "numeric",                 # group 2    :  events, total
                               "numeric", "numeric",                 # group 2    :  mean, se
                               "numeric"))                           # total N
  # some cleanup
  data$N = NULL # this colums has all NA's
  data$outcome.measure[grep("&#223;-agonist", data$outcome.measure)] = "beta-agonist"
  data$outcome.measure[grep("Hedge[s ]*&#180", data$outcome.measure)] = "Hedges' g"
  
  return(data)
}

# Reads the data
pb.readData2 = function(path, file) {
  
  data = read.csv(file.path(path, file),
                  colClasses=c("character",                          # id
                               "integer", "character", "character",  # comparison :  nr, name, id
                               "integer", "character", "character",  # outcome    :  nr, name, measure
                               "character", "character",             #               id, flag
                               "integer",  "character", "character", # subgroup   :  nr, name, id
                               "character","character","character",  # study      :  id, name, year
                               "character",                          #               data_source  
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
  data$study.year = as.numeric(data$study.year)
  
  return(data)
}


# Cleans the raw data with some reg expressions
pb.clean = function(data) {
  
  names = unique(data$outcome.measure)
  data$outcome.measure.new = data$outcome.measure
  
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
  aliases[[7]] = grep("^odds ratio[s]*|[^peto] odds ratio", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[8]] = grep("peto", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[9]] = grep("hazard|^HR$|Survival HR|HR and variance", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[10]] = grep("rate ratio", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[11]] = grep("risk dif[ference]*$|risk dif[ference]*\\s+.[^%]", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[12]] = grep("^rate dif", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  aliases[[13]] = grep("prevented fraction", names, perl = TRUE, ignore.case = TRUE, value = TRUE)
  
  
  
  for (i in 1:length(aliases)) {
    idx =  data$outcome.measure.new %in% aliases[[i]] # find all aliases
    data$outcome.measure.new[idx] = aliases[[i]][1] # replace consistently with 1st element
  }
  
  return(list(data, aliases))
}

# Adds pool.nr. We need to know what outcomes and subgroups can be combined in a pooled analysis.
# Searches for duplicates in {file.nr, comparison.nr, outcome.nr, subgroup.nr}.
pb.pool = function(data) {
  
  # create table to check for duplicates, these get the same pool nr
  tab = cbind(data$file.nr, data$comparison.nr, data$outcome.nr, data$subgroup.nr)
  idx = duplicated(tab)
  
  data$pool.nr = rep(NA, nrow(data))
  from = which(!idx)
  to   = c(from[2:length(from)] - 1, length(idx))
  
  # iterate through list of duplicates and assign nr
  c = 1
  for (i in 1:length(from)) {
    data$pool.nr[from[i]:to[i]] = c
    c = c + 1
  }
  return(data)
}

## Creates a database of reviews (each line one review).
## Fetches the title, year for each review
## Adds additional variables
#library(roadoi)
pb.createReviews = function(data) {
  file_fetch = 'oadoi_fetch.RData'
  table = data.frame(file.nr = data$file.nr, doi=data$doi,
                     rev.title = NA, rev.year = NA)
  table$doi = as.character(table$doi)
  table = table[!duplicated(table),]
  
  # do not fetch if file already exists, takes ~30 min to get 5,000 titles
  if (file.exists(file.path(PATH_RESULTS, file_fetch))) {
    load(file.path(PATH_RESULTS, file_fetch))
  } else {
    fetch = oadoi_fetch(dois = table$doi, email = "simon.schwab@uzh.ch")
    save(fetch, file = file.path(PATH_RESULTS, file.fetch))
  }
  
  # merge fetched data with database
  idx = match(tolower(fetch$doi), tolower(table$doi))
  table$rev.title = rep(NA, nrow(table))
  table$rev.year = rep(NA, nrow(table))
  table$rev.title[idx] = fetch$title
  table$rev.year[idx] = fetch$year
  
  # populate review database with important variables
  
  # add number of studies per review and pool variables
  table$nr.studies = rep(NA, nrow(table))
  table$pool.nr1 = rep(NA, nrow(table)) 
  table$pool.outName1 = rep(NA, nrow(table))
  table$pool.count1 = rep(NA, nrow(table)) # nr of outcomes that can be pooled
  table$pool.nr2 = rep(NA, nrow(table))
  table$pool.outName2 = rep(NA, nrow(table)) 
  table$pool.count2 = rep(NA, nrow(table))
  
  for (i in 1:nrow(table)) {
    d = subset(data, data$file.nr == table$file.nr[i])
    table$nr.studies[i] = length(unique(d$study.name))
    s = sort(table(d$pool.nr), decreasing = TRUE) # count and sort pool.nr
    
    table$pool.nr1[i] = names(s[1])
    table$pool.nr2[i] = names(s[2])
    
    table$pool.count1[i] = s[1]
    table$pool.count2[i] = s[2]
    
    # get autcome name
    on = d$outcome.name[d$pool.nr == as.numeric(names(s[1]))]
    assert(length(unique(on)) == 1) # check all the same outcomes
    table$pool.outName1[i] = on[1]
    
    on = d$outcome.name[d$pool.nr == as.numeric(names(s[2]))]
    assert(length(unique(on)) == 1) # check all the same outcomes
    table$pool.outName2[i] = on[2]
  }
  
  return(table)
}

## Search keyword in database and returns all rows that match.
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

library(rvest)
pb.crawlCochrane = function(dois) {
  
  mylist = list() # stores a list of data frames per review with all the studies
  
  # iterate through reviews
  for (m in 1:length(dois)) {
    
    print(paste("review nr:", m))
    url = paste0("https://www.cochranelibrary.com/cdsr/doi/", dois[m], "/references")
    src = read_html(url)
    #if (m %% 8 == 0) { Sys.sleep(40) } # or we get blacklisted for a while
    Sys.sleep(7)
    studies = html_nodes(css="div.bibliographies.references_includedStudies", x=src)
    n = length(studies)
    
    studies.df = data.frame(refShort = rep(NA, n), authors = rep(NA, n), title = rep(NA, n),
                            journal = rep(NA, n), year = rep(NA, n), doi = rep(NA, n),
                            timesCited = rep(NA, n), refLong = rep(NA, n), error = rep(NA, n))
    
    for (i in 1:n) { # iterate through studies of a meta-analysis
      
      #print(i)
      refShort = html_text(html_nodes(css="h4.title", x=studies[i]))
      refShort = sub("\\s[\\{\\+].*", "", refShort)
      # often multiple references are listed for each study called records
      records = html_nodes(css="div.bibliography-section", x=studies[i])
      
      idx_rec = 1 # take fist record if the is only one
      
      # if there are more than one
      if (length(records) > 1) {
        # select correct record and remove all special characters like Malmström -> Malmstr.m and 2018a -> 2018
        idx_rec = grep(paste0("^\\s*[the]*\\s*", gsub("[^A-Za-z0-9]+|[a-z]$", ".*", refShort)),
                       html_text(records), ignore.case = TRUE)
      }
      
      # if none has been found
      if (length(idx_rec) == 0) {  # relax on year only (due Country codes like Chile, Europe as First Author name)
        idx_rec = grep(sub("(.*)(19[0-9]{2}|20[0-9]{2})(.*)", "\\2", refShort, perl = TRUE), # match year 19** or 20**
                       html_text(records), ignore.case = TRUE)
      }
      
      # if there are still more than 1...
      if (length(idx_rec) > 1) { # if more than one record, take the one that is cited more often, or has more links to original work
        moreCites = rep(NA, length(idx_rec))
        moreLinks = rep(NA, length(idx_rec))
        for (j in 1:length(idx_rec)) { # for each duplicate with same shortRef author/year get citations and links
          links = html_nodes("a", x = records[idx_rec[j]])
          moreLinks[j] = length(links)
          idx = grep("Web of Science", html_text(links), ignore.case = TRUE) # determine link times cited
          timesCited = s(html_text(html_nodes(css="span", x=links[idx])))
          moreCites[j] = as.numeric(sub(".*:", "", timesCited))
        }
        if (all(is.na(moreCites))) {
          idx_rec = idx_rec[which.max(moreLinks)]
        } else {
          idx_rec = idx_rec[which.max(moreCites)]
        }
      }
      
      # if idx_rec is still not 1 then something went wrong.
      if (length(idx_rec) != 1) {
        authors = title = journal = year = doi = timesCited = refLong = NA
        error = 1
        
      } else { # success: fetch and fill data
        
        refLong  = html_text(html_nodes(css="div", x=records[idx_rec]))
        authors  = sub("\\..*", "", refLong)
        title    = s(html_text(html_nodes(css="span.citation-title", x=records[idx_rec])))
        journal  = gsub("^\\s|\\s$", "", s(html_text(html_nodes(css="span.citation", x=records[idx_rec]))))
        year     = s(html_text(html_nodes(css="span.pubYear", x=records[idx_rec])))
        
        links    = html_nodes("a", x = records[idx_rec]) # select
        idx = grep("Link to article", html_text(links), ignore.case = TRUE) # determine link doi
        doi = s(html_attr(name="href", x=links[idx]))
        
        idx = grep("Web of Science", html_text(links), ignore.case = TRUE) # determine link times cited
        timesCited = s(html_text(html_nodes(css="span", x=links[idx])))
        timesCited = sub(".*:", "", timesCited)
        error = 0
        
      }
      # n x 9 data frame
      studies.df[i,] = c(refShort, authors, title, journal, year, doi, timesCited, refLong, error)
    }
    
    mylist[[m]] = studies.df
  }
  
  return(mylist)
}

require(tidyverse)
require(meta)

pb.bias.bin <- function(data){
  metadat <- data %>% filter(outcome.measure == "Risk Ratio" | outcome.measure == "Odds Ratio") %>% 
    filter(file.nr != 3014) %>% #file.nr 3014 gibt eine seltsame Fehlermeldung, muss ausgeschlossen werden
    filter(events1 > 0 | events2 > 0) %>% #comparisons mit events  = 0 ausschliessen
    filter(total1 - events1 > 0 | total2 - events2 > 0) %>% #comparisons mit events = total ausschliessen - vertr??gt der alg. irgendwie nicht.
    group_by(file.nr, outcome.nr, subgroup.nr) %>% #Die nachfolgenden Rechnungen werden jeweils f??r diese Gruppen gemacht.
    mutate(n = n()) %>% filter(n > 9) %>% #Nur meta-analysen mit mehr als 9 comparisons behalten
    summarize(doi = unique(doi), n = n(), #Summarize fasst alle daten in einer Gruppe zusammen
              pval.peters.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "peters")$p.val,
              pval.harbord.bin = metabias(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR"), method = "score")$p.val,
              trim.bin = trimfill(metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2))$k0 / n(),
              Q.bin = metabin(event.e = events1, n.e = total1, event.c = events2, n.c = total2, sm = "OR")$Q)
  return(metadat)
}     

pb.bias.cont <- function(data){
  metadat <- data %>% filter(outcome.measure == "Mean Difference" | outcome.measure == "Std. Mean Difference") %>%
    filter(sd1 > 0 & sd2 > 0 ) %>% filter(!is.na(sd1) & !is.na(sd2)) %>% #Nur comparisons mit sinnvollen sd's behalten
    filter(mean1 != 0 | mean2 != 0 ) %>% filter(!is.na(mean1) & !is.na(mean2)) %>% #Nur comparisons behalten die nicht beide means  = 0 haben
    group_by(file.nr, outcome.nr, subgroup.nr) %>% #Siehe oben
    mutate(n = n()) %>% filter(n > 9) %>% 
    summarize(doi = unique(doi), n = n(),
              pval.egger.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                                         method = "linreg")$p.val,
              pval.thomsom.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                                           method = "mm")$p.val,
              pval.begg.cont = metabias(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2), 
                                        method = "rank")$p.val,
              trim.cont = trimfill(metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2))$k0 / n(),
              Q.cont = metacont(n.e = total1, mean.e = mean1, sd.e = sd1, n.c = total2, mean.c = mean2, sd.c = sd2)$Q) 
  return(metadat)
}

# Returns NA if string has length 0
s <- function(object) {
  if (length(object) == 0) {
    return(NA)
  } else {
    return(object)
  }
}
