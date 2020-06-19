## by Simon Schwab 2019, 2020

library(rvest)
library(httr)
# get cochrane library data
get.cola = function(doi, path) {
  
  file_log = file.path(path, "cola.log")
  
  for (i in 1:length(doi)) {
    id = strsplit(doi[i], split = "\\.")[[1]][3]
    file_url = paste0(id, "StatsDataOnly.rm5")
    file_out = paste0(id, ".xml")
    
    # if file is not yet in path, get it.
    if (!file.exists(file.path(path, file_out))) {
      url = paste0("https://www.cochranelibrary.com/cdsr/doi/", doi[i])
      url_data = paste0(url, "/media/CDSR/", id, "/table_n/", file_url)
      
      # create www session
      #ua = user_agent("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36")
      ua = user_agent("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_14_5) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/79.0.3945.130 Safari/537.36")
      session_with_ua = html_session(url, ua)
      if (!session_with_ua$response$status_code == 200) {
        cat(paste(Sys.time(), "error: url not found", doi[i], "\n"), file = file_log, append = TRUE)
      } else {
        session_data = jump_to(session_with_ua, url_data) # jump to dataset
        
        # write to file if session is ok
        if (session_data$response$status_code == 200) {
          bin = content(session_data$response, as = "raw")
          writeBin(bin, file.path(path, file_out))
          cat(paste(Sys.time(), "ok: download successful", doi[i], "\n"), file = file_log, append = TRUE)
        } else {
          cat(paste(Sys.time(), "error: data not found", doi[i], "\n"), file = file_log, append = TRUE)
        }
      }
      Sys.sleep(2)
    } else {
      cat(paste(Sys.time(), "ok: file already exists", doi[i], "\n"), file = file_log, append = TRUE)
      Sys.sleep(0.1) # this is needed or we get blacklisted
    }
  }
  
}

parse.studies = function(file, path) {
  # Structure of the Cochrane library
  # COMPARISON - NAME  
  #            - IV_/CONT_|DICH_OUTCOME (k) - NAME
  #                                         - CONT_DATA (k)
  #                                         - CONT_SUBGROUP - NAME (1)
  #                                                         - CONT_DATA (k)
  
  # Examples
  # empty: file="CD005507.xml"
  # html:  file="CD005602.xml"
  # xml:   file="CD005601.xml"
  
  # check file exits
  if (!file.exists(file.path(path, file))) {
    stop(paste("File not found", file))
  }
  
  # check file is empty
  info = file.info(file.path(path, file))
  if (info$size == 0) {
    stop(paste("Empty file", file))
  }
  
  # check file is XML
  con = file(file.path(path, file), "r")
  first_line = readLines(con, n=1)
  close(con)
  if (!grepl("xml version", substr(first_line,1,20))) {
    stop(paste("File is not XML", file))
  }
  
  xmlsrc = read_xml(file.path(path, file))
  
  # get total number of outcomes and create data frame
  l = length(html_nodes(css=c("DICH_DATA,CONT_DATA,IV_DATA,IPD_DATA"), x=xmlsrc))
  
  if (l == 0) { # in case there is no outcome data in the file
    stop(paste("No outcome data found for", file))
  }
  
  tab = data.frame(id=rep(NA, l), comparison.nr = rep(NA, l), comparison.name = rep(NA, l), comparison.id = rep(NA, l),
                   outcome.nr = rep(NA, l), outcome.name = rep(NA, l), outcome.measure = rep(NA, l), outcome.id = rep(NA, l), outcome.flag = rep(NA, l),
                   subgroup.nr = rep(NA, l), subgroup.name = rep(NA, l), subgroup.id = rep(NA, l), 
                   study.id = rep(NA, l), study.name = rep(NA, l), study.year = rep(NA, l), study.data_source = rep(NA, l),
                   effect = rep(NA, l), se = rep(NA, l),
                   events1 = rep(NA, l), total1 = rep(NA, l), mean1 = rep(NA, l), sd1 = rep(NA, l),
                   events2 = rep(NA, l), total2 = rep(NA, l), mean2 = rep(NA, l), sd2 = rep(NA, l))
  
  id = sub('\\.xml$', '', file)
  
  # 1. Table with study information
  xml.study = html_nodes(css=c("STUDIES STUDY"), x=xmlsrc)
  j = length(xml.study)
  table.studies = data.frame(id=rep(NA,j), name=rep(NA,j), year=rep(NA,j), data_source=rep(NA,j))
  table.studies$id = html_attr(x = xml.study, name="ID")
  table.studies$name = html_attr(x = xml.study, name="NAME")
  table.studies$year = gsub("\\.", "", html_attr(x = xml.study, name="YEAR")) # 2006. -> 2006
  table.studies$data_source = html_attr(x = xml.study, name="DATA_SOURCE")
  table.studies$id = gsub("--", "-", table.studies$id) # fix study id
  table.studies$name = gsub("  ", " ", table.studies$name) # fix study id
  
  # 2. Parse single study data from comparisons, outcomes, subgroups, and data nodes
  #    We have to loop across comparison, outcome and subgroups as the data from single
  #    studies have no ID but are hierarchically structured.
  k = 1 # counter
  xml.comparison = html_nodes(css=c("COMPARISON"), x=xmlsrc)
  for (c in 1:length(xml.comparison)) {
    # print(paste("Comparison:", c))
    comparison.nr = html_attr(x = xml.comparison[c], name="NO")
    comparison.name = cleanstr(html_text(html_node(css=c("NAME"), x=xml.comparison[c])))
    comparison.id = html_attr(x = xml.comparison[c], name="ID")
    
    xml.outcome = html_nodes(css=c("DICH_OUTCOME,CONT_OUTCOME,IV_OUTCOME,IPD_OUTCOME"), x=xml.comparison[c])
    
    # there are five type of outcomes, see
    # 
    # CONT_OUTCOME - continuous data with mean and SD
    # DICH_OUTCOME - dichotomous data with no. of events
    # IV_OUTCOME   - inverse variance meta-analysis in case summary data for groups
    #                is unavailable. Inludes effects and SE only.
    # IPD_OUTCOME  - individual patient data
    # OTHER_OUTCOME - Additional informationn, notes, etc.
    if (length(xml.outcome) > 0) { # handle if xml node is empty
      for (o in 1:length(xml.outcome)) {
        # print(paste("Outcome:", o))
        outcome.nr = html_attr(x = xml.outcome[o], name="NO")
        outcome.name = cleanstr(html_text(html_node(css=c("NAME"), x=xml.outcome[o])))
        outcome.measure = html_attr(x = xml.outcome[o], name="EFFECT_MEASURE")
        outcome.id = html_attr(x = xml.outcome[o], name="ID")
        outcome.flag = NA
        if (grepl("^<CONT_OUTCOME", xml.outcome[o]))  { outcome.flag = "CONT" }
        if (grepl("^<DICH_OUTCOME", xml.outcome[o]))  { outcome.flag = "DICH" }
        if (grepl("^<IV_OUTCOME", xml.outcome[o]))    { outcome.flag = "IV" }
        if (grepl("^<IPD_OUTCOME", xml.outcome[o]))   { outcome.flag = "IPD" }
        if (grepl("^<OTHER_OUTCOME", xml.outcome[o])) { outcome.flag = "OTHER" }
        
        if (is.na(outcome.measure)) { # DICH and CONT have a attribute, IV and IPD have a node for effect measure.
          outcome.measure = html_text(html_nodes(css=c("EFFECT_MEASURE"), x=xml.outcome[o]))
        }
        
        has.subgroup = FALSE
        xml.subgroup = 1
        subgroup.nr = 0
        subgroup.name = ""
        subgroup.id = NA
        if (html_attr(name = "SUBGROUPS", x = xml.outcome[o]) == "YES") {
          has.subgroup = TRUE
          xml.subgroup = html_nodes(css=c("DICH_SUBGROUP,CONT_SUBGROUP,IV_SUBGROUP,IPD_SUBGROUP"), x=xml.outcome[o])
        }
        # print(paste("Subgroup:", has.subgroup))
        for (s in 1:length(xml.subgroup)) {
          # print(paste("Subgroup:", s))
          if (has.subgroup) {
            subgroup.nr = html_attr(x = xml.subgroup[s], name="NO")
            subgroup.name = cleanstr(html_text(html_node(css=c("NAME"), x=xml.subgroup[s])))
            xml.data = html_nodes(css=c("DICH_DATA,CONT_DATA,IV_DATA,IPD_DATA"), x=xml.subgroup[s])
            subgroup.id = html_attr(x = xml.subgroup[s], name="ID")
          } else {
            xml.data = html_nodes(css=c("DICH_DATA,CONT_DATA,IV_DATA,IPD_DATA"), x=xml.outcome[o])
          }
          
          if (length(xml.data) > 0) {
            for (d in 1:length(xml.data)) {
              # print(paste("Data:", d, "k:", k))
              tab$id[k] = sub('\\.xml$', '', file)
              
              tab$comparison.nr[k] = comparison.nr
              tab$comparison.name[k] = comparison.name
              tab$comparison.id[k] = comparison.id
              
              tab$outcome.nr[k] = outcome.nr
              tab$outcome.name[k] = outcome.name
              tab$outcome.measure[k] = outcome.measure
              tab$outcome.id[k] = outcome.id
              tab$outcome.flag[k] = outcome.flag
              
              tab$subgroup.nr[k] = subgroup.nr
              tab$subgroup.name[k] = subgroup.name
              tab$subgroup.id[k] = subgroup.id
              
              tab$study.id[k] = gsub("--", "-", html_attr(x = xml.data[d], name="STUDY_ID")) # replce "--" with "-" in study id
              idx = which(table.studies$id %in% tab$study.id[k])[1]
              tab$study.name[k] = table.studies$name[idx]
              tab$study.year[k] = table.studies$year[idx]
              tab$study.data_source[k] = table.studies$data_source[idx]
              
              tab$effect[k]  = as.numeric(html_attr(x = xml.data[d], name="EFFECT_SIZE"))
              tab$se[k]      = as.numeric(html_attr(x = xml.data[d], name="SE"))
              tab$events1[k] = as.numeric(html_attr(x = xml.data[d], name="EVENTS_1"))
              tab$total1[k]  = as.numeric(html_attr(x = xml.data[d], name="TOTAL_1"))
              tab$mean1[k]   = as.numeric(html_attr(x = xml.data[d], name="MEAN_1"))
              tab$sd1[k]     = as.numeric(html_attr(x = xml.data[d], name="SD_1"))
              tab$events2[k] = as.numeric(html_attr(x = xml.data[d], name="EVENTS_2"))
              tab$total2[k]  = as.numeric(html_attr(x = xml.data[d], name="TOTAL_2"))
              tab$mean2[k]   = as.numeric(html_attr(x = xml.data[d], name="MEAN_2"))
              tab$sd2[k]     = as.numeric(html_attr(x = xml.data[d], name="SD_2"))
              
              k = k + 1
            } # data loop
          }
        } # subgroup loop
      } # outcome loop
    }
  } # comparison loop
  return(tab)
}

parse.ma = function(file, path) {
  
  # check file exits
  if (!file.exists(file.path(path, file))) {
    stop(paste("File not found", file))
  }
  
  # check file is empty
  info = file.info(file.path(path, file))
  if (info$size == 0) {
    stop(paste("Empty file", file))
  }
  
  # check file is XML
  con = file(file.path(path, file), "r")
  first_line = readLines(con, n=1)
  close(con)
  if (!grepl("xml version", substr(first_line,1,20))) {
    stop(paste("File is not XML", file))
  }
  
  xmlsrc = read_xml(file.path(path, file))
  
  # get total number of outcomes and create data frame
  xml.overall = html_nodes(css=c("DICH_OUTCOME,CONT_OUTCOME,IV_OUTCOME,IPD_OUTCOME,
                                  DICH_SUBGROUP,CONT_SUBGROUP,IV_SUBGROUP,IPD_SUBGROUP"), x=xmlsrc)
  
  l = length(xml.overall)
  
  if (l == 0) { # in case there is no outcome data in the file
    stop(paste("No outcome data found for", file))
  }
  
  tab = data.frame(id=rep(NA, l), 
                   outcome.id = rep(NA, l), name = rep(NA, l), outcome.flag = rep(NA, l),
                   hasSubgroups = rep(NA, l), isSubgroup = rep(FALSE, l),
                   effect = rep(NA, l), ci_start = rep(NA, l), ci_end = rep(NA, l),
                   z = rep(NA, l), p_z = rep(NA, l), estimable = rep(NA, l), studies = rep(NA, l),
                   hetero.chi2 = rep(NA, l), hetero.df = rep(NA, l), hetero.p_chi2 = rep(NA, l),
                   I2 = rep(NA, l), tau2 = rep(NA, l), # tau is missing in XML, e.g. CD012282 CMP-004.01 
                   total1 = rep(NA, l), total2 = rep(NA, l))
  
  for (o in 1:length(xml.overall)) {
    
    tab$id[o] = sub('\\.xml$', '', file)
    tab$outcome.id[o] = html_attr(x = xml.overall[o], name="ID") # both outcome and subgroup id
    tab$name[o] = cleanstr(html_text(html_node(css=c("NAME"), x=xml.overall[o]))) # both outcome and subgroup names
    tab$outcome.flag[o] = NA
    if (grepl("^<CONT_(OUTCOME|SUBGROUP)", xml.overall[o]))  { tab$outcome.flag[o] = "CONT" }
    if (grepl("^<DICH_(OUTCOME|SUBGROUP)", xml.overall[o]))  { tab$outcome.flag[o] = "DICH" }
    if (grepl("^<IV_(OUTCOME|SUBGROUP)", xml.overall[o]))    { tab$outcome.flag[o] = "IV" }
    if (grepl("^<IPD_(OUTCOME|SUBGROUP)", xml.overall[o]))   { tab$outcome.flag[o] = "IPD" }
    if (grepl("^<OTHER_(OUTCOME|SUBGROUP)", xml.overall[o])) { tab$outcome.flag[o] = "OTHER" }
    
    tab$hasSubgroups[o] = html_attr(x = xml.overall[o], name="SUBGROUPS") == "YES"
    if (grepl("^<(DICH|CONT|IV|IPD)_SUBGROUP", xml.overall[o])) { 
      tab$isSubgroup[o] = TRUE
    }
    
    tab$effect[o] = as.numeric(html_attr(x = xml.overall[o], name="EFFECT_SIZE"))
    tab$ci_start[o] = as.numeric(html_attr(x = xml.overall[o], name="CI_START"))
    tab$ci_end[o] = as.numeric(html_attr(x = xml.overall[o], name="CI_END"))
    tab$z[o] = as.numeric(html_attr(x = xml.overall[o], name="Z"))
    tab$p_z[o] = as.numeric(html_attr(x = xml.overall[o], name="P_Z"))
    tab$estimable[o] = html_attr(x = xml.overall[o], name="ESTIMABLE")
    
    tab$studies[o] = as.numeric(html_attr(x = xml.overall[o], name="STUDIES"))
    tab$hetero.chi2[o] = as.numeric(html_attr(x = xml.overall[o], name="CHI2"))
    tab$hetero.df[o] = as.numeric(html_attr(x = xml.overall[o], name="DF"))
    tab$hetero.p_chi2[o] = as.numeric(html_attr(x = xml.overall[o], name="P_CHI2"))
    
    tab$I2[o] = as.numeric(html_attr(x = xml.overall[o], name="I2"))
    tab$tau2[o] = as.numeric(html_attr(x = xml.overall[o], name="TAU2"))
    
    tab$total1[o] = as.numeric(html_attr(x = xml.overall[o], name="TOTAL_1"))
    tab$total2[o] = as.numeric(html_attr(x = xml.overall[o], name="TOTAL_2"))
    
  }
  
  # cleaning
  tab$hasSubgroups[is.na(tab$hasSubgroups)] = FALSE
  
  return(tab)
}

cleanstr = function(string) {
  return(gsub("\n|\"|\\\\|[[:space:]]*$", "", string))
}

# asserts to vectors of numericals are same
# in diagnose mode, returns unequal pairs with corresponding row index
assertnum = function(x, y, p = 10^-7, diagnose=FALSE) {
  rnames = 1:length(x)
  idx = !is.na(x) & !is.na(y)
  x_=x[idx]
  y_=y[idx]
  rnames_ = rnames[idx]
  if(diagnose) {
    uneq = abs(x_ - y_) > p
    tab = cbind(x_[uneq], y_[uneq])
    rownames(tab) = rnames_[uneq]
    return( tab )
    
  } else {
    assert( abs(x_ - y_) < p )
  }
}

parseReferences = function(cddoi) {
  
  url = paste0("https://www.cochranelibrary.com/cdsr/doi/", cddoi, "/references")
  src = read_html(url)
  id = strsplit(cddoi, split = "\\.")[[1]][3]
  # check if biblipgraphy is empty
  if (length(html_nodes(css="div.bibliographies", src)) == 0)
    return(1)
  else {
    
    studies = html_nodes(css="div.bibliographies.references_includedStudies", x=src)
    characteristics = html_nodes(css="section.characteristicIncludedStudiesContent", x=src)
    tables = html_nodes(css="div.table", x=characteristics)
    
    # sanity checks
    assert(length(studies) == length(tables), fact = "Study references and table characteristics mismatch.")
    check.studies = trimws(gsub("\\{.*", "", word(html_text(html_nodes(css="h4.title", x=studies)), 1, 2)))
    check.tables = trimws(gsub("\\{.*", "", word(html_text(html_nodes(css="span.table-title", x=tables)), 1, 2)))
    idx.studies = order(check.studies)
    idx.tables  = order(check.tables)
    # References and tables mismatch may mismatch in few of the reviews
    # cbind(check.studies, check.tables, check.studies == check.tables, idx.studies, idx.tables)
    assert(all(check.studies[idx.studies] == check.tables[idx.tables]), fact = "Study references and table characteristics order mismatch.")
    studies = studies[idx.studies]
    tables = tables[idx.tables]
    
    n = length(studies)
    
    df.ref = data.frame(id = rep(NA, n), refShort = rep(NA, n), authors = rep(NA, n), title = rep(NA, n),
                        journal = rep(NA, n), year = rep(NA, n), pmid = rep(NA, n),
                        gschol = rep(NA, n), refLong = rep(NA, n))
    df.char = data.frame(char.methods = rep(NA, n), char.participants = rep(NA, n), char.interventions = rep(NA, n),
                         char.outcomes = rep(NA, n), char.notes = rep(NA, n))
    df.bias = data.frame(bias.random_sequence = rep(NA, n), bias.allocation_concealment = rep(NA, n), 
                         bias.blinding_participants = rep(NA, n), bias.blinding_outcomes = rep(NA, n),
                         bias.incomplete_outcome = rep(NA, n), bias.selective_reporting= rep(NA, n))
    
    for (i in 1:n) { # iterate through studies of a review and get reference and characteristics
      
      #print(i)
      refShort = html_text(html_nodes(css="h4.title", x=studies[i]))
      refShort = sub("\\s[\\{\\+].*", "", refShort)
      # often multiple references are listed for each study called records
      records = html_nodes(css="div.bibliography-section", x=studies[i])
      
      idx_rec = 1 # take fist record if the is only one
      
      # if there are more than one record
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
      
      # if there are still more than 1 take the ref with most links
      if (length(idx_rec) > 1) { # if more than one record, take the one that is cited more often, or has more links to original work
        moreLinks = rep(NA, length(idx_rec))
        for (j in 1:length(idx_rec)) { # for each duplicate with same shortRef author/year get links
          links = html_nodes(".citation-link", x = records[idx_rec[j]])
          moreLinks[j] = length(links)
        }
        idx_rec = idx_rec[which.max(moreLinks)]
      }
      
      # if idx_rec is still not 1 then something went wrong.
      if (length(idx_rec) != 1) {
        refLong = authors = title = journal = year = pmid = gschol =  NA
        
      } else { # success: fetch and fill data
        
        refLong  = s(html_text(html_nodes(css="div", x=records[idx_rec])))
        authors  = sub("\\..*", "", refLong)
        title    = s(html_text(html_nodes(css="span.citation-title", x=records[idx_rec])))
        journal  = gsub("^\\s|\\s$", "", s(html_text(html_nodes(css="span.citation", x=records[idx_rec]))))
        year     = s(html_text(html_nodes(css="span.pubYear", x=records[idx_rec])))
        
        # Get links to PubMed and Google Scholar
        links    = html_nodes("a", x = records[idx_rec]) # select
        idx = grep("PubMed", html_text(links), ignore.case = TRUE) # Pubmed
        pmid = s(html_attr(name="href", x=links[idx]))
        idx = grep("Google Scholar", html_text(links), ignore.case = TRUE) # Google Scholar
        gschol = s(html_attr(name="href", x=links[idx]))
        
      }
      # n x 8 data frame
      df.ref[i,] = c(id, refShort, authors, title, journal, year, pmid, gschol, refLong)
      
      # add study characteristics
      rows = html_nodes("tr", x = tables[i]) 
      for (k in 1:length(rows)) { # iterate through rows of study characteristic table from study i
        td = html_nodes("td", x = rows[k]) # extract table data
        
        if ( grepl("^Methods$", trimws(html_text(td[1]))) ) { df.char$char.methods[i] = html_text(td[2])}
        if ( grepl("^Participants$", trimws(html_text(td[1]))) ) { df.char$char.participants[i] = html_text(td[2])}
        if ( grepl("^Interventions$", trimws(html_text(td[1]))) ) { df.char$char.interventions[i] = html_text(td[2])}
        if ( grepl("^Outcomes$", trimws(html_text(td[1]))) ) { df.char$char.outcomes[i] = html_text(td[2])}
        if ( grepl("^Notes$", trimws(html_text(td[1]))) ) { df.char$char.notes[i] = html_text(td[2])}
        
        if ( grepl("Random sequence generation", html_text(td[1])) ) { df.bias$bias.random_sequence[i] = html_text(td[2])}
        if ( grepl("Allocation concealment", html_text(td[1])) ) { df.bias$bias.allocation_concealment[i] = html_text(td[2])}
        if ( grepl("Blinding of participants", html_text(td[1])) ) { df.bias$bias.blinding_participants[i] = html_text(td[2])}
        if ( grepl("Blinding of outcome", html_text(td[1])) ) { df.bias$bias.blinding_outcomes[i] = html_text(td[2])}
        if ( grepl("Incomplete outcome", html_text(td[1])) ) { df.bias$bias.incomplete_outcome[i] = html_text(td[2])}
        if ( grepl("Selective reporting", html_text(td[1])) ) { df.bias$bias.selective_reporting[i] = html_text(td[2])}
        
      }
      
    } # close for each study
    return(cbind(df.ref, df.char, df.bias))
  }
}