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


# transforms logORs so that we get a estimate and SE for all the levels (and not k-1 levels).
translogOR = function(theta, SEtheta, CovTheta) {
  
  n = length(theta) + 1
  theta = c(0, theta)
  SEtheta = c(0, SEtheta)
  
  CovTheta = rbind(0,(cbind(0,CovTheta)))
  
  assert(SEtheta == sqrt(diag(CovTheta)))
  thetaStar = theta - mean(theta)
  
  # centering matrix
  A = diag((n-1)/n, n)
  A[upper.tri(A)] = -1/n
  A[lower.tri(A)] = -1/n
  
  assertnum((A %*% theta), thetaStar)
  
  CovThetaStar = (A %*% CovTheta %*% t(A))
  SEthetaStar = sqrt(diag(CovThetaStar))
  
  tab = data.frame(theta = theta, SEtheta = SEtheta,
                   thetaStar = thetaStar, SEthetaStar = SEthetaStar)
  return(tab)
}

tabOR = function(M, numeric = FALSE) {
  assert(ncol(M) == 2)
  logOr = M[,1]
  SElogOr = M[,2]
  
  if (!numeric) {
    tab = data.frame(OR = sprintf("%.2f",exp(logOr)),
                     CI = sprintf("from %.2f to %.2f", 
                                  exp(logOr-1.96*SElogOr),
                                  exp(logOr+1.96*SElogOr)),
                     z = sprintf("%.2f",logOr/SElogOr),
                     p = formatPval(2*pnorm(-abs(logOr/SElogOr)))
    )
    colnames(tab)[2:4] = c("95%-CI", "z-value", "p-value")
  } else if (numeric) {
    tab = data.frame(theta = logOr,
                     lower = logOr-1.96*SElogOr,
                     upper = logOr+1.96*SElogOr
    )
  }
  return(tab)
}

tabLm = function(M) {
  theta = M[,1]
  SEtheta = M[,2]
  tab = data.frame(estimate = sprintf("%.3f", theta),
                   CI = sprintf("from %.3f to %.3f", 
                                theta-1.96*SEtheta,
                                theta+1.96*SEtheta),
                   t = sprintf("%.2f",theta/SEtheta),
                   p = formatPval(2*pnorm(-abs(theta/SEtheta)))
  )
  colnames(tab)[2:4] = c("95%-CI", "t-value", "p-value")
  
  return(tab)
}


# converts OR to a SMD (Hedges' g)
or2smd = function(or, se, n1, n2) {
  
  d = log(or)*sqrt(3)/pi # Cohen's d
  Vd = (se^2)*3/pi^2 # Variance, Cooper p. 233
  
  # Hedges' g bias correction
  df = n1 + n2 - 2
  j = 1 - 3/(4*df - 1) # Cooper, p. 213
  j = pmin(j, 1) # for cases df is 0, -1 or -2 
  
  g = d*j # Hedges' g 
  SEg = sqrt(j^2*Vd)
  
  return(data.frame(g=g, SEg=SEg))
}

# converts SMD to r (Cohens' r)
smd2r = function(d, vd, n1, n2) {
  
  a = (n1 + n2)^2 / (n1*n2) # Cooper p. 234
  r = d / sqrt(d^2 + a)
  vr = (a^2*vd)/(d^2 + a)^3
  SEr = sqrt(vr)
  
  z = 0.5*log((1 + r) / (1 - r))
  SEz = sqrt(1 / sqrt(n1 + n2 - 3))
  
  return(data.frame(r=r, SEr=SEr, z=z, SEz=SEz))
}

# delta AIC
dAIC = function(fit) {
  AIC = drop1(fit)
  tab = data.frame(deltaAIC = AIC$AIC - AIC$AIC[1])
  rownames(tab) = rownames(AIC)
  return(tab)
}

# categorize bayes factors
catbf = function(bf) {
  cat = rep(NA, length(bf))
  
  idx = bf <= 1 & bf > 1/3
  cat[idx] = "weak"
  idx = bf <= 1/3 & bf > 1/10
  cat[idx] = "moderate"
  idx = bf <= 1/10 & bf > 1/30
  cat[idx] = "substantial"
  idx = bf <= 1/30 & bf > 1/100
  cat[idx] = "strong"
  idx = bf <= 1/100 & bf > 1/300
  cat[idx] = "very strong"
  idx = bf <= 1/300
  cat[idx] = "decisive"
  
  cat = factor(cat, levels = c("weak", "moderate", "substantial", 
                               "strong", "very strong", "decisive"))
  return(cat)
}

# categorize bayes factors with less categories
catbf.less = function(bf) {
  cat = rep(NA, length(bf))
  
  idx = bf <= 1 & bf > 1/10
  cat[idx] = "weak to moderate"
  idx = bf <= 1/10 & bf > 1/100
  cat[idx] = "substantial to strong"
  idx = bf <= 1/100
  cat[idx] = "very strong to decisive"
  
  cat = factor(cat, levels = c("weak to moderate", "substantial to strong", "very strong to decisive"))
  return(cat)
}