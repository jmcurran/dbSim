processTextDB = function(Lines){
  nLines = length(Lines)
  Lines = Lines[nchar(Lines)>0]
  Lines = gsub("^[[:space:]]+(.*)$", "\\1", Lines)
  warningMessages = NULL
  
  if(length(Lines) < nLines){
    cat(paste("Dropped", nLines - length(Lines), "empty lines\n"))
    nLines = length(Lines)
  }
  
  alleleLineNumber = grepl("^[Aa]l*lell*e[as]*,", Lines)
  
  if(!any(alleleLineNumber)){
    errMsg = paste("The input file must contain the allele column", "and be labelled 'Allele'", sep = "\n")
    stop(errMsg)
  }else{
    alleleLineNumber = which(alleleLineNumber)
    
    if(length(alleleLineNumber) > 1){
      msg = "There is more than one line that has allele information. Using the first"
      warningMessages = c(warningMessages, msg)
      warning(msg)
      alleleLineNumber = alleleLineNumber[1]
    }
  }
  
  
  countLineNumber = grepl('^2*[Nn]',Lines)
  
  
  if(!any(countLineNumber)){
    cat(paste("Last line: ", Lines[nLines], "\n"))
    stop("One of the last line of the input file must be labelled n or N")
  }else{
    countLineNumber = which(countLineNumber)
    
    if(length(countLineNumber)>1){
      msg = "There is more than one line that has count information"
      warningMessages = c(warningMessages, msg)
      warning(msg)
      countLineNumber = countLineNumber[1]
    }
  }
  
  locusLine = Lines[alleleLineNumber]
  countLine = Lines[countLineNumber]
  
  if(alleleLineNumber > 1){
    Lines = Lines[-(1:alleleLineNumber)]
    nDropped = length(1:alleleLineNumber)
    countLineNumber = countLineNumber - nDropped
    nLines = nLines - nDropped
  }else{
    Lines = Lines[-1]
    countLineNumber = countLineNumber - 1
    nLines = nLines - 1
  }
  
  if(countLineNumber < nLines){
    Lines = Lines[-(countLineNumber:nLines)]
  }else{
    Lines = Lines[-nLines]
  }
  
  if(length(Lines) < 5){
    stop(paste("There doesn't appear to be any data in the file:", fileName))
  }
  
  ## drop any lines without a valid allele
  i = !grepl("^[0-9.]+,", Lines)
  
  if(any(i)){
    idx = which(i)
    msg = paste("Dropping lines:", paste(idx, collapse = ","), "because the allele names are invalid")
    warningMessages = c(warningMessages, msg)
    warning(msg)
    Lines = Lines[-idx]
  }
  
  locusLine = sub("^[Aa][^,]*,","", locusLine)
  
  div = 1
  if(grepl("^2", countLine)){
    div = 2
  }
  
  countLine = sub("^2*[Nn][^,]*,","", countLine)
  
  Loci = unlist(strsplit(locusLine,","))
  Loci = Loci[nchar(Loci) > 0]
  nLoci = length(Loci)
  
  if(is.null(multiplex)){
    Loci = cleanLoci(Loci)
  }else{
    Loci = cleanLoci(Loci, multiplex)
  }
  
  unusedLoci = Loci$unused
  
  #fileName = sub('^.*/([^/]*$)', '\\1', fileName)
  
  if(any(nchar(Loci)==0))
    stop("Cannot have empty locus labels")
  
  
  
  Counts = ceiling(as.numeric(unlist(strsplit(countLine, ",")))[1:nLoci] / div)
  
  if(any(Counts<=0)){
    stop("All counts must be >= 0\n")
  }
  
  nLoci = length(Loci$loci[Loci$keep])
  freqs = vector(length = nLoci, mode = "list")
  names(freqs) = Loci$loci[Loci$keep]
  
  for(i in 1:nLoci)
    freqs[[i]] = list(a = NULL, f = NULL)
  
  ## strip out anything that is not a number, period or comma
  
  Lines = gsub("[^0-9.,]+", "", Lines)
  
  lineCtr = 2
  keep = Loci$keep
  for(line in Lines){
    Tokens = unlist(strsplit(line, ','))
    Allele = as.numeric(Tokens[1])
    Tokens = Tokens[-1]
    Tokens = Tokens[keep]
    
    if(length(Tokens) != nLoci){
      errMsg = paste("Incorrect number of fields on line ", lineCtr)
      stop(errMsg)
    }
    
    idx = which(regexpr("^[0-9.]+", Tokens)!=-1)
    
    for(i in idx){
      fx = as.numeric(Tokens[i])
      
      if(fx < 0 | fx > 1)
        stop("Frequencies must be between 0 and 1")
      
      if(fx > 0){
        freqs[[i]]$a = c(freqs[[i]]$a, Allele)
        freqs[[i]]$f = c(freqs[[i]]$f, fx)
      }
    }
    
    lineCtr = lineCtr + 1
  }
  
  invisible(list(Loci = Loci$loci[Loci$keep], Counts = Counts[keep], 
                 Freqs = freqs,
                 warnings = warningMessages, 
                 unusedLoci = unusedLoci
  ))
}

readDB = function(fileName){
  f1 = file(fileName, "r")
  
  if(!isOpen(f1))
    stop(paste0("Couldn't open ", fileName, " for reading\n"))
  
  Lines = readLines(f1)
  close(f1)
  
  return(processTextDB(Lines))
}
