library(RCurl)
library(XML)

minmax <- function(X, minX, maxX){
  (X-minX)/(maxX-minX)
}

shift0 <- function(X){
  X - min(X)
}

catCount <- function(X){
  sapply(levels(as.factor(X)), function(cat){ sum(X == cat) })
}

getSMILES.CIDs <- function(CIDs){
  uri <- 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/IsomericSMILES/XML'
  params <- c(cid=paste(CIDs, collapse=','))
  body <- paste('cid=', paste(CIDs, collapse=','), sep='')
  res <- postForm(uri, .params=params, style='post', .contentEncodeFun=curlPercentEncode)
  doc <- xmlParse(res)
  root <- xmlRoot(doc)
  props <- xmlElementsByTagName(root, 'Properties', recursive=TRUE)
  pnames <- c()
  smis <- c()
  for(p in props){
    cid <- xmlElementsByTagName(p, 'CID', recursive=TRUE)[[1]]
    smi <- xmlElementsByTagName(p, 'IsomericSMILES', recursive=TRUE)[[1]]
    pnames <- c(pnames,xmlValue(cid))
    smis <- c(smis, xmlValue(smi))
  }
  names(smis) <- pnames
  return(smis)
}
