# Built by Damon Poburko, June 2020, Simon Fraser University
# Use case: 
# - search a selection of genpept formatted text files for a peptide motif

#build your motif search term (srchTerm):
# Basic regex searches can be emploed

# CODING AMINO ACID GROUPS:
# hydroph"0"bic residues: 0 (zero) = "(A|I|L|M|F|W|Y|V)"
# acidic residues: @ = "(D|E)"
# basic residues: b = "(R|H|K)"
# single wild card:X (this get converted to regex wildcard ".") 
# i{n1,n2}i occurs n1 - n2 times in sequence
# see this link for peptide nomenclature https://www.sciencedirect.com/science/article/pii/S0014579301032951   "Normalization of nomenclature for peptide motifs as ligands of modular protein domains"
# see this page for regex search tips    https://cran.r-project.org/web/packages/stringr/vignettes/regular-expressions.html    
# adding ^ restricts a sequence to be search from the end of a string, e.g. "^X{6,15}AXR0D", AXR0DX must be within 6-15 residues from the N-terminus
# adding $ restricts a sequence to be search from the end of a string, e.g. "AXR0DX{6,15}$", AXR0DX must be within 6-15 residues from the C-terminus

#srchTerm = c("(R|K)(L|V|I)X{4,5}(H|Q)(L|A)")  
#srchTerm = c("@X{2,3}00")  
#srchTerm = c("AXR0DX{6,15}$") # putative VNUT C-tail motifputativeUniqueVNUTtail


# a comma separated list of motifs to search for
srchTerm = c("(G|V)T(C|G)(L|V)LY(C|S)","@X{2,3}00") #pupativeVNUT-1stTM

# a comman separated list of matched motif labels
fNameSuffix = c("putativeLDCV","putativeDileucine")  

# ===================== End of user input ======================================================================
#Prompts user to select genpept formatted text file from NCBI
fList <- choose.files(default = "", caption = "Select genpept formatted files",multi = TRUE, filters = Filters, index = nrow(Filters))

plot.progress <- function(...)	{
  
  vectOfBar0 <- c(...)
  progressTitle <- vectOfBar0[1]
  vectOfBar <- as.numeric(vectOfBar0[2:length(vectOfBar0)])*100
  tleft = c(0,0,0)
  numOfBar <- length(vectOfBar)
  barNames <- c("files","motif(s)")
  plot(c(0,100), c(0,numOfBar), type='n', xlab='', ylab='', yaxt='n', mar=c(3,3,3,3))
  for(i in 1:numOfBar) {
    rect(0, 0.1+i-1, vectOfBar[i], 0.9+i-1, col=cm.colors(numOfBar)[i])
    text(0.5, 0.5+i-1, paste( barNames[i], ': ', round(vectOfBar[i],2), '%', sep=''), adj=0)
  }
  title(paste0(progressTitle,'...'))
}

# Function to search a sequence (subject) for a pattern (querry)
aasearch <-function(querry, subject) {

  #Manual tested for trouble shooting
  subjectO <- subject #make copy of subject to reference original sequence before replacements
  querry <- gsub("0","(A|I|L|M|F|W|Y|V)",querry)
  querry <- gsub("@","(D|E)",querry)
  querry <- gsub("B","(R|H|K)",querry)
  querry <- gsub("X",".",querry)

  c <- gregexpr(querry,subject ) #look for querry in subject
  loci <- c[[1]]      #extract starting positions of matches from list c
  locus <- c[[1]][1]  
  matchedLengths <-attr(c[[1]],"match.length",exact=FALSE)
  
  if (locus==-1) {
        ss <- c("none")   #grab the sequence that matches the pattern
        sn <- c(paste0(0,"..",0))         #grab the positions that match the pattern
        fis = 0
  } else {
    for (i in 1:length(c[[1]])) {
      fi <- c[[1]][i]                      #first index of matched string
      li <- c[[1]][i] + matchedLengths[i]-1  #last index of matched string
      
      if (i==1) {
        ss <- c(substring(subjectO,fi,li))   #grab the sequence that matches the pattern
        sn <- c(paste0(fi,"..",li))         #grab the positions that match the pattern
        fis <- fi
      } else {
        ss <- rbind(ss, c(substring(subjectO,fi,li)))
        sn <- rbind(sn,c(paste0(fi,"..",li)))
        fis <-rbind(fis,fi)
      }
    }
  }
  rr <- cbind(ss,sn,fis)  
  return(rr) 
}

# Function to read amino acid sequence from a genpep formatted text file
addToOrganism = FALSE

parseFile = function(pathIn) {

    #pathIn <- fList[2]
    con <- file(pathIn, "r")
    seq <- as.character("")
    addToSeq <- FALSE
    an <- "ABC12345" #accession number
    sr <- "Genus species" #source #organism
    dfn <- "sequence definition"
    sl <- "1" #sequence length 
    taxon <- "00000"
    orgTxt <- as.character("")
    lineNum = 0
    
  while ( TRUE ) {
    lineNum <- lineNum+1

    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    #print(paste0("line ",lineNum,": ",line))
    if (unlist(gregexpr(pattern ='ACCESSION',line)) == 1 ) {
      an <- unlist(strsplit(line,"  "))
      an <- an[length(an)]
    }
    if (unlist(gregexpr(pattern ='SOURCE',line)) == 1 ) {
      sr <- unlist(strsplit(line,"  "))
      sr <- sr[length(sr)]
    }
    if (unlist(gregexpr(pattern ='DEFINITION',line)) == 1 ) {
      #line = "Vesicular glutamate transporter 1, partial [Ophiophagus hannah]."
      dfn <- unlist(strsplit(line,"  "))
      dfn <- dfn[length(dfn)]
      cutAt <- gregexpr("\\[",dfn)[[1]]
      dfn <- substr(dfn,1,cutAt-2)
      dfn <- gsub(",","",dfn)
    }
    if (unlist(gregexpr(pattern ='LOCUS',line)) == 1 ) {
      b = strsplit(strsplit(line," aa ")[[1]][1],"  ")[[1]]
      sl = gsub(" ","",b[length(b)])
    }
    if (unlist(gregexpr(pattern ='taxon:',line)) >= 1 ) {
      bb = strsplit(strsplit(line,":")[[1]][2],"  ")[[1]]
      taxon =  gsub("[^[:digit:].]", "", bb)
    }

    #toggle collection of taxonomy "off" if line ends with "." or subsequent keywords
    if (unlist(gregexpr(pattern ='\\.$',line)) == 1 ||  unlist(gregexpr(pattern ='REFERENCE',line)) == 1
        || unlist(gregexpr(pattern ='COMMENT',line)) == 1 || unlist(gregexpr(pattern ='FEATURES',line)) == 1) {
      addToOrganism = FALSE
    }
    #Copy taxonomy of organisms to output
    if (addToOrganism == TRUE) {
      #orgTxt = paste0(orgTxt,as.character(line))
      orgTxt = gsub(" ","", paste0(orgTxt,as.character(line)))
    }
    #toggle collection of taxonomy "on" for next lap(s) if line contains "ORGANISM
    if (unlist(gregexpr(pattern ='ORGANISM',line)) > 0 ) {
      addToOrganism = TRUE
    }

    #add current line to growing sequence if you have passed the "ORIGIN" line in the file
    if(addToSeq == TRUE) {
      seq <- paste0(seq, as.character(line))
    }   
    
    #check if you have passed the "ORIGIN" line in the file. Next loop will start to build sequence.
    if (unlist(gregexpr(pattern ='ORIGIN',line)) == 1 ) {
      addToSeq <- TRUE
    }
    if (unlist(gregexpr(pattern ='//',line)) == 1 ) {
      addToSeq <- FALSE
    }
    #End of current file
    seq <- gsub("\\d+","",seq)        #remove numbers from seq
    seq <- toupper(gsub(" ","",seq))  #remove spaces and cast to upper case 
    seq <- toupper(gsub("//","",seq))  #remove spaces and cast to upper case 
    
  } # close while(TRUE) block
    
  close(con)
  return(c(seq,an,sr,dfn,sl,orgTxt,taxon))  
}

rm(searchResult)
rm(searchResultsAll)

for (st in 1:length(srchTerm)) {
  
  srchTerm1 <- toupper(srchTerm[st]) #ensure that srchTerm is upper case
  fNameSuffix1 <- fNameSuffix[st]
  

  headers <-c("aa","positions","aaStart","aaFromEnd","accession","spp","definition","seq length","taxonomy","taxon")
  searchResultsAll <- headers
  
  #generate a fileName and path for the output
  pdList <- gregexpr("\\\\.",fList[1]) #split file path to vector separated by "\\"
  currFile <- fList[1]
  fName <- gsub("\\\\","", substr( currFile,pdList[[1]][length(pdList[[1]])], nchar(currFile)  ) )
  srchTermOut <- gsub("B","[RKH]",gsub("0" ,"[phi]" ,gsub("@","[DE]",srchTerm1)))  #replace search code with file friendly version
  srchTermOut <- gsub("\\|" ,"" ,srchTermOut)  #replace search code with file friendly version
  fileOut <- paste0(Sys.Date(),"_",srchTermOut,fNameSuffix1,".csv")
  pathOut <- paste0(substr( currFile,1, pdList[[1]][length(pdList[[1]])]),fileOut)

  #delete previous versions of the same file
  if (file.exists(pathOut)) file.remove(pathOut)
  
    
  for (i in 1:length(fList)) {
      plot.progress(srchTerm1,i/length(fList),st/length(srchTerm))
      pf <- parseFile(fList[i])[1]
      currFile <- fList[i]
      pf <- parseFile(currFile)
      rm(searchResult)
      searchResult <- aasearch(srchTerm1, unlist(pf[1]))
      aaFromEnd <- as.numeric(matrix( data = unlist(pf[5]), nrow =length(as.numeric(searchResult[,3]) ) , ncol=1, byrow=FALSE))  - as.numeric(searchResult[,3])
      searchResult <- cbind(searchResult,aaFromEnd, unlist(strsplit(pf[2], " "))[2],unlist(pf[3]),unlist(pf[4]),unlist(pf[5]),unlist(pf[6]),unlist(pf[7]))   #add accession numbers and species common names to results table
      searchResultsAll <- rbind(searchResultsAll,searchResult)
      colnames(searchResult) <- headers
      if (i==1) writeColNames = TRUE
      if (i>1) writeColNames = FALSE
      write.table(searchResult, pathOut, append = TRUE, quote = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = writeColNames, qmethod = c("escape", "double"))
  }
  colnames(searchResultsAll) <- headers
  #write.table(searchResultsAll[2:nrow(searchResultsAll),], pathOut, append = TRUE, quote = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
  doneMsg <- paste0(length(fList)," files searched for ", srchTerm1, ". Search saved.")
  print(searchResultsAll[,1:6])
  print(doneMsg)

}


