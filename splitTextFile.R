# made by Damon Poburko
# input: multi-entry genpept(full) or genebank(full) saved file from NCBI protein or nucleotide search
# output: a folder of .txt files for each individual accession number in the composite input file
# options: set cdsOnly as TRUE or FALSE to have nucleotide sequence output restricted to the CDS listed under features. 

library("readtext")
library("stringr")

## ======= USER INPUT ========================================================================
#setwd("C:\\Users\\Damon\\Documents\\splitTextFile")
#dir <- getwd()

fList <- choose.files(default = "", caption = "Select genpept formatted files",multi = TRUE, filters = Filters, index = nrow(Filters))
cdsOnly <- TRUE

pdList <- gregexpr("\\\\.",fList[1]) #split file path to vector separated by "\\"
currFile <- fList[1]
fName <- gsub("\\\\","", substr( currFile,pdList[[1]][length(pdList[[1]])], nchar(currFile)  ) )
pathOut <- paste0(substr( currFile,1, pdList[[1]][length(pdList[[1]])]))
fType = ""
if (gregexpr(".gb$",currFile)[[1]][1] > -1) fType = ".gb"
if (gregexpr(".gp$",currFile)[[1]][1] > -1) fType = ".gp"
  

# Read in all files from a folder
pathIn <- currFile

# create output folder
nFolderSuffix <- "_nucleotide"
if (cdsOnly == TRUE) nFolderSuffix <- paste0(nFolderSuffix,"-CDS")
if (fType == ".gp") oFolderPath <- paste0(pathOut,gsub(".gp","",fName),"_peptide","\\")
if (fType == ".gb" ) oFolderPath <- paste0(pathOut,gsub(".gb","",fName),nFolderSuffix,"\\")
if (!dir.exists(oFolderPath)){
  dir.create(oFolderPath, recursive = TRUE)
}

#check each line for new accession # or Source value to make name

# fName from Accession + source

an <- "ABC12345" #accession number
sr <- "Genus species" #source organism
fLines <-  c(noquote(""))
#fLines <- rbind(fLines,noquote("def"))
#write.table(fLines, pathOut, append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
# rbind approach works, but need to sort out how nut to write quotation marks around text
CDSnucs = ""


getCDS <- function(someText, cdsStart, cdsEnd) {
  #someText <- fLines
  #cdsStart <-1
  #cdsEnd <-1896
  seq <- ""
  r2r <- c() #rowsToReplace
  startNTS <- FALSE
  originLine <-0
  fwdSlashLine <- 0
  rowsKept <- matrix(data = TRUE, nrow = nrow(someText), ncol = 1 )
  
  for (j in 1:nrow(someText)) {
    if (unlist(gregexpr(pattern ='//',someText[j])) == 1 ) {
      startNTS = FALSE #check if reached end of sequence
      fwdSlashLine <- j
    } 
    if (startNTS == TRUE) {
      seq = paste0(seq,someText[j])
      r2r = cbind(r2r,j)
      rowsKept[j] <- FALSE
    } 
    #add this toggle after collection of NTS code block so that it starts on the line after "ORIGIN"
    if (unlist(gregexpr(pattern ='ORIGIN',someText[j])) == 1 ) {
      startNTS = TRUE
      originLine <-j
    } 
  }
  print(paste0("seq: ",seq))
  seq <- gsub("\\d+","",seq)        #remove numbers from seq
  seq <- gsub(" ","",seq)  #remove spaces and cast to upper case 
  seq <- substr(seq,cdsStart, cdsEnd) # get substring representing the CDS
  
  seqSplit <- c()
  nc <- nchar(seq)
  nr <- ceiling(nc/60)
  for (i in 1:nr) {
    subSeq <- substr(seq,1, min(60,nchar(seq)))
    seq <- substr(seq,61,nchar(seq))
    seqNum <- as.character(60*(i-1)+1)
    prefix <- paste0( strrep(" ",10-nchar(seqNum)),as.character(seqNum)," ")
    #add a space every 10 letters
    seqSplit <- rbind(seqSplit, paste0(prefix,gsub("(([a-z]){10})","\\1 ", subSeq))) 
  }
  
  #textOut = rbind(someText[1:originLine,1],seqSplit, someText[fwdSlashLine:nrow(someText),1])
  part1 <- someText[1:originLine]
  part2 <- someText[fwdSlashLine:nrow(someText)]
  textOut = c(part1, seqSplit, part2)  
  return(textOut)
}


processGenpept = function(pathIn) {

  con = file(pathIn, "r")
  hasCDS = FALSE
  hasORIGIN = FALSE
  passedFEATURES = FALSE
  filesSaved = 0
  cdsF <-0
  cdsL <-0
  
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    #print(line)
    if (unlist(gregexpr(pattern ='ACCESSION',line)) == 1 ) {
      an <- unlist(strsplit(line,"  "))
      an <- an[length(an)]
    }
    if (unlist(gregexpr(pattern ='SOURCE',line)) == 1 ) {
      sr <- unlist(strsplit(line,"  "))
      sr <- sr[length(sr)]
    }
    #line = "CDS             1..1896"
    if (unlist(gregexpr(pattern ='FEATURES',line)) == 1 ) {
      passedFEATURES = TRUE
    }
    
    
    
    if (unlist(gregexpr(pattern ='   CDS',line)) > 1 & passedFEATURES==TRUE ) {
      hasCDS = TRUE
      cds0 <- unlist(strsplit(line,"  "))
      cds0 <- cds0[length(cds0)]
      cds0 <- unlist(strsplit(cds0,"\\."))
      cdsF <- as.numeric(cds0[1])
      cdsL <- as.numeric(cds0[length(cds0)])
    }
    if (unlist(gregexpr(pattern ='ORIGIN',line)) == 1 ) {
      hasORIGIN = TRUE
    }
    
    fLines <- rbind(fLines,noquote(line)) #add current line to growing new file content
    
    #End of current file
    if (line == "") {

      # in the future, need to create optional function to read fLines and replace ORIGIN section with only the CDS
      # borrow getting of ORIGIN from motifSearch, 
    
      if (fType == ".gb" & cdsOnly==TRUE) fLines <- getCDS(fLines,cdsF,cdsL)

      pathOut <- paste0(oFolderPath,"/",an,"_",sr,".txt") 
      print(paste0("saving ",pathOut))
      write.table(fLines, pathOut, append = FALSE, quote = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
      filesSaved <- filesSaved+1

      fLines <- c("") #reset fLines for next file
      an <- "ABC12345" #accession number
      sr <- "Genus species" #source organism
      
    }
  
  }
  
  close(con)
}

 processGenpept(pathIn)

print("Done extracting files!")
