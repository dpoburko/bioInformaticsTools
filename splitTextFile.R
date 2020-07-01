
# made by Damon Poburko

library("readtext")
library("stringr")

## ======= USER INPUT ========================================================================
#setwd("C:\\Users\\Damon\\Documents\\splitTextFile")
dir <- getwd()

# Read in all files from a folder
#pathIn <- paste0(dir,"/sequence.txt")
#pathIn <- paste0(dir,"/solute carrier family 17 member 9 .gp")
pathIn <- paste0(dir,"/sialin_spp.gp")
#txt <- readtext(paste0(dir,"/sequence.txt"))

# create output folder
oFolderName <- "SLC17A5"
oFolderPath <- paste0(dir,"/",oFolderName,"/")
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

processFile = function(pathIn) {
  con = file(pathIn, "r")
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
    #Need some function to write lines to file. 
    fLines <- rbind(fLines,noquote(line))
    
    #End of current file
    if (line == "" ) {
      #reset ab and sr
      #save to file
      #reset table
      
      pathOut <- paste0(oFolderPath,"/",an,"_",sr,".txt") 
      print(paste0("saving ",pathOut))
      write.table(fLines, pathOut, append = FALSE, quote = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
      fLines <- c("")
      an <- "ABC12345" #accession number
      sr <- "Genus species" #source organism
    }
  
  }
  
  close(con)
}

processFile(pathIn)

