# Built by Damon Poburko, June 2020, Simon Fraser University
# Use case: 
# - search a selection of genpept formatted text files for a peptide motif

#build your motif search term (srchTerm):
# We are working on incorporating regex searches
# coding amino acid groups:
# hydroph"0"bic residues: 0 (zero) = "(A|I|L|M|F|W|Y|V)"
# acidic residues: @ = "(D|E)"
# basic residues: b = "(R|H|K)"
# single wild card:X (this get converted to regex wildcard ".") 
# i{n1,n2}i occurs n1 - n2 times in sequence
# see this page for regex search tips https://cran.r-project.org/web/packages/stringr/vignettes/regular-expressions.html    

#srchTerm = c("(R|K)(L|V|I)XXXX(H|Q)(L|A)")  
#srchTerm = c("PXXXXX(R|K)X{2,4}(A|I|L)")  
#srchTerm = c("(R|K)(L|V|I)X{4,5}(H|Q)(L|A)")  
#srchTerm = c("B{2,3}X{9,11}B{3,4}")  
#srchTerm = c("YXX0")
srchTerm = c("EXXXXX00")  

fNameSuffix = "_dileucineLike"

plot.progress <- function(...)	{
  vectOfBar <- c(...)*100
  tleft = c(0,0,0)
  numOfBar <- length(vectOfBar)
  plot(c(0,100), c(0,numOfBar), type='n', xlab='', ylab='', yaxt='n', mar=c(3,3,3,3))
  for(i in 1:numOfBar) {
    rect(0, 0.1+i-1, vectOfBar[i], 0.9+i-1, col=cm.colors(numOfBar)[i])
    text(0.5, 0.5+i-1, paste( i, ': ', round(vectOfBar[i],2), '%', sep=''), adj=0)
  }
  title('Progress...')
}

# "Normalization of nomenclature for peptide motifs as ligands of modular protein domains"
# https://www.sciencedirect.com/science/article/pii/S0014579301032951

#Prompt user to select genpept formatted text file from NCBI
fList <- choose.files(default = "", caption = "Select genpept formatted files",multi = TRUE, filters = Filters, index = nrow(Filters))

#Example sequences for working with
#vglut2<- as.character("MESVKQRILAPGKEGIKNFAGKSLGQIYRVLEKKQDNRETIELTEDGKPLEVPEKKAPLCDCTCFGLPRRYIIAIMSGLGFCISFGIRCNLGVAIVDMVNNSTIHRGGKVIKEKAKFNWDPETVGMIHGSFFWGYIITQIPGGYIASRLAANRVFGAAILLTSTLNMLIPSAARVHYGCVIFVRILQGLVEGVTYPACHGIWSKWAPPLERSRLATTSFCGSYAGAVIAMPLAGILVQYTGWSSVFYVYGSFGMVWYMFWLLVSYESPAKHPTITDEERRYIEESIGESANLLGAMEKFKTPWRKFFTSMPVYAIIVANFCRSWTFYLLLISQPAYFEEVFGFEISKVGMLSAVPHLVMTIIVPIGGQIADFLRSKQILSTTTVRKIMNCGGFGMEATLLLVVGYSHTRGVAISFLVLAVGFSGFAISGFNVNHLDIAPRYASILMGISNGVGTLSGMVCPIIVGAMTKNKSREEWQYVFLIAALVHYGGVIFYALFASGEKQPWADPEETSEEKCGFIHEDELDEETGDITQNYINYGTTKSYGATSQENGGWPNGWEKKEEFVQESAQDAYSYKDRDDYS")
#vglut1 <- as.character("MEFRQEEFRKLAGRALGKLHRLLEKRQEGAETLELSADGRPVTTQTRDPPVVDCTCFGLPRRYIIAIMSGLGFCISFGIRCNLGVAIVSMVNNSTTHRGGHVVVQKAQFSWDPETVGLIHGSFFWGYIVTQIPGGFICQKFAANRVFGFAIVATSTLNMLIPSAARVHYGCVIFVRILQGLVEGVTYPACHGIWSKWAPPLERSRLATTAFCGSYAGAVVAMPLAGVLVQYSGWSSVFYVYGSFGIFWYLFWLLVSYESPALHPSISEEERKYIEDAIGESAKLMNPLTKFSTPWRRFFTSMPVYAIIVANFCRSWTFYLLLISQPAYFEEVFGFEISKVGLVSALPHLVMTIIVPIGGQIADFLRSRRIMSTTNVRKLMNCGGFGMEATLLLVVGYSHSKGVAISFLVLAVGFSGFAISGFNVNHLDIAPRYASILMGISNGVGTLSGMVCPIIVGAMTKHKTREEWQYVFLIASLVHYGGVIFYGVFASGEKQPWAEPEEMSEEKCGFVGHDQLAGSDDSEMEDEAEPPGAPPAPPPSYGATHSTFQPPRPPPPVRDY")
#vnutrt <- as.character("MPSQRSSLMQPIPEETRKTPSAAAEDKRWSRPECQLWTGMLLLGTCLLYCTRVTMPVCTVAMSQDFGWNKKEAGIVLSSFFWGYCLTQVVGGHLGDRIGGEKVILLSASAWGFITVTTPLLAHLGSGHLAFVTFSRILTGLLQGVYFPALTSLLSQRVQESERSFTYSTVGAGSQVGTLVTGGIGSVLLDRCGWQSVFYFSGGLTLLWVYYVYKYLLDEKDLVLALGVLAQGLPVTRPSKVPWRQLFRKASVWAVICSQLSSACSFFILLSWLPTFFKETFPHSKGWVFNVVPWLLAIPASLFSGFISDRLISQGYRVITVRKFMQVMGLGLSSIFALCLGHTTSFLKSMIFASASIGFQTFNHSGISVNIQDLAPSCAGFLFGVANTAGALAGVVGVCLGGYLIETTGSWTCVFHLVAIVSNLGLGTFLVFGKAQRVDLVPTHEDL")

# Function to search a sequence (subject) for a pattern (querry)
aasearch <-function(querry, subject) {

  #Manual tested for trouble shooting
  #subject <- c("RKMYAEYMKRYSSL")
  subjectO <- subject #make copy of subject to reference original sequence before replacements
  
  #querry = ("YXX0")
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
  } else {
    for (i in 1:length(c[[1]])) {
      fi <- c[[1]][i]                      #first index of matched string
      li <- c[[1]][i] + matchedLengths[i]-1  #last index of matched string
      
      if (i==1) {
        ss <- c(substring(subjectO,fi,li))   #grab the sequence that matches the pattern
        sn <- c(paste0(fi,"..",li))         #grab the positions that match the pattern
      } else {
        ss <- rbind(ss, c(substring(subjectO,fi,li)))
        sn <- rbind(sn,c(paste0(fi,"..",li)))
      }
    }
  }
  rr <- cbind(ss,sn)  
  return(rr) 
}

# Function to read amino acid sequence from a genpep formatted text file

parseFile = function(pathIn) {

    #pathIn <- fList[1]
    con <- file(pathIn, "r")
    seq <- as.character("")
    addToSeq <- FALSE
    an <- "ABC12345" #accession number
    sr <- "Genus species" #source organism
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
    #fLines <- rbind(fLines,noquote(line))

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
    if (line == "" ) {

      # pathOut <- paste0(oFolderPath,"/",an,"_",sr,".txt") 
      # print(paste0("saving ",pathOut))
      # write.table(fLines, pathOut, append = FALSE, quote = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
      # fLines <- c("")
      # an <- "ABC12345" #accession number
      # sr <- "Genus species" #source organism
    }
    seq <- gsub("\\d+","",seq)        #remove numbers from seq
    seq <- toupper(gsub(" ","",seq))  #remove spaces and cast to upper case 
    
  }
  close(con)
  return(c(seq,an,sr))  
}
rm(searchResult)
rm(searchResultsAll)
searchResultsAll <- c("aa","positions","accession","spp")

#generate a fileName and path for the output
pdList <- gregexpr("\\\\.",fList[1]) #split file path to vector separated by "\\"
fName <- gsub("\\\\","", substr( currFile,pdList[[1]][length(pdList[[1]])], nchar(currFile)  ) )
srchTermOut <- gsub("B","[RKH]",gsub("0" ,"[phi]" ,gsub("@","[DE]",srchTerm)))  #replace search code with file friendly version
srchTermOut <- gsub("\\|" ,"" ,srchTermOut)  #replace search code with file friendly version
fileOut <- paste0(Sys.Date(),"_",srchTermOut,fNameSuffix,".csv")
pathOut <- paste0(substr( currFile,1, pdList[[1]][length(pdList[[1]])]),fileOut)

for (i in 1:length(fList)) {
    plot.progress(i/length(fList))
    pf <- parseFile(fList[i])[1]
    currFile <- fList[i]
    pf <- parseFile(currFile)
    rm(searchResult)
    searchResult <- aasearch(srchTerm, unlist(pf[1]))
    searchResult <- cbind(searchResult,unlist(strsplit(pf[2], " "))[2],unlist(pf[3]))   #add accession numbers and species common names to results table
    searchResultsAll <- rbind(searchResultsAll,searchResult)

    write.table(searchResult, pathOut, append = TRUE, quote = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))
}

doneMsg <- paste0(length(fList)," files searched for ", srchTerm, ". Search saved.")
print(searchResultsAll)
print(doneMsg)



