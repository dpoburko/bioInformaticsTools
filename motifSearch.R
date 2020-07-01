# Built by Damon Poburko, June 2020, Simon Fraser University
# build your search pattern
# 0 (zero) hydroph"0"bics, @ for acids, b for basic, "." wildcard
srchTerm = c("@X00")  

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

## ====
#phi =  V, I, L, F, W, Y and M
# Aromatic amino acids are represented by ‘Ω’ (Omega) and they include F, W, Y.
# Hydrophilic amino acids ζ = E, D, K, R and H
# Large aliphatic side chains ‘Ψ’ (Psi);V, I, L and M
# Amino acids with small chains π = P, G, A and S
# uncharged: 9 =  S, T, N, Q

## ====
hydrophobics <- c("A","I","L","M","F","W","Y","V") # hydrophobic: 0 =  A, I, L, M, F, W, Y, V 
acids <- c("D","E") # negative/acid charge: @ = D, E 
bases <- c("R","H", "K") # positive charge: + = R, H, K

#Prompt user to select genpept formatted text file from NCBI
fList <- choose.files(default = "", caption = "Select files",multi = TRUE, filters = Filters, index = nrow(Filters))

#Example sequences for working with
#vglut2<- as.character("MESVKQRILAPGKEGIKNFAGKSLGQIYRVLEKKQDNRETIELTEDGKPLEVPEKKAPLCDCTCFGLPRRYIIAIMSGLGFCISFGIRCNLGVAIVDMVNNSTIHRGGKVIKEKAKFNWDPETVGMIHGSFFWGYIITQIPGGYIASRLAANRVFGAAILLTSTLNMLIPSAARVHYGCVIFVRILQGLVEGVTYPACHGIWSKWAPPLERSRLATTSFCGSYAGAVIAMPLAGILVQYTGWSSVFYVYGSFGMVWYMFWLLVSYESPAKHPTITDEERRYIEESIGESANLLGAMEKFKTPWRKFFTSMPVYAIIVANFCRSWTFYLLLISQPAYFEEVFGFEISKVGMLSAVPHLVMTIIVPIGGQIADFLRSKQILSTTTVRKIMNCGGFGMEATLLLVVGYSHTRGVAISFLVLAVGFSGFAISGFNVNHLDIAPRYASILMGISNGVGTLSGMVCPIIVGAMTKNKSREEWQYVFLIAALVHYGGVIFYALFASGEKQPWADPEETSEEKCGFIHEDELDEETGDITQNYINYGTTKSYGATSQENGGWPNGWEKKEEFVQESAQDAYSYKDRDDYS")
#vglut1 <- as.character("MEFRQEEFRKLAGRALGKLHRLLEKRQEGAETLELSADGRPVTTQTRDPPVVDCTCFGLPRRYIIAIMSGLGFCISFGIRCNLGVAIVSMVNNSTTHRGGHVVVQKAQFSWDPETVGLIHGSFFWGYIVTQIPGGFICQKFAANRVFGFAIVATSTLNMLIPSAARVHYGCVIFVRILQGLVEGVTYPACHGIWSKWAPPLERSRLATTAFCGSYAGAVVAMPLAGVLVQYSGWSSVFYVYGSFGIFWYLFWLLVSYESPALHPSISEEERKYIEDAIGESAKLMNPLTKFSTPWRRFFTSMPVYAIIVANFCRSWTFYLLLISQPAYFEEVFGFEISKVGLVSALPHLVMTIIVPIGGQIADFLRSRRIMSTTNVRKLMNCGGFGMEATLLLVVGYSHSKGVAISFLVLAVGFSGFAISGFNVNHLDIAPRYASILMGISNGVGTLSGMVCPIIVGAMTKHKTREEWQYVFLIASLVHYGGVIFYGVFASGEKQPWAEPEEMSEEKCGFVGHDQLAGSDDSEMEDEAEPPGAPPAPPPSYGATHSTFQPPRPPPPVRDY")
#vnutrt <- as.character("MPSQRSSLMQPIPEETRKTPSAAAEDKRWSRPECQLWTGMLLLGTCLLYCTRVTMPVCTVAMSQDFGWNKKEAGIVLSSFFWGYCLTQVVGGHLGDRIGGEKVILLSASAWGFITVTTPLLAHLGSGHLAFVTFSRILTGLLQGVYFPALTSLLSQRVQESERSFTYSTVGAGSQVGTLVTGGIGSVLLDRCGWQSVFYFSGGLTLLWVYYVYKYLLDEKDLVLALGVLAQGLPVTRPSKVPWRQLFRKASVWAVICSQLSSACSFFILLSWLPTFFKETFPHSKGWVFNVVPWLLAIPASLFSGFISDRLISQGYRVITVRKFMQVMGLGLSSIFALCLGHTTSFLKSMIFASASIGFQTFNHSGISVNIQDLAPSCAGFLFGVANTAGALAGVVGVCLGGYLIETTGSWTCVFHLVAIVSNLGLGTFLVFGKAQRVDLVPTHEDL")

# Function to search a sequence (subject) for a pattern (querry)
aasearch <-function(querry, subject) {
  
  subjectO <- subject #make copy of subject to reference original sequence before replacements
  querry <- gsub("X",".",querry)
  
  # find and replace hydrophobic residues with "0"
  if (length((grep("0", querry)))==1) {   
    for (i in 1:length(hydrophobics)) {
      subject <- gsub(hydrophobics[i],"0",subject)
    }
  }
  # find and replace acidic residues with "@"
  if (length((grep("@", querry)))==1) {
    for (i in 1:length(acids)) {
      subject <- gsub(acids[i],"@",subject)
    }
  }
  
  c <- gregexpr(querry,subject) #look for querry in subject
  loci <- c[[1]]      #extract starting positions of matches from list c
  locus <- c[[1]][1]  
  
  for (i in 1:length(c[[1]])) {
    fi = c[[1]][i]                      #first index of matched string
    li = c[[1]][i] + nchar(querry)-1  #last index of matched string
    str
    if (i==1) {
      ss <- c(substring(subjectO,fi,li))   #grab the sequence that matches the pattern
      sn <- c(paste0(fi,"..",li))         #grab the positions that match the pattern
    } else {
      ss <- rbind(ss, c(substring(subjectO,fi,li)))
      sn <- rbind(sn,c(paste0(fi,"..",li)))
    }
  }
  rr = cbind(ss,sn)  
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
      #reset ab and sr
      #save to file
      #reset table
      
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

for (i in 1:length(fList)) {
  
  plot.progress(i/length(fList))
  pf <- parseFile(fList[i])[1]
  currFile = fList[i]

  pf <- parseFile(currFile)
  rm(searchResult)
  searchResult <- aasearch(srchTerm, unlist(pf[1]))
  searchResult <- cbind(searchResult,unlist(strsplit(pf[2], " "))[2],unlist(pf[3]))   #add accession numbesr and species common names to results table
  
  #generate a fileName and path for the output
  pdList <- gregexpr("\\\\.",currFile) #split file path to vector separated by "\\"
  fName <- gsub("\\\\","", substr( currFile,pdList[[1]][length(pdList[[1]])], nchar(currFile)  ) )
  srchTermOut <- gsub("0" ,"[phi]" ,gsub("@","[DE]",srchTerm))  #replace search code with file friendly version
  fileOut <- paste0(Sys.Date(),"_",srchTermOut,".csv")
  pathOut <- paste0(substr( currFile,1, pdList[[1]][length(pdList[[1]])]),fileOut)
  #print(paste0("saving ",pathOut))
  write.table(searchResult, pathOut, append = TRUE, quote = FALSE, sep = "|", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"))
  
}
doneMsg = paste0(length(fList)," files searched for ", srchTerm, ". Search saved.")
print(doneMsg)




## ===== list of asci character ===============

# \\dDigit, 0,1,2 ... 9
# \\DNot Digit
# \\sSpace
# \\SNot Space
# \\wWord
# \\WNot Word
# \\tTab
# \\nNew line
# ^Beginning of the string
# $End of the string
# \Escape special characters, e.g. \\ is "\", \+ is "+"
# |Alternation match. e.g. /(e|d)n/ matches "en" and "dn"
# •Any character, except \n or line terminator
# [ab]a or b
# [^ab]Any character except a and b
# [0-9]All Digit
# [A-Z]All uppercase A to Z letters
# [a-z]All lowercase a to z letters
# [A-z]All Uppercase and lowercase a to z letters
# i+i at least one time
# i*i zero or more times
# i?i zero or 1 time
# i{n}i occurs n times in sequence
# i{n1,n2}i occurs n1 - n2 times in sequence
# i{n1,n2}?non greedy match, see above example
# i{n,}i occures >= n times
# [:alnum:]Alphanumeric characters: [:alpha:] and [:digit:]
# [:alpha:]Alphabetic characters: [:lower:] and [:upper:]
# [:blank:]Blank characters: e.g. space, tab
# [:cntrl:]Control characters
# [:digit:]Digits: 0 1 2 3 4 5 6 7 8 9
# [:graph:]Graphical characters: [:alnum:] and [:punct:]
# [:lower:]Lower-case letters in the current locale
# [:print:]Printable characters: [:alnum:], [:punct:] and space
# [:punct:]Punctuation character: ! " # $ % & ' ( ) * + , - . / : ; < = > ? @ [ \ ] ^ _ ` { | } ~
# [:space:]Space characters: tab, newline, vertical tab, form feed, carriage return, space
# [:upper:]Upper-case letters in the current locale
# [:xdigit:]Hexadecimal digits: 0 1 2 3 4 5 6 7 8 9 A B C D E F a b c d e f

## ________________


