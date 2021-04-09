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
#srchTerm = c("@X{2,3}00")  
#To Test
#YX{2,3}00 dileucine-like
#B{2,3}X{9,11}B{3,4}_nucBipartite
#(R|K)(L|V|I)X{4,5}(H|Q)(L|A)
#H@@0 (vnut tail?)
#2020-07-13_PXXXXX(RK)X{2,4}(AILV)
#EX{4,5}(L|I)(L|I)
#@X{3,5}00 dileucine like
#@X00
#srchTerm = c("@X{3,5}00") #dileucine-like  
#srchTerm = c("@@@@S@S@") # VMAT2 LDCV localization https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2897747/ 
#srchTerm = c("@XX(S|T)XXXDL$") # putative VNUT C-tail motif

#http://elm.eu.org/elms/elmPages/TRG_PTS1.html [SAPTC][KRH][LMFI]$)|([KRH][SAPTC][NTS][LMFI]
#srchTerm = c("(S|A|P|T|C)(K|R|H)(L|M|F|I)$")   #PST1 - confirmed syntax works. None in 295 VNUT
#srchTerm = c("(K|R|H)(S|A|P|T|C)(N|T|S)(L|M|F|I)$")   #PST1 - confirmed syntax works. None in 295 VNUT
#srchTerm = c("B@@0X$") # akin to HEDL. Purely mammalian
#srchTerm = c("B@@0X{1,25}$") # akin to HEDL at c-term but a little ways away - none in SLC17A9  
#srchTerm = c("@XXX00") # putative VNUT C-tail motif
fNameSuffix = "putativeGolgiRetention"
#srchTerm = c("AXR0DX{6,15}") # putative VNUT C-tail motifputativeUniqueVNUTtail
#srchTerm = c("@X{2,3}00") # putative VNUT C-tail motif
#srchTerm = c("(G|V)T(C|G)(L|V)LY(C|S)") pupativeVNUT-1stTM
#srchTerm = c("(D|E)(D|E)S(D|E)S(D|E)") # "putativeLDCV"
srchTerm = c("(F|L)(L|I)(I|V)XX(R|K)") # "putativeGolgiRetention"

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
srchTerm <- toupper(srchTerm) #ensure that srchTerm is upper case

#Example sequences for working with
#vglut2<- as.character("MESVKQRILAPGKEGIKNFAGKSLGQIYRVLEKKQDNRETIELTEDGKPLEVPEKKAPLCDCTCFGLPRRYIIAIMSGLGFCISFGIRCNLGVAIVDMVNNSTIHRGGKVIKEKAKFNWDPETVGMIHGSFFWGYIITQIPGGYIASRLAANRVFGAAILLTSTLNMLIPSAARVHYGCVIFVRILQGLVEGVTYPACHGIWSKWAPPLERSRLATTSFCGSYAGAVIAMPLAGILVQYTGWSSVFYVYGSFGMVWYMFWLLVSYESPAKHPTITDEERRYIEESIGESANLLGAMEKFKTPWRKFFTSMPVYAIIVANFCRSWTFYLLLISQPAYFEEVFGFEISKVGMLSAVPHLVMTIIVPIGGQIADFLRSKQILSTTTVRKIMNCGGFGMEATLLLVVGYSHTRGVAISFLVLAVGFSGFAISGFNVNHLDIAPRYASILMGISNGVGTLSGMVCPIIVGAMTKNKSREEWQYVFLIAALVHYGGVIFYALFASGEKQPWADPEETSEEKCGFIHEDELDEETGDITQNYINYGTTKSYGATSQENGGWPNGWEKKEEFVQESAQDAYSYKDRDDYS")
#vglut1 <- as.character("MEFRQEEFRKLAGRALGKLHRLLEKRQEGAETLELSADGRPVTTQTRDPPVVDCTCFGLPRRYIIAIMSGLGFCISFGIRCNLGVAIVSMVNNSTTHRGGHVVVQKAQFSWDPETVGLIHGSFFWGYIVTQIPGGFICQKFAANRVFGFAIVATSTLNMLIPSAARVHYGCVIFVRILQGLVEGVTYPACHGIWSKWAPPLERSRLATTAFCGSYAGAVVAMPLAGVLVQYSGWSSVFYVYGSFGIFWYLFWLLVSYESPALHPSISEEERKYIEDAIGESAKLMNPLTKFSTPWRRFFTSMPVYAIIVANFCRSWTFYLLLISQPAYFEEVFGFEISKVGLVSALPHLVMTIIVPIGGQIADFLRSRRIMSTTNVRKLMNCGGFGMEATLLLVVGYSHSKGVAISFLVLAVGFSGFAISGFNVNHLDIAPRYASILMGISNGVGTLSGMVCPIIVGAMTKHKTREEWQYVFLIASLVHYGGVIFYGVFASGEKQPWAEPEEMSEEKCGFVGHDQLAGSDDSEMEDEAEPPGAPPAPPPSYGATHSTFQPPRPPPPVRDY")
#vnutrt <- as.character("MPSQRSSLMQPIPEETRKTPSAAAEDKRWSRPECQLWTGMLLLGTCLLYCTRVTMPVCTVAMSQDFGWNKKEAGIVLSSFFWGYCLTQVVGGHLGDRIGGEKVILLSASAWGFITVTTPLLAHLGSGHLAFVTFSRILTGLLQGVYFPALTSLLSQRVQESERSFTYSTVGAGSQVGTLVTGGIGSVLLDRCGWQSVFYFSGGLTLLWVYYVYKYLLDEKDLVLALGVLAQGLPVTRPSKVPWRQLFRKASVWAVICSQLSSACSFFILLSWLPTFFKETFPHSKGWVFNVVPWLLAIPASLFSGFISDRLISQGYRVITVRKFMQVMGLGLSSIFALCLGHTTSFLKSMIFASASIGFQTFNHSGISVNIQDLAPSCAGFLFGVANTAGALAGVVGVCLGGYLIETTGSWTCVFHLVAIVSNLGLGTFLVFGKAQRVDLVPTHEDL")

# Function to search a sequence (subject) for a pattern (querry)
aasearch <-function(querry, subject) {

  #Manual tested for trouble shooting
  #subject <- c("RKMYAEYMKRYSSLPKMN")
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
headers <-c("aa","positions","aaStart","aaFromEnd","accession","spp","definition","seq length","taxonomy","taxon")
searchResultsAll <- headers

#generate a fileName and path for the output
pdList <- gregexpr("\\\\.",fList[1]) #split file path to vector separated by "\\"
currFile <- fList[1]
fName <- gsub("\\\\","", substr( currFile,pdList[[1]][length(pdList[[1]])], nchar(currFile)  ) )
srchTermOut <- gsub("B","[RKH]",gsub("0" ,"[phi]" ,gsub("@","[DE]",srchTerm)))  #replace search code with file friendly version
srchTermOut <- gsub("\\|" ,"" ,srchTermOut)  #replace search code with file friendly version
fileOut <- paste0(Sys.Date(),"_",srchTermOut,fNameSuffix,".csv")
pathOut <- paste0(substr( currFile,1, pdList[[1]][length(pdList[[1]])]),fileOut)

#delete previous versions of the same file
if (file.exists(pathOut)) file.remove(pathOut)


for (i in 1:length(fList)) {
    plot.progress(i/length(fList))
    pf <- parseFile(fList[i])[1]
    currFile <- fList[i]
    pf <- parseFile(currFile)
    rm(searchResult)
    searchResult <- aasearch(srchTerm, unlist(pf[1]))
    #STOPPED HERE FIGURE WHY TAXONOMY not being listed
    
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
doneMsg <- paste0(length(fList)," files searched for ", srchTerm, ". Search saved.")
print(searchResultsAll[,1:6])
print(doneMsg)



