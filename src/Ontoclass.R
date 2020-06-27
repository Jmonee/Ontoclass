## for each uniprotID, get all of the GO numbers who are ancesstos to transporter activity..
#check if there ia a chebi mapping already existant
# if not get chebi
# get Chebi number for each term.
ws=Ontoclassdir
Gomapping<- read.csv(paste0(ws,"/mapping/GO2C.csv") ,stringsAsFactors=FALSE)
tcdbmapping<- read.csv(paste0(ws,"/mapping/TCDBMapping_substrates.csv") ,stringsAsFactors=FALSE)

dbpath=paste0(ws,"/db/tcdb")

############################################### Functions   ###############################################

getConciseCategory<- function(categoryTemp,cellValue)
{
  
  directCat<- NA
  categoryIndex<- which(unlist(lapply(categoryTemp, function(x)  !is.empty(x))))
  if(is.empty(categoryIndex))
  {
    directCat<- names(unlist(apply(directMappingdf, 2, function(x)  grep(paste0("^",AllGoChebiNames[[i]][j],"$"),x ,ignore.case=TRUE))))
    return(c(NA,directCat))
  }
  if(length(categoryIndex) ==1 ) #mapped to only one category
  {
    if(categoryIndex==1)
      directCat<- names(which(unlist(apply(directMappingdf, 2, function(x)  any(grepl(paste0("^",AllGoChebiNames[[i]][j],"$"),x ,ignore.case=TRUE)) ))))
    
    # print(directCat)
    #print(paste0(AllGoChebiIds[[i]][j],"  got mapped to one category: ", CategoryChebiMapping$Category[categoryIndex]))
    nCellValue=CategoryChebiMapping$Category[categoryIndex]
    return(c(nCellValue,directCat))
  } else{ #more than one category, get the most consice ones
    #print(paste0(AllGoChebiIds[[i]][j],"  got mapped to multiple categories: ", paste0(CategoryChebiMapping$Category[categoryIndex],collapse=" , ")))
    #print("optimizing...............")
    
    AllChebi<-unlist(categoryTemp[categoryIndex]) #all chebi that are related to [[i]][j] go term
    possiblePairs<- combn( unlist(categoryTemp[categoryIndex]), m=2) #get all possbile pairs
    pairsInTheSameCategory<-apply(possiblePairs, 2, function(x)  length(unique(CategoryIDAncestors[x])) == 1 )
    evaluateTheseColPairs<- which(!pairsInTheSameCategory ) #coloumIndex
    
    ignoreMe<- c() #contains all of the Chebi that were in a ansestor of other chebi (we dont need them because we found a better more concise representative)
    for( z in 1:length(evaluateTheseColPairs)) #length of pairs that do not belong to same category
    {
      #print(paste0("z",z))
      
      pair1<- possiblePairs[1,evaluateTheseColPairs[z]]
      pair2<- possiblePairs[2,evaluateTheseColPairs[z]]
      pair1Ancestors<- CategoryAncestors[pair1] 
      pair2Ancestors<- CategoryAncestors[pair2]
      print(pair1)
      print(pair2)
      
      if(any(c(pair1,pair2) %in% ignoreMe ))
        next
      
      if(pair1 %in% unlist(pair2Ancestors))         # if i am existant in any category parent remove me...
      {
        
       # print(paste0("ignore ",pair1))
        ignoreMe<-c(ignoreMe, pair1)
      }
      if(pair2 %in% unlist(pair1Ancestors))
      {
       # print(paste0("ignore ",pair2))
        ignoreMe<- c(ignoreMe, pair2)
      }
      
      
    }
    latestAbstractCategoryIndx<- unname(unlist(CategoryIDAncestors[AllChebi[which(! AllChebi %in% ignoreMe )]]))
    latestAbstractCategories<- paste0(unique(CategoryChebiMapping$Category[latestAbstractCategoryIndx]), collapse = ";")
    #print(paste0(AllGoChebiIds[[i]][j]," after optimizing  got mapped to multiple categories: ",latestAbstractCategories))
    if((length(unique(CategoryChebiMapping$Category[latestAbstractCategoryIndx]))> 1) | unique(latestAbstractCategoryIndx) ==1 ) #   #check if there is a multiple categories (UDP is mapped to nuclebase and sugar) 
      #and we can furthur limit it to only one using direct mapping  file 
    {
      print(paste0(" checking whether (",AllGoChebiNames[[i]][j], ") has direct mapping"))
      print(latestAbstractCategories)
      directCat<- names(unlist(apply(directMappingdf, 2, function(x)  grep(paste0("^",AllGoChebiNames[[i]][j],"$"),x ,ignore.case=TRUE))))
      print(directCat )
      if(length(directCat)> 1)
        print(paste0("Ops!  (",AllGoChebiNames[[i]][j], ") is mapped to multiple direct categories! ", directCat))
      
    }
    
    # if(cellValue == ""){
    #   cellValue<-latestAbstractCategories
    # } else{
    #   # remove duplicates
    # 
    #   cellValue<- paste(cellValue,latestAbstractCategories,sep =  ",")
    #  # print(cellValue)
    #   cellValue<- paste0(unique(trimws(unlist(strsplit(unlist(strsplit(cellValue,",")),";")))), collapse = ";")
    #   #print(cellValue)
    # }
    
    return(c(latestAbstractCategories,directCat))
  }
  
  
}

decideWhatToDo<- function(GoChebiCategory)
{
  # adjust according what you decide with the professor on what caregory to use (direct)
  
  Chebicategory<- unlist(strsplit(GoChebiCategory,","))
  
  
  decision<- rep(NA,length(Chebicategory))
  importantCat<- rep(NA,length(Chebicategory))
  for(y in 1:length(Chebicategory))
  {
    #decisionA<- getFinalCategory(Chebicategory[[x]][y],ChebiWdirectMapping[[x]][y], ChebiRole[[x]][y]) #if there is direct mapping, always choose it
    decisionB<- getRuleOfthumbCategeory(Chebicategory[y]) # basic rules of thumb
    decision[y]<-decisionB
    #dfGoChebiCategory$decision[x]<- getFinalCategory(dfGoChebiCategory$Chebicategory[x],dfGoChebiCategory$ChebiWdirectMapping[x], dfGoChebiCategory$ChebiRole[x])
  }
  
  
  
  return( getImportantCategory( unique(decision)))
}

getsubstrateCategory<- function(GoNumber)
{
  category<-NA
  index<- which(Gomapping$GO_Number %in% GoNumber)
  if(length(index)>=1) # this Go number is relate
  {
    category<- Gomapping$decision[index]
  }
  
  
  return (category)
  
}

getRuleOfthumbCategeory<- function(decisionB)
{
  #amino acid derivative;carboxylic acid;amide -> amino acid derivative
  
  # Category
  # 1.A Nonselective
  # 1.B water
  # 1.C inorganic cation 
  # 1.D inorganic anion 
  # 1.E Other inorganic
  # 2.A organic anion
  # 2.B organic cation
  # 3.A monosaccharide and derivative
  # 3.B oligosaccharide and derivative
  # 3.C polysaccharide and derivative
  # 4.A monocarboxylic acid
  # 4.B dicarboxylic acid
  # 4.C tricarboxylic acid
  # 5.A amino acid
  # 5.B amino acid derivative
  # 5.C peptide
  # 5.D amine
  # 5.E polyamin
  # 5.F protein
  # 5.G Other oganic amino compounds
  # 6.A nucleobase
  # 6.B nucleoside
  # 6.C nucleic acid
  # 6.D nucleotide
  # 7.A polyol
  # 7.B organophosphorus 
  # 7.C amide
  # 7.D Other organic
  # monosaccharide and derivative;nucleotide ->nucleotide
  #   monosaccharide and derivative;organophosphorus -> organophosphorus
  #other carbs;polyol
 # print(174)
  splitted<- trimws(unlist(strsplit(decisionB,";")))
  
  if(length(splitted)>1 &  ("1.A Nonselective" %in% splitted) )
    splitted=splitted[-which(splitted=="1.A Nonselective")]
  
  if(all(c("7.B organophosphorus","3.A monosaccharide and derivative") %in% splitted))
  {
    splitted=splitted[-which(splitted=="3.A monosaccharide and derivative")]
  }
  if(all(c("6.D nucleotide","3.A monosaccharide and derivative") %in% splitted))
  {
    splitted=splitted[-which(splitted=="3.A monosaccharide and derivative")]
  }
  if(all(c("7.A polyol","7.D Other organic") %in% splitted))
  {
    splitted=splitted[-which(splitted=="7.D Other organic")]
  }
  if(all("5.B amino acid derivative" %in% splitted, any(c("4.A monocarboxylic acid", "4.B dicarboxylic acid", "4.C tricarboxylic acid") %in% splitted))) 
  {
    splitted[-which(splitted=="carboxylic acid")]
    splitted=splitted[- which(splitted %in% c("4.A monocarboxylic acid", "4.B dicarboxylic acid", "4.C tricarboxylic acid"))]
  }
  if(all(c("5.B amino acid derivative","7.C amide") %in% splitted))
  {
    splitted=splitted[-which(splitted=="5.B amino acid derivative")]
  }
  if(all(c("5.B amino acid derivative","5.A amino acid") %in% splitted))
  {
    splitted=splitted[-which(splitted=="5.A amino acid")]
  }
  if(all(c("5.A amino acid","7.C amide") %in% splitted))
  {
    splitted=splitted[-which(splitted=="5.A amino acid")]
  }
  return(paste0(splitted, collapse = ";"))
}
getImportantCategory<- function(decision)
{
  splitted<- trimws(unlist(strsplit(decision,",")))
  done=FALSE
  if(length(splitted)>1 &  ("1.A Nonselective" %in% splitted) )
    splitted=splitted[-which(splitted=="1.A Nonselective")]
  
  while(! done)
  {
    
    
    if(length(  splitted )>=2 )
    {
      if(all("unspecificed anion" %in% splitted, any(c("1.D inorganic anion", "2.A organic anion") %in% splitted))) 
      {
        splitted=splitted[- which(splitted =="unspecificed anion")]
        next
      }
      
      if(all("unspecificed cation" %in% splitted, any(c("1.C inorganic cation") %in% splitted))) 
      {
        splitted=splitted[- which(splitted =="unspecificed cation")]
        next
      }
      
      if(all(splitted %in% c("unspecificed cation", "1.D inorganic anion")))
      {
        splitted=splitted[-which(splitted %in% c("unspecificed cation", "1.D inorganic anion"))]
        if(length(splitted)==0)
        {
          splitted="MC"
          done= TRUE
        }
        next
      }
      
      
      if(all(splitted %in% c("unspecificed anion", "1.C inorganic cation")))
      {
        splitted=splitted[-which(splitted %in% c("unspecificed anion", "1.C inorganic cation"))]
        if(length(splitted)==0)
        {
          splitted="MC"
          done= TRUE
        }
        next
      }
      if(any(splitted %in% c("1.C inorganic cation", "1.D inorganic anion")))
      {
        splitted=splitted[-which(splitted %in% c("1.C inorganic cation", "1.D inorganic anion"))]
        if(length(splitted)==0)
        {
          splitted="MC"
          done= TRUE
        }
        next
        
        
      }
      
      if("7.D Other organic" %in% splitted )
        splitted=splitted[-which(splitted =="7.D Other organic")]
      done= TRUE
      
    }
    else
      done=TRUE
  }
  return(paste0(splitted, collapse = ","))
  
}

getOntoclassOutput<- function(uniprotID)
{
Seqdetails<-list(UniProt.ID=NA,AC=NA, OC=NA,OS=NA,KW=NA,GoMF=NA,Evidance=NA,SCL=NA, Fasta=NA)

Seqdetails$UniProt.ID=uniprotID
murl<- paste0("https://www.uniprot.org/uniprot/",uniprotID,".txt")
print(murl)
data<- readLines(url(murl))
listm<- strsplit(data,"  ")

# xm <-lapply(listm, function(x) if (any(c("AC","KW","DR","OS", "OC", "CC") %in% x) ){ print(x)})
names(listm)<- lapply(listm, function(x)  x[1])
#AC numbers
numbers<- c()
for (i in which(names(listm)=="AC"))
{
  numbers<-c(numbers,unlist(listm[[i]])[2])
}
ACnumbers<- paste0(numbers, collapse = "")
if (ACnumbers=="") Seqdetails$AC<- NA else Seqdetails$AC<- ACnumbers

OC<- strsplit(listm[[min(which (names(listm) =="OC"))]],",")[[2]][1]
Seqdetails$OC<- OC
OS<- strsplit(listm[[min(which (names(listm) =="OS"))]],",")[[2]][1]
Seqdetails$OS<- OS
# print(OC)
#print(OS)
# for extracting keywords
#print(309)
Keywords<- c()
for (i in which(names(listm)=="KW"))
{
  details<- unlist(listm[[i]])[2]
  Keywords<-c(Keywords,details)
}
Keywords<- paste0(Keywords, collapse = "")
if (Keywords=="") Seqdetails$KW<- NA else Seqdetails$KW<- Keywords
#print(Keywords)


# for extracting GO MF
GOMFNoIEA<-  c()
Evidance<- c()
GoMF<-c()
for (i in which(names(listm)=="DR"))
{
  details<- unlist(listm[[i]])[2]
  detailsSplited <- unlist(strsplit(details,";"))
  if(" GO" %in% detailsSplited )
  {
    if(grepl("F:",detailsSplited[3],fixed=TRUE))
    {
      GoNumber<- detailsSplited[2]
      GoDescription<- detailsSplited[3]
      EvidanceCode<- detailsSplited[4]
      Evidance<- c(Evidance,EvidanceCode)
      
      #print(GoNumber)
      #print(GoDescription)
      #print(EvidanceCode)
      GoMF<- c(GoMF, paste0(GoNumber,";",GoDescription,";",EvidanceCode))
    }
  }
}
GOMFNoIEA<- paste0(GOMFNoIEA,collapse = "," )
GoMF<- paste0(GoMF, collapse = ",")
Evidance<- paste0(Evidance, collapse = ",")
if (GoMF=="") Seqdetails$GoMF<- NA else Seqdetails$GoMF<- GoMF
if (Evidance=="") Seqdetails$Evidance<- NA else Seqdetails$Evidance<- Evidance




#CC   -!- SUBCELLULAR LOCATION:
SCL<- c()
inSCL<-FALSE
for (i in which(names(listm)=="CC"))
{
  # print(listm[[i]])
  details<- unlist(listm[[i]])[2]
  
  if (inSCL)
  {
    #  print("still in SCL")
    #  print(unlist(listm[[i]])[4])
    SCL<- c(SCL,unlist(listm[[i]])[4])
  }
  
  if(grepl("-!-", details,fixed=TRUE))
  {
    # print("first if")
    if(grepl("-!- SUBCELLULAR LOCATION:",details,fixed=TRUE))
    { 
      inSCL=TRUE
      SCL <- c(SCL,unlist(strsplit(details,"-!- SUBCELLULAR LOCATION:"))[2])
      #print("second if ")
    } else {
      # print("else")
      inSCL= FALSE
    }
  }
}
SCL<- paste0(SCL[-(length(SCL))], collapse = "")
if (SCL=="") Seqdetails$SCL<- NA else Seqdetails$SCL<- SCL
# print(SCL)


Seqdetails$Fasta<- str_replace_all(string=trimws(paste(unlist(listm[seq(which(names(listm)=="SQ")+1,     (which(names(listm)=="//")-1))]),collapse = "")),pattern=" ", repl="")
#print(391)

write.fasta(Seqdetails$Fasta,Seqdetails$UniProt.ID,paste0(ws,Seqdetails$UniProt.ID,".fasta"))
closeAllConnections()

closeAllConnections()
#print(df)
#print(Seqdetails$GoMF)

if(! is.na(Seqdetails$GoMF))
Goterms<-unlist(strsplit(Seqdetails$GoMF,"GO:"))
numberofGOterms<- length(Goterms)-1
categories<- c()
categoriesEvidance<- c()
for(j in 2:(numberofGOterms+1)) # each go term in that sequence
{
  substrateCateg<- getsubstrateCategory(paste0("GO:",str_extract(Goterms[j], "[0-9]+")))
  
  #print(407)
  
  if (!is.na(substrateCateg))
  {
    evidance<- trimws(sub(",","",unlist(strsplit(Goterms[j],";"))[3]))
    categories<- c(categories,substrateCateg)
    categoriesEvidance<- c(categoriesEvidance, evidance)
  }
  
  
}
if(!is.null(categories))
{
  
  uniqueIndex<- (which(!duplicated(categories)))
  if (length(uniqueIndex)> 1 & ("1.A. Nonselective" %in% trimws(categories[uniqueIndex]) ))
  {
    uniqueIndex<- uniqueIndex [ - which(trimws(categories[uniqueIndex])  %in% "1.A. Nonselective")]
    print("***")
    print(categories[uniqueIndex])
  }
  
  Seqdetails$category<- paste0(categories[uniqueIndex],"*",categoriesEvidance[uniqueIndex] ,collapse = ";")
  Seqdetails$Chebicategory<- paste0(categories[uniqueIndex],collapse = ",")
  
}
if(length(Seqdetails$Chebicategory)>0)
Seqdetails$decision<-decideWhatToDo(Seqdetails$Chebicategory)
############################################### TCDB Mapping for the output   ###############################################
#makeblastdb -dbtype prot -in  tcdb    -parse_seqids

system(paste0("blastp -comp_based_stats 1 -db ", dbpath,  " -query  ", ws,Seqdetails$UniProt.ID,".fasta -outfmt 6  -num_alignments 1  -out  ",ws,Seqdetails$UniProt.ID,"TCDB_blast.txt"))

blastresults<- read.table(paste0(ws,Seqdetails$UniProt.ID,"TCDB_blast.txt"))
blastresults<- blastresults[which(!duplicated(blastresults[,1])),]
Substrate<- rep("", length(blastresults[,1]))
AllQueryname<- sub("\\|.*","",sub(".+?\\|","",blastresults[,1]))
SIndex<- which(grepl(".+?\\.", AllQueryname))
Substrate[SIndex] <- sub("\\..*","", AllQueryname[SIndex])
Substrate[SIndex] <- sub("\\T","",Substrate[SIndex])
AllQueryname[SIndex]<-sub(".+?\\.","", AllQueryname[SIndex])


HitTCDBID<-  sub("\\|.*","",sub(".+?\\|","",blastresults[,2]))
Exact<- rep(FALSE, length(blastresults[,1]))
Hit<- rep(FALSE, length(blastresults[,1]))
Exact[which(blastresults$V3==100.00 & blastresults$V11==0)]<- TRUE
Hit[which(blastresults$V3> 40 & blastresults$V11<1e-20)]<- TRUE

id <- read.table(paste0(ws,"db/acc2tcid.tsv"),sep="\t", quote = "")
TCDBFam<-  rep(NA, length(blastresults[,1]))
TCDBFam[Hit]<- as.vector(id[match(HitTCDBID[Hit], id[,1]),2])
Substrate= tcdbmapping$Substrates[which(tcdbmapping$TCDBID==TCDBFam)]

Seqdetails$TCDBexact<-Exact

Seqdetails$TCDBHit<-Hit

Seqdetails$TCDBFam=TCDBFam
Seqdetails$tcdbSubstrate=paste(Substrate,collapse = ";")
closeAllConnections()
#remove intermediate files
file.remove(paste0(ws,Seqdetails$UniProt.ID,".fasta"))
file.remove(paste0(ws,Seqdetails$UniProt.ID,"TCDB_blast.txt"))
return (Seqdetails)
}

