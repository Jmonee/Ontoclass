#! /usr/bin/Rscript

###################################################
## Name: main.R
## input: the uniport identifiers of the proteins of intrest in .csv one per line.
## output: automatically assigned class/es (if it finds one) for the substate of transporters, and other related information.
## reference: Alballa, Munira, and Gregory Butler. "Ontology-based transporter substrate annotation for benchmark datasets." 2019 IEEE International Conference on Bioinformatics and Biomedicine (BIBM). IEEE, 2019.
## Author: Munira Alballa
##################################################


args <- commandArgs(trailingOnly=TRUE)

terminate <- FALSE

out <- "."
Ontoclassdir <- "."
db<-"./db/"
for(i in args){
  arg = strsplit(i, "=")[[1]];
  
  switch(arg[1],
         "-query"={
           query <- arg[2]
         },
         "-out"={
           out <- arg[2]
         },
         "-Ontoclass"={
           Ontoclassdir <- normalizePath(arg[2])
         },

         "-help"={
           cat("Ontoclass v2.0 (June. 2020)\n")
           cat("\n")
           cat("Usage: Ontoclass -query=<input> [-Ontoclass=<Ontoclassdir>] [-out=<outdir>] [-db=<database path>]\n")
           cat("\n")
           cat("\t<input> is your sequence input file in fasta format\n")
           cat("\t<out> is the output directory where you want the predicted results, formatted as csv\n")
           cat("\t\t<out> defaults to '",out,"'\n")
           cat("\t<Ontoclassdir> is the directory where the base TooT-SC files are located")
           cat("\t\t<Ontoclassdir> defaults to '",Ontoclassdir,"'\n")
           cat("\n")
           terminate <- TRUE
           break
         }
  )
}

if(!terminate) {
  
  if(!exists("query")) {
    stop("-query has not been passed")
  }

  
  
  test_file <- normalizePath(path.expand(query))
  resultspath <- paste0(normalizePath(path.expand(out)),"/")
  
  
  library("Biostrings")
  library("stringr")
  library("rjson")
  library("seqinr")
  library("ontologyIndex")
  wd=normalizePath(path.expand(".")) # change the the tool directory
  
 
  
#  Ontoclassdir="/Users/muniraalballa/Dropbox/gitTest/OntoClass/"
#  resultspath=paste0(Ontoclassdir,"output/")
#  testfile=paste0(Ontoclassdir,"test.csv")
  
  #testing data with unknown substrates
  source(paste0(Ontoclassdir,"/src/Ontoclass.R"))
  IDs<-read.csv(testfile, header = T, stringsAsFactors = F)
  
  output=data.frame(UniProt.ID=IDs[,1], AC=NA,OC=NA, KW=NA, GoMF=NA, Evidance=NA, SCL=NA, Fasta=NA, OS=NA, category=NA, Chebicategory=NA, decision=NA, CDBexact=NA
                    ,TCDBHit=NA ,TCDBFam=NA,TCDBSubstrate=NA, stringsAsFactors = F)
  
  
for(i in 1: length(output$UniProt.ID))
{
  
 output[i,]<-  getOntoclassOutput(output$UniProt.ID[i])
}
  

  write.csv(output,paste0(resultspath,"out.csv"))
  
}

