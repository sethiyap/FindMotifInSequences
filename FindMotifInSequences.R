##Script to make multiple combinatioon of motifs with wildcards
##Define the wildcards in hash
##Run script with motif sequence defined keeping the rest same as defined below.
##combination("TATATAR",1,c(),fastaSequence = "Sc_Stable_fiveprimeUTR.fa")
##1. input the motif
##2. fastasequence input file to search the motif

require("hash")
require("seqinr") 


          combination =  function(Sequenc,pos,prefix,fastaSequence){
          
                    
                    #Define the hash with wild cards
                     h = hash(W=c("A","T"), X=c("A","G","C"), Y=c("T","G"), R=c("A","G"),M=c("C","A"))
          
                     sequ=unlist(strsplit(Sequenc,""))
                     
          
                  
                    if(pos==length(sequ)+1){
                    
                    
                              motifcombination = paste(prefix,collapse = "")
                              print(motifcombination)
                              
                              
                              ##Calling function to find the motif combination in the fasta file
                              
                              findMotif(motifcombination,fastaSequence)
                              return()
                    
                    }
          
                    base = sequ[pos]
                    
          
                              if(any(keys(h)==base)){
                                        value = h[[base]]
                                       
                                        for(i in 1:length(value)){
                                                  
                                                  newSeq = append(prefix,value[i])
                                                  
                                                  combination(sequ,pos+1,newSeq,fastaSequence)
                                        }
                              }
                              else{
                                        newSeq = append(prefix, base)
                                        combination(sequ,pos+1,newSeq,fastaSequence)
                              }
}



          findMotif = function(motif,input_file){
                    
                    
                    myseq <- read.fasta(input_file, as.string = TRUE,forceDNAtolower = FALSE)
                    
                    
                    for(i in 1:length(myseq)){
                              
                              ##Condition to find the motif in fasta sequence
                              
                              if(grepl(motif,myseq[i], perl = T)==T){
                            
                                        
                               id = (attr(myseq[i],which="name"))
                               
                               ## Write the outputfile with ID's containing given motif in sequence
                               
                               write.table(id, file=paste("outfile",".txt"),quote=FALSE,append = TRUE, row.names = F, col.names = F)
                               
                               }
                    }
          }



