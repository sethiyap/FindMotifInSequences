###############
## Author : Pooja Sethiya
## Institute : Chris Lab Faculty of Health Science / University of Macau. 
## Date : 18_08_2017
## Email : yb57662@umac.mo
###############
##Script to make multiple combinatioon of motifs with N
##only for N i.e A,T,C,G
##Run script with motif sequence defined keeping the rest same as defined below.
##getMoitfCombination("TATATAn")
##1. input the motif
##############################################




getMoitfCombination=function(seq){

          #seq = c("TTCnnnnnnnTTC")
          index <- grep ("n", s2c(seq))
          
          
          x <- expand.grid(rep(list(c('A','T', 'G','C')),length(index)))
          
          
          seqR = as.data.frame(rep(seq,nrow(x)))
          head(as.data.frame(seqR),100)
          
                    ss = lapply(as.character(seqR[,1]),function(i){
                              
                              sc = s2c(i)
                              return(sc)
                    })
                    
          
          ss = do.call("rbind",ss)
          
                    for(i in 1:length(index)){
                              ss[,index[i]] <- as.character(x[,i])
                              
                    }
          
          
          
                    motifCombination = apply(ss,1,function(i) c2s(i))
                    
          motifCombination = data.frame(motifCombination)
          
          print(head(motifCombination))    
          
          write.table(motifCombination, file=paste(seq,"motif_combination.txt", sep="_"), quote = F,row.names = F)

}



