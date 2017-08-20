###############
## Author : Pooja Sethiya
## Institute : Chris Lab Faculty of Health Science / University of Macau. 
## Email : yb57662@umac.mo
## Date: 20/08/2017
###############
## This script gives you occurence of your motif of interest in multi fasta file with location and frequency
## Requirement: Fasta file of multiple sequence
## mymotifs.txt (Directly use output of getMotifCombination.R )
## motifCombination
# TTCAAGAAAAAAAAAGAA
# TTCTAGAAAAAAAAAGAA
# TTCGAGAAAAAAAAAGAA
# TTCCAGAAAAAAAAAGAA
##

rm(list=ls())


MotifInGene = function(SequenceFile, mymotifs, outfile){
  
          ##### Required Packages #####
          require(Biostrings)
          library("Rgraphviz")
          library(psych)
          require(ggplot2)
          require(reshape2)
          require(IRanges)
          require(kebabs)
          
          ###### Input sequence file
          CgPromoter=readDNAStringSet(SequenceFile,format = "fasta")
          print(head(CgPromoter))
          mymotifs = read.table(mymotifs,header=TRUE)
          mymotifs = as.matrix(mymotifs)
          print(head(mymotifs))
          #Convert DNAStringSet to DNAString

          dnaString <- lapply(CgPromoter,as.character)
          
          dnaString <- unlist(dnaString)
          dnaString <- c2s(dnaString)
          seq = DNAString(dnaString)
          
          ######### Make Combination of Kmers (can make combination more than 15 and maximum 21)
          kmer <- spectrumKernel(k=nchar(mymotifs[1]), normalized=FALSE)
          
          ########Calculate KmerFrequency in given sequence file
          
          system.time(kmerCount <- drop(getExRep(seq, kmer)))
          dat = as.data.frame(kmerCount)
          kmerCount = 0 #Empty the memory consuming object
     
          
          ###### subset frequency of motifs of interest
          
          
          
          FreqMyMotif <- subset(dat,rownames(dat) %in% mymotifs)
          
          FreqMyMotif.True =subset(FreqMyMotif,kmerCount>0)
          TrueMyMotifs = as.matrix(rownames(FreqMyMotif.True))
          
          cat("Genes with Motif:", length(TrueMyMotifs))

          ##Plot frequency
          
          bar1 = ggplot(FreqMyMotif.True,aes(rownames(FreqMyMotif.True),kmerCount))+
                    geom_bar(stat = "identity")+
                    geom_text(stat='identity',aes(label=kmerCount),vjust=-1,fontface="bold")+
                    theme(axis.text.x= element_text(face="bold", colour="black", size=12,angle=45,vjust=0.8),
                          axis.text.y = element_text(face="bold", color="black",size=12),
                          axis.title.y=element_text(face="bold", color="black",size=12))+
                    labs(x="",y="No. of Genes")
          
                                                    
          print(bar1)
          #########################################
   ll <- list()
   tt  <-  list()
             for(i in seq_along(TrueMyMotifs)){
                       # i=1
                       
                       mi0 <- vmatchPattern(TrueMyMotifs[i], CgPromoter)     
                       
                       
                       coords = as.data.frame(mi0)
                       
                       nmatch_per_seq <- elementNROWS(mi0)
                       
                       pos = which(nmatch_per_seq>0)
                       
                       tt[[i]] = table(nmatch_per_seq)
                       
                       Freq = nmatch_per_seq[coords$group]
                       
                       genes = names(mi0)[coords$group]
                       
                       start=coords[,3]
                       
                       ll[[i]] = as.data.frame(cbind(genes,Freq,start))
                      
             }
             
   names(ll)=TrueMyMotifs[,1]
   names(tt)=TrueMyMotifs[,1]
   
   melt_tt = melt(tt,id.vars=names(tt))
   
   dat.melt = melt_tt[(melt_tt$nmatch_per_seq >0),]
   head(dat.melt)
   
   
   library(dplyr)
   totals <- dat.melt %>%
             group_by(L1) %>%
             summarize(total = sum(value))
   
             bar2=ggplot(dat.melt)+ geom_bar(stat="identity") +  
                       aes(x = (L1), y = value, label = value, fill = factor(nmatch_per_seq),vjust=-0.8,fontface="bold")+
                       theme() +
                       geom_text(aes(L1, total, label = total, fill = NULL), data = totals)+
                       theme(axis.text.x= element_text(face="bold", colour="black", size=12,angle=70,vjust=0.8,hjust=0.59),
                             axis.text.y = element_text(face="bold", color="black",size=12),
                             axis.title.y=element_text(face="bold", color="black",size=12),
                             legend.title=element_text(face="bold", color="black",size=12),
                             legend.text=element_text(face="bold", color="black",size=10))+
                       labs(x="",y="No. of Genes")+
                       guides(fill = guide_legend(title = "Frequency per Sequence"))
   
            print(bar2)
   
            ##write in outfile
             for(i in 1:length(names(ll))){
                       
                       
                       print(names(ll[i]))
                       # write.xlsx(ll[i],file="Ca_TTCnnGAAnnTTC.xlsx",sheetName=names(ll[i]),row.names=F,append=T,col.names=T)
                       
                       write.table(ll[i], file=paste(outfile,".txt",sep=""),sep="\t",row.names=FALSE,quote =FALSE,append = T)
                       
             }
   
}  
   
