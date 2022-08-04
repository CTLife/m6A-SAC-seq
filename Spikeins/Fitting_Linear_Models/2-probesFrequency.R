#! /usr/bin/env Rscript
###########################################################################################################
### 2-probesFrequency.R: Statistic absolute frequency and Mutation frequency.
### Author: Yong Peng
### Run "./2-probesFrequency.R  -h"  or "Rscript 2-probesFrequency.R  -h" to get some help.
###########################################################################################################





## To run the script in command lines.
#################################################
suppressPackageStartupMessages( library(optparse) )      

getParameters_f <- function() {
  option_list_Local <- list(   # Options list with associated default value.   
      optparse::make_option(opt_str=c("-f", "--file"),
			default="spikein_1_only-spikein.txt",
			type="character",   dest="file",
			help="Input file. [default: %default]."),

      optparse::make_option(opt_str=c("-o", "--outDir"),
			default="spikein_1_only-spikein",
			type="character",   dest="outDir",
			help="Path or dir name for output files. [default: %default].")
  )

  # now parse the command line to check which option is given and get associated values
  parser_Local <- optparse::OptionParser(usage="usage: %prog [options]",
		option_list=option_list_Local, 
		description="2-probesFrequency.R: Statistic absolute frequency2 and Mutation frequency2 from the result of script Four_spikeins.pl. Junuary 2, 2020.",                             
		epilogue="For comments, bug reports etc..., please contact Yong Peng <yongp@outlook.com>"
  )

  opt_Local <- optparse::parse_args(parser_Local, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options
  return(opt_Local)   
}


opt_g = getParameters_f()
inFile_g     <- opt_g$file
outDir_g     <- opt_g$outDir
# inFile_g     <- "1-merge/WT/merged.spikein-5.seq"
# outDir_g     <- "2-probesFrequency"
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g,   recursive = TRUE) }

options(digits=10)
continue_on_error <- function() {
  print( "NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())' " )
}
options( error=continue_on_error )  ## This option is very important.


library(tidyverse)
library(stringr)
library(ComplexHeatmap)  
library(corrplot) 
library(gplots)
library(spgs)
#############################################





#############################################
matrix1 = read.table(file= inFile_g, header = FALSE, sep = character(0), quote = "\"'", dec = ".")
dim(matrix1)

my_4mer = c()
my_4mer_name = c()
length(my_4mer)
length(my_4mer_name)
bases4 = c("A", "G", "C", "T")
index1 = 0
for(i1 in bases4) {
  for(i2 in bases4) { 
    for(i3 in bases4) { 
      for(i4 in bases4) {
        index1 = index1 + 1
        my_4mer[ index1 ] = paste(i1, i2, "[AGCT]", i3, i4, sep="")
        my_4mer_name[ index1 ] = paste(i1, i2, "A", i3, i4, sep="")  
      }
    }
  }
}
my_4mer
my_4mer_name 
length(my_4mer)
length(my_4mer_name)

vector1 = str_sub(string=matrix1[,1], start = 11, end = 15)
## vector1 = reverseComplement(x=vector1, content="dna", case="upper" )    ##RC  ##RC

my_list = vector(mode = "list", length = length(my_4mer_name) )
length(vector1)
length(my_list)


for(i in c(1:length(my_4mer_name)) ) {
  my_list[[i]] = vector1[str_detect(string = vector1, pattern = my_4mer[i])]
}
## str_detect(string = vector1, pattern = my_4mer[1])
## my_list

absolutefrequency2 = matrix(nrow=4, ncol=length(my_4mer))
rownames(absolutefrequency2) = c("A", "G", "C", "T")
colnames(absolutefrequency2) = my_4mer_name

Mutationfrequency2 = absolutefrequency2
num_eachKmer = Mutationfrequency2[1,]

for(i in c(1:length(my_4mer_name)) ) {
  temp1 = my_list[[i]]
  num_eachKmer[i] = length(temp1)  
  absolutefrequency2[1,i] = sum(str_sub(string=temp1, start = 3, end = 3) =='A')
  absolutefrequency2[2,i] = sum(str_sub(string=temp1, start = 3, end = 3) =='G')
  absolutefrequency2[3,i] = sum(str_sub(string=temp1, start = 3, end = 3) =='C')
  absolutefrequency2[4,i] = sum(str_sub(string=temp1, start = 3, end = 3) =='T')  
  Mutationfrequency2[1,i] = absolutefrequency2[1,i]/num_eachKmer[i]
  Mutationfrequency2[2,i] = absolutefrequency2[2,i]/num_eachKmer[i]
  Mutationfrequency2[3,i] = absolutefrequency2[3,i]/num_eachKmer[i]
  Mutationfrequency2[4,i] = absolutefrequency2[4,i]/num_eachKmer[i]
  mysum = Mutationfrequency2[1,i] + Mutationfrequency2[2,i] + Mutationfrequency2[3,i] + Mutationfrequency2[4,i]
  print(mysum)
}
#sum(str_sub(string=my_list[[1]], start = 3, end = 3) =='A')
#sum(str_sub(string=my_list[[1]], start = 3, end = 3) =='G')
#sum(str_sub(string=my_list[[1]], start = 3, end = 3) =='C')
#sum(str_sub(string=my_list[[1]], start = 3, end = 3) =='T')
print( sum(absolutefrequency2) )

write.table(x=absolutefrequency2, file = paste(outDir_g, "1_Absolute_Frequency.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

write.table(x=Mutationfrequency2, file = paste(outDir_g, "2_Mutation_Frequency.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

write.table(x=num_eachKmer,       file = paste(outDir_g, "3_num_eachKmer.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")


my_col1=colorRampPalette( c( "blue", "yellow3" , "pink", "red", "red"  ), bias = 1 )
Mutationfrequency22 = Mutationfrequency2*100
Mutationfrequency22 = round(x=Mutationfrequency22, digits = 0)

pdf( file=paste(outDir_g, "4_Mutation_Frequency.heatmap.pdf", sep="/"),  height=5,   width=60 )
  heatmap.2(x=Mutationfrequency22,           
           # dendrogram control
           dendrogram = "none", 
           Rowv = FALSE,
           Colv="Rowv" , 
           symm = FALSE,           
           # data scaling
           scale =  "none" ,
           na.rm=TRUE,           
           # colors
           col=my_col1(10),
           trace = "none",            
           # cell labeling
           cellnote= Mutationfrequency22 ,
           notecex=1,
           notecol="white",
           na.color=par("bg")            
)
dev.off() 


pdf( file=paste(outDir_g, "5_Mutation_Frequency.heatmap2.pdf", sep="/"),  height=5,   width=60 )
heatmap.2(x=Mutationfrequency22,          
          # dendrogram control
          dendrogram = "none", 
          Rowv = FALSE,
          Colv="Rowv" , 
          symm = FALSE,          
          # data scaling
          scale =  "none" ,
          na.rm=TRUE,         
          # colors
          col=my_col1(10),
          trace = "none", 
          key = F,
          keysize = 0.01,          
          # cell labeling
          cellnote= Mutationfrequency22 ,
          notecex=1,
          notecol="white",
          na.color=par("bg")        
)
dev.off() 



mutationFre = 1-Mutationfrequency2[1,]

mutationFre2 = sort(mutationFre, decreasing = T) 
mutationFre3 = round(x=mutationFre2*100, digits = 0)
mynames3 = names(mutationFre3)

pdf(paste(outDir_g, "6_Mutation_Frequency.pdf", sep="/"),  height=5,   width=20)
xx <- barplot(mutationFre3[1:80], xaxt = 'n',  xlab = '5-mer motifs',  ylim = c(0,110),  ylab = " Mutation Frequency (%)", col="red",  width = 1 )
## Add text at top of bars
text(x = xx, y = mutationFre3[1:80], label = mutationFre3[1:80], pos = 3, cex = 0.8, col = "black")
## Add x-axis labels 
axis(1, at=xx, labels=mynames3[1:80], tick=FALSE, las=2, line=0, cex.axis=0.5, cex=3 )
dev.off()



pdf(paste(outDir_g, "7_Mutation_Frequency.pdf", sep="/"),  height=5,   width=50)
xx <- barplot(mutationFre3, xaxt = 'n',  xlab = '5-mer motifs',  ylim = c(0,110),  ylab = " Mutation Frequency (%)", col="red",  width = 1 )
## Add text at top of bars
text(x = xx, y = mutationFre3, label = mutationFre3, pos = 3, cex = 0.8, col = "black")
## Add x-axis labels 
axis(1, at=xx, labels=mynames3, tick=FALSE, las=2, line=0, cex.axis=0.5, cex=3 )
dev.off()
























