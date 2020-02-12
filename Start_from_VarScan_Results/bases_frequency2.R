#! /usr/bin/env Rscript
###########################################################################################################
### bases_frequency2.R: Statistic absolute frequency and relative frequency from the result of script Four_spikeins.pl.
### Author: Yong Peng
### Run "./bases_frequency2.R  -h"  or "Rscript bases_frequency2.R  -h" to get some help.
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
		description="bases_frequency2.R: Statistic absolute frequency2 and relative frequency2 from the result of script Four_spikeins.pl. Junuary 2, 2020.",                             
		epilogue="For comments, bug reports etc..., please contact Yong Peng <yongp@outlook.com>"
  )

  opt_Local <- optparse::parse_args(parser_Local, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options
  return(opt_Local)   
}


opt_g = getParameters_f()
inFile_g     <- opt_g$file
outDir_g     <- opt_g$outDir
# inFile_g     <- "spikein_1_only-spikein.txt"
# outDir_g     <- "spikein_1_only-spikein.5mers"
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g,   recursive = TRUE) }

options(digits=10)
continue_on_error <- function() {
  print( "NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())' " )
}
options( error=continue_on_error )  ## This option is very important.
#############################################





#############################################
matrix1 = read.table(file= inFile_g, header = FALSE, sep = character(0), quote = "\"'", dec = ".")
dim(matrix1)

my_4mer = c()
length(my_4mer)
my_4mer_name = c()
length(my_4mer_name)
bases4 = c("A", "G", "C", "T")
index1 = 0
for(i1 in bases4) {
  for(i2 in bases4) { 
    for(i3 in bases4) { 
      for(i4 in bases4) {
        index1 = index1 + 1
        my_4mer[ index1 ] = paste(i1, i2, "[AGCT]", i3, i4, sep="")
        my_4mer_name[ index1 ] = paste(i1, i2, "N", i3, i4, sep="")
      }
    }
  }
}
my_4mer
length(my_4mer)
my_4mer_name 
length(my_4mer_name)


library(tidyverse)
library(stringr)

vector1 = str_sub(string=matrix1[,1], start = 20, end = 24)
length(vector1)
my_list = vector(mode = "list", length = length(my_4mer_name) )
length(my_list)

for(i in c(1:length(my_4mer_name)) ) {
  my_list[[i]] = vector1[str_detect(string = vector1, pattern = my_4mer[i])]
}


absolutefrequency2 = matrix(nrow=4, ncol=length(my_4mer))
rownames(absolutefrequency2) = c("A", "G", "C", "T")
colnames(absolutefrequency2) = my_4mer_name
relativefrequency2 = absolutefrequency2
num_eachKmer = relativefrequency2[1,]

for(i in c(1:length(my_4mer_name)) ) {
  temp1 = my_list[[i]]
  num_eachKmer[i] = length(temp1)
  
  absolutefrequency2[1,i] = sum(str_sub(string=temp1, start = 3, end = 3) =='A')
  absolutefrequency2[2,i] = sum(str_sub(string=temp1, start = 3, end = 3) =='G')
  absolutefrequency2[3,i] = sum(str_sub(string=temp1, start = 3, end = 3) =='C')
  absolutefrequency2[4,i] = sum(str_sub(string=temp1, start = 3, end = 3) =='T')
  
  relativefrequency2[1,i] = absolutefrequency2[1,i]/num_eachKmer[i]
  relativefrequency2[2,i] = absolutefrequency2[2,i]/num_eachKmer[i]
  relativefrequency2[3,i] = absolutefrequency2[3,i]/num_eachKmer[i]
  relativefrequency2[4,i] = absolutefrequency2[4,i]/num_eachKmer[i]
}




write.table(x=absolutefrequency2, file = paste(outDir_g, "absolute_Frequency.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

write.table(x=relativefrequency2, file = paste(outDir_g, "relative_Frequency.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

write.table(x=num_eachKmer,       file = paste(outDir_g, "num_eachKmer.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")





suppressPackageStartupMessages( library(ComplexHeatmap) ) 
suppressPackageStartupMessages( library(corrplot) )
library(gplots)

reset_outliers2 <- function(x, na.rm = TRUE ) {
  qnt <- quantile(x, probs=c(0.1, 0.9) , type=1,  na.rm = na.rm )  
  y <- x
  y[x < qnt[1] ] <- qnt[1]
  y[x > qnt[2] ] <- qnt[2]    
  y
}

myScaleMatrix2 <- function( matrix_temp8, upper_temp8 = 1, lower_temp8 = -1 ) {
  rawMatrix_2 = reset_outliers2(matrix_temp8)  
  rawMatrix_2 = lower_temp8 + (upper_temp8 - lower_temp8) * ( rawMatrix_2 - min(rawMatrix_2) )/( max(rawMatrix_2)- min(rawMatrix_2) )
  return(rawMatrix_2)
}

MyHeatmaps_1_f <- function(matrix2,  path2,   fileName2,   height2=30,   width2=30, is.corr2=TRUE,  my_col2 ) {                                                             
  #matrix2 = myScaleMatrix2( matrix2 )  
  pdf( file = paste(path2, fileName2, sep="/"),  width=width2, height=height2  )
  print( corrplot(matrix2, method = "circle", type = "full",  title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
  print( corrplot(matrix2, method = "circle", type = "upper", title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )                                       
  print( corrplot(matrix2, method = "number", type = "upper", title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
  print( corrplot(matrix2, method = "pie",    type = "upper", title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
  print( corrplot(matrix2, method = "color",  type = "upper", title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
  print( corrplot(matrix2, method = "color",  type = "full",  title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
  dev.off()
  
}  

MyHeatmaps_2_f <- function(matrix2,  path2,   fileName2,   height2=30,   width2=30 ) { 
  #matrix2 = myScaleMatrix2( matrix2 )  
  pdf( file = paste(path2, fileName2, sep="/"),  width=width2, height=height2  )
  print( heatmap(x = matrix2,  col = colorRampPalette(c("blue", "white", "red"))(20),       symm = FALSE,  scale =  "none"  , Rowv = NA, Colv = NA     )  )
  print( heatmap(x = matrix2,  col = colorRampPalette(c("green4", "white", "purple"))(20),  symm = FALSE,  scale =  "none"  , Rowv = NA, Colv = NA     )  )
  print( heatmap(x = matrix2,  col = colorRampPalette(c("cyan", "white", "red"))(20),       symm = FALSE,  scale =  "none"  , Rowv = NA, Colv = NA     )  )
  print( heatmap(x = matrix2,  col = colorRampPalette(c("cyan", "white", "purple"))(20),    symm = FALSE,  scale =  "none"  , Rowv = NA, Colv = NA     )  )
  print( heatmap(x = matrix2,  col = colorRampPalette(c("blue", "white", "purple"))(20),    symm = FALSE,  scale =  "none"  , Rowv = NA, Colv = NA     )  )
  dev.off()
  
}  

MyHeatmaps_3_f <- function(matrix2,  path2,   fileName2,   height2=30,   width2=30 ) { 
  #matrix2 = myScaleMatrix2( matrix2 )  
  pdf( file = paste(path2, fileName2, sep="/"),  width=width2, height=height2  )
  print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("blue", "white", "red"))(20)        , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )          
  print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("blue", "white", "red"))(20)        , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
  print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("green4", "white", "purple"))(20)   , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
  print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("green4", "white", "purple"))(20)   , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
  print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("cyan", "white", "red"))(20)        , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
  print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("cyan", "white", "red"))(20)        , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
  print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("cyan", "white", "purple"))(20)     , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
  print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("cyan", "white", "purple"))(20)     , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
  print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("blue", "white", "purple"))(20)     , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
  print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("blue", "white", "purple"))(20)     , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
  dev.off()
  
}  


my_col1=colorRampPalette( c( "blue", "yellow3" , "pink", "red", "red"  ), bias = 1 )

relativefrequency22 = relativefrequency2*100
relativefrequency22 = round(x=relativefrequency22, digits = 1)

pdf( file=paste(outDir_g, "relative_Frequency.heatmap.pdf", sep="/"),  height=5,   width=60 )
heatmap.2(x=relativefrequency22,
           
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
           cellnote= relativefrequency22 ,
           notecex=1,
           notecol="white",
           na.color=par("bg") 
           
)
dev.off() 





relativefrequency22 = round(x=relativefrequency22, digits = 0)
pdf( file=paste(outDir_g, "relative_Frequency.heatmap2.pdf", sep="/"),  height=5,   width=60 )
heatmap.2(x=relativefrequency22,
          
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
          cellnote= relativefrequency22 ,
          notecex=1,
          notecol="white",
          na.color=par("bg") 
          
)
dev.off() 





