#! /usr/bin/env Rscript
###########################################################################################################
### bases_frequency.R: Statistic absolute frequency and relative frequency from the result of script Four_spikeins.pl.
### Author: Yong Peng
### Run "./bases_frequency.R  -h"  or "Rscript bases_frequency.R  -h" to get some help.
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
		description="bases_frequency.R: Statistic absolute frequency and relative frequency from the result of script Four_spikeins.pl. Junuary 2, 2020.",                             
		epilogue="For comments, bug reports etc..., please contact Yong Peng <yongp@outlook.com>"
  )

  opt_Local <- optparse::parse_args(parser_Local, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options
  return(opt_Local)   
}


opt_g = getParameters_f()
inFile_g     <- opt_g$file
outDir_g     <- opt_g$outDir
# inFile_g     <- "spikein_1_only-spikein.txt"
# outDir_g     <- "spikein_1_only-spikein"
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
vector1 = unlist( apply(X=matrix1, MARGIN=1, FUN=strsplit, split="" ) )
length(vector1)
nrow(matrix1) * 33
matrix2 = matrix( data=vector1, nrow = nrow(matrix1), byrow = TRUE )
dim(matrix2)

#table(matrix2[,1])
#table(matrix2[,2])
#table(matrix2[,3])
#table(matrix2[,22])
#length(which(matrix2[,22]=='T'))
#sum(matrix2[,22]=='T')

absoluteFrequency = matrix(nrow=4, ncol=ncol(matrix2))
rownames(absoluteFrequency) = c("A", "G", "C", "T")
relativeFrequency = absoluteFrequency
#absoluteFrequency 
#relativeFrequency

ncol(matrix2)
nrow(matrix2)

for(i in c(1:ncol(matrix2)) ) {  
  absoluteFrequency[1,i] = sum(matrix2[,i]=='A')
  absoluteFrequency[2,i] = sum(matrix2[,i]=='G')
  absoluteFrequency[3,i] = sum(matrix2[,i]=='C')
  absoluteFrequency[4,i] = sum(matrix2[,i]=='T')
  
  relativeFrequency[1,i] = absoluteFrequency[1,i]/nrow(matrix2)
  relativeFrequency[2,i] = absoluteFrequency[2,i]/nrow(matrix2)
  relativeFrequency[3,i] = absoluteFrequency[3,i]/nrow(matrix2)
  relativeFrequency[4,i] = absoluteFrequency[4,i]/nrow(matrix2)
  
}


 

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

relativeFrequency2 = relativeFrequency*100
relativeFrequency2 = round(x=relativeFrequency2, digits = 2)

pdf( file=paste(outDir_g, "heatmap-2.pdf", sep="/"),  height=5,   width=10 )
heatmap.2 (x=relativeFrequency2[,c(20:24)],
           
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
           cellnote= relativeFrequency2[,c(20:24)] ,
           notecex=2,
           notecol="white",
           na.color=par("bg") 
           
)
dev.off() 


relativeFrequency2 = round(x=relativeFrequency2, digits = 0)

pdf( file=paste(outDir_g, "heatmap-1.pdf", sep="/"),  height=5,   width=35 )
heatmap.2 (x=relativeFrequency2,
           
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
           cellnote= relativeFrequency2 ,
           notecex=2,
           notecol="white",
           na.color=par("bg") 
           
)
dev.off() 











