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
  matrix2 = myScaleMatrix2( matrix2 )  
  pdf( file = paste(path2, fileName2, sep="/"),  width=width2, height=height2  )
      print( corrplot(matrix2, method = "circle", type = "full",  title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "circle", type = "upper", title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )                                       
      print( corrplot(matrix2, method = "number", type = "upper", title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "pie",    type = "upper", title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "color",  type = "upper", title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "color",  type = "full",  title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
  dev.off()
  
}  



MyHeatmaps_1A_f <- function(matrix2,  path2,   fileName2,   height2=30,   width2=30, is.corr2=TRUE,  my_col2 ) {                                                             
  ## matrix2 = myScaleMatrix2( matrix2 )  
  pdf( file = paste(path2, fileName2, sep="/"),  width=width2, height=height2  )
      print( corrplot(matrix2, method = "circle", type = "full",  title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "circle", type = "upper", title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", tl.srt = 45, col = my_col2 )  )                                       
      print( corrplot(matrix2, method = "number", type = "full", title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "pie",    type = "full", title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "color",  type = "full", title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "color",  type = "full",  title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
  dev.off()
  
}  



MyHeatmaps_2_f <- function(matrix2,  path2,   fileName2,   height2=30,   width2=30 ) { 
  matrix2 = myScaleMatrix2( matrix2 )  
  pdf( file = paste(path2, fileName2, sep="/"),  width=width2, height=height2  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("blue", "white", "red"))(20),       symm = FALSE,  scale =  "none"     )  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("green4", "white", "purple"))(20),  symm = FALSE,  scale =  "none"     )  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("cyan", "white", "red"))(20),       symm = FALSE,  scale =  "none"     )  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("cyan", "white", "purple"))(20),    symm = FALSE,  scale =  "none"     )  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("blue", "white", "purple"))(20),    symm = FALSE,  scale =  "none"     )  )
  dev.off()
  
}  


MyHeatmaps_3_f <- function(matrix2,  path2,   fileName2,   height2=30,   width2=30 ) { 
  matrix2 = myScaleMatrix2( matrix2 )  
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



#################################################
rawMatrix_1 <- read.table("Intervene_pairwise_frac_matrix.txt", header=TRUE, sep="\t", quote = "", comment.char = "") 
rawMatrix_1 = rawMatrix_1[,-1]
rownames(rawMatrix_1)
colnames(rawMatrix_1)
dim(rawMatrix_1)


outPath = "figures"
if( ! file.exists(outPath)  ) { dir.create(outPath,  recursive = TRUE)  }



rawMatrix_2 <- rawMatrix_1
dim(rawMatrix_2)


for(i in c(1:nrow(rawMatrix_2))) { 
     rawMatrix_2[i,i] = 0
}



my_col1=colorRampPalette( c(   "black", "black", "black", "pink", "red", "red",   "red4") )


MyHeatmaps_1A_f(matrix2=as.matrix(rawMatrix_2),  path2=outPath,   fileName2="heatmap.pdf",   
                height2=9,   width2=6, is.corr2=FALSE,   my_col2=my_col1(10)  )

 



rawMatrix_2 = rawMatrix_2*100
rawMatrix_3 = round(rawMatrix_2, digits = 0)
rawMatrix_3


my_col2=colorRampPalette( c( "grey80",  "yellow3", "red"  ), bias = 1.5 )


pdf( file="heatmap-3.pdf",  height=5,   width=5 )
heatmap.2 (x=as.matrix(rawMatrix_3),
           
           # dendrogram control
           dendrogram = "none", 
           Rowv = FALSE,
           Colv="Rowv" , 
           symm = FALSE,
           
           # data scaling
           scale =  "none" ,
           na.rm=TRUE,
           
           # colors
           col=my_col2(40),
           trace = "none", 
           
           # cell labeling
           cellnote= as.matrix(rawMatrix_3) ,
           notecex=1.0,
           notecol="white",
           na.color=par("bg") 
            
)
dev.off() 







