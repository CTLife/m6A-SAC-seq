
 
#################################################

outDir_g = "histogram"
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g,   recursive = TRUE) }

options(digits=10)
continue_on_error <- function() {
  print( "NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())' " )
}
options( error=continue_on_error )  ## This option is very important.
#############################################






#############################################
matrix1 = read.table(file= "number.txt", header = F, sep = character(0), quote = "\"'", dec = ".")
dim(matrix1) 
 
  
ratio1 = as.numeric( matrix1[,1] )
col_names1 = matrix1[,2]
  
  
pdf(paste(outDir_g, "number_of_peaks.pdf", sep="/"),  height=5,   width=10)

xx1 <- barplot( ratio1, xaxt = 'n',  xlab = 'samples',   ylab = "number", col="red",  width = 1 )
text(x = xx1, y = ratio1, label = ratio1, pos = 3, cex = 0.9, col = "black")
axis(1, at=xx1, labels=col_names1, tick=FALSE, las=2, line=0, cex.axis=0.6, cex=5 )

dev.off()

  
 
 