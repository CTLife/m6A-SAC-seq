
## To run the script in command lines.
#################################################

outDir_g = "5-histogram"
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g,   recursive = TRUE) }

options(digits=10)
continue_on_error <- function() {
  print( "NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())' " )
}
options( error=continue_on_error )  ## This option is very important.
#############################################






#############################################
matrix1 = read.table(file= "2-probesFrequency/spikein-1.seq/1_Absolute_Frequency.txt", header = TRUE, sep = character(0), quote = "\"'", dec = ".")
matrix2 = read.table(file= "2-probesFrequency/spikein-2.seq/1_Absolute_Frequency.txt", header = TRUE, sep = character(0), quote = "\"'", dec = ".")
matrix3 = read.table(file= "2-probesFrequency/spikein-3.seq/1_Absolute_Frequency.txt", header = TRUE, sep = character(0), quote = "\"'", dec = ".")
matrix4 = read.table(file= "2-probesFrequency/spikein-4.seq/1_Absolute_Frequency.txt", header = TRUE, sep = character(0), quote = "\"'", dec = ".")
matrix5 = read.table(file= "2-probesFrequency/spikein-5.seq/1_Absolute_Frequency.txt", header = TRUE, sep = character(0), quote = "\"'", dec = ".")

cat("##########################\n")
dim(matrix1)
dim(matrix2)
dim(matrix3)
dim(matrix4)
dim(matrix5)
cat("##########################\n")

col_names1 = colnames(matrix1)
col_names2 = colnames(matrix2)
col_names3 = colnames(matrix3)
col_names4 = colnames(matrix4)
col_names5 = colnames(matrix5)

# cat("##########################\n")
# col_names1 == col_names2
# col_names1 == col_names3
# col_names1 == col_names4
# col_names1 == col_names5
# cat("##########################\n")

 

numbers1 = as.vector( colSums(matrix1) )
numbers2 = as.vector( colSums(matrix2) )
numbers3 = as.vector( colSums(matrix3) )
numbers4 = as.vector( colSums(matrix4) )
numbers5 = as.vector( colSums(matrix5) )

 

ratio1 = 100*as.vector( unlist(matrix1[1,]) )/numbers1
ratio2 = 100*as.vector( unlist(matrix2[1,]) )/numbers2
ratio3 = 100*as.vector( unlist(matrix3[1,]) )/numbers3
ratio4 = 100*as.vector( unlist(matrix4[1,]) )/numbers4
ratio5 = 100*as.vector( unlist(matrix5[1,]) )/numbers5


ratio1 = 100-ratio1
ratio2 = 100-ratio2
ratio3 = 100-ratio3
ratio4 = 100-ratio4
ratio5 = 100-ratio5

ratio1 = round(x=ratio1, digits = 0) 
ratio2 = round(x=ratio2, digits = 0) 
ratio3 = round(x=ratio3, digits = 0) 
ratio4 = round(x=ratio4, digits = 0) 
ratio5 = round(x=ratio5, digits = 0) 





pdf(paste(outDir_g, "all_motifs.pdf", sep="/"), height=20, width=40 )
par(mfrow = c(5, 1) )

xx1 <- barplot( ratio1, xaxt = 'n',  xlab = '5-mer motifs (0% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx1, y = ratio1, label = ratio1, pos = 3, cex = 1, col = "black")
axis(1, at=xx1, labels=col_names1, tick=FALSE, las=2, line=0, cex.axis=1, cex = 9 )

xx2 <- barplot( ratio2, xaxt = 'n',  xlab = '5-mer motifs (25% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx2, y = ratio2, label = ratio2, pos = 3, cex = 1, col = "black")
axis(1, at=xx2, labels=col_names2, tick=FALSE, las=2, line=0, cex.axis=1, cex = 9 )

xx3 <- barplot( ratio3, xaxt = 'n',  xlab = '5-mer motifs (50% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx3, y = ratio3, label = ratio3, pos = 3, cex = 1, col = "black")
axis(1, at=xx3, labels=col_names3, tick=FALSE, las=2, line=0, cex.axis=1, cex = 9 )

xx4 <- barplot( ratio4, xaxt = 'n',  xlab = '5-mer motifs (75% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx4, y = ratio4, label = ratio4, pos = 3, cex = 1, col = "black")
axis(1, at=xx4, labels=col_names4, tick=FALSE, las=2, line=0, cex.axis=1, cex = 9 )

xx5 <- barplot( ratio5, xaxt = 'n',  xlab = '5-mer motifs (100% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx5, y = ratio5, label = ratio5, pos = 3, cex = 1, col = "black")
axis(1, at=xx5, labels=col_names5, tick=FALSE, las=2, line=0, cex.axis=1, cex = 9 )

dev.off()
 
  
 


DRACH_index = grepl(pattern = "[AGT][AG]AC[ACT]", x=col_names1, ignore.case = FALSE, perl = TRUE,  fixed = FALSE, useBytes = FALSE)  

matrix1 = matrix1[, DRACH_index]
matrix2 = matrix2[, DRACH_index]
matrix3 = matrix3[, DRACH_index]
matrix4 = matrix4[, DRACH_index]
matrix5 = matrix5[, DRACH_index]

col_names1 = colnames(matrix1)
col_names2 = colnames(matrix2)
col_names3 = colnames(matrix3)
col_names4 = colnames(matrix4)
col_names5 = colnames(matrix5)

numbers1 = as.vector( colSums(matrix1) )
numbers2 = as.vector( colSums(matrix2) )
numbers3 = as.vector( colSums(matrix3) )
numbers4 = as.vector( colSums(matrix4) )
numbers5 = as.vector( colSums(matrix5) )

ratio1 = 100*as.vector( unlist(matrix1[1,]) )/numbers1
ratio2 = 100*as.vector( unlist(matrix2[1,]) )/numbers2
ratio3 = 100*as.vector( unlist(matrix3[1,]) )/numbers3
ratio4 = 100*as.vector( unlist(matrix4[1,]) )/numbers4
ratio5 = 100*as.vector( unlist(matrix5[1,]) )/numbers5

ratio1 = 100-ratio1
ratio2 = 100-ratio2
ratio3 = 100-ratio3
ratio4 = 100-ratio4
ratio5 = 100-ratio5

ratio1 = round(x=ratio1, digits = 0) 
ratio2 = round(x=ratio2, digits = 0) 
ratio3 = round(x=ratio3, digits = 0) 
ratio4 = round(x=ratio4, digits = 0) 
ratio5 = round(x=ratio5, digits = 0) 


pdf(paste(outDir_g, "DRACH.pdf", sep="/") , height=20, width=20)
par(mfrow = c(5, 1) )

xx1 <- barplot( ratio1, xaxt = 'n',  xlab = 'DRACH motifs (0% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx1, y = ratio1, label = ratio1, pos = 3, cex = 3, col = "black")
axis(1, at=xx1, labels=col_names1, tick=FALSE, las=2, line=0, cex.axis=2, cex = 20 )

xx2 <- barplot( ratio2, xaxt = 'n',  xlab = 'DRACH motifs (25% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx2, y = ratio2, label = ratio2, pos = 3, cex = 3, col = "black")
axis(1, at=xx2, labels=col_names2, tick=FALSE, las=2, line=0, cex.axis=2, cex = 20 )

xx3 <- barplot( ratio3, xaxt = 'n',  xlab = 'DRACH motifs (50% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx3, y = ratio3, label = ratio3, pos = 3, cex = 3, col = "black")
axis(1, at=xx3, labels=col_names3, tick=FALSE, las=2, line=0, cex.axis=2, cex = 20 )

xx4 <- barplot( ratio4, xaxt = 'n',  xlab = 'DRACH motifs (75% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx4, y = ratio4, label = ratio4, pos = 3, cex = 3, col = "black")
axis(1, at=xx4, labels=col_names4, tick=FALSE, las=2, line=0, cex.axis=2, cex = 20 )

xx5 <- barplot( ratio5, xaxt = 'n',  xlab = 'DRACH motifs (100% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx5, y = ratio5, label = ratio5, pos = 3, cex = 3, col = "black")
axis(1, at=xx5, labels=col_names5, tick=FALSE, las=2, line=0, cex.axis=2, cex = 20 )

dev.off()
 
 









matrix1 = read.table(file= "2-probesFrequency/spikein-1.seq/1_Absolute_Frequency.txt", header = TRUE, sep = character(0), quote = "\"'", dec = ".")
matrix2 = read.table(file= "2-probesFrequency/spikein-2.seq/1_Absolute_Frequency.txt", header = TRUE, sep = character(0), quote = "\"'", dec = ".")
matrix3 = read.table(file= "2-probesFrequency/spikein-3.seq/1_Absolute_Frequency.txt", header = TRUE, sep = character(0), quote = "\"'", dec = ".")
matrix4 = read.table(file= "2-probesFrequency/spikein-4.seq/1_Absolute_Frequency.txt", header = TRUE, sep = character(0), quote = "\"'", dec = ".")
matrix5 = read.table(file= "2-probesFrequency/spikein-5.seq/1_Absolute_Frequency.txt", header = TRUE, sep = character(0), quote = "\"'", dec = ".")


DRACH_index = grepl(pattern = "[AGT][AG]AC[ACT]", x=col_names1, ignore.case = FALSE, perl = TRUE,  fixed = FALSE, useBytes = FALSE)  

matrix1 = matrix1[, -DRACH_index]
matrix2 = matrix2[, -DRACH_index]
matrix3 = matrix3[, -DRACH_index]
matrix4 = matrix4[, -DRACH_index]
matrix5 = matrix5[, -DRACH_index]

col_names1 = colnames(matrix1)
col_names2 = colnames(matrix2)
col_names3 = colnames(matrix3)
col_names4 = colnames(matrix4)
col_names5 = colnames(matrix5)

numbers1 = as.vector( colSums(matrix1) )
numbers2 = as.vector( colSums(matrix2) )
numbers3 = as.vector( colSums(matrix3) )
numbers4 = as.vector( colSums(matrix4) )
numbers5 = as.vector( colSums(matrix5) )

ratio1 = 100*as.vector( unlist(matrix1[1,]) )/numbers1
ratio2 = 100*as.vector( unlist(matrix2[1,]) )/numbers2
ratio3 = 100*as.vector( unlist(matrix3[1,]) )/numbers3
ratio4 = 100*as.vector( unlist(matrix4[1,]) )/numbers4
ratio5 = 100*as.vector( unlist(matrix5[1,]) )/numbers5

ratio1 = 100-ratio1
ratio2 = 100-ratio2
ratio3 = 100-ratio3
ratio4 = 100-ratio4
ratio5 = 100-ratio5

ratio1 = round(x=ratio1, digits = 0) 
ratio2 = round(x=ratio2, digits = 0) 
ratio3 = round(x=ratio3, digits = 0) 
ratio4 = round(x=ratio4, digits = 0) 
ratio5 = round(x=ratio5, digits = 0) 

 




pdf(paste(outDir_g, "nonDRACH.pdf", sep="/"), height=20, width=40 )
par(mfrow = c(5, 1) )

xx1 <- barplot( ratio1, xaxt = 'n',  xlab = '5-mer motifs (0% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx1, y = ratio1, label = ratio1, pos = 3, cex = 1, col = "black")
axis(1, at=xx1, labels=col_names1, tick=FALSE, las=2, line=0, cex.axis=1, cex = 9 )

xx2 <- barplot( ratio2, xaxt = 'n',  xlab = '5-mer motifs (25% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx2, y = ratio2, label = ratio2, pos = 3, cex = 1, col = "black")
axis(1, at=xx2, labels=col_names2, tick=FALSE, las=2, line=0, cex.axis=1, cex = 9 )

xx3 <- barplot( ratio3, xaxt = 'n',  xlab = '5-mer motifs (50% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx3, y = ratio3, label = ratio3, pos = 3, cex = 1, col = "black")
axis(1, at=xx3, labels=col_names3, tick=FALSE, las=2, line=0, cex.axis=1, cex = 9 )

xx4 <- barplot( ratio4, xaxt = 'n',  xlab = '5-mer motifs (75% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx4, y = ratio4, label = ratio4, pos = 3, cex = 1, col = "black")
axis(1, at=xx4, labels=col_names4, tick=FALSE, las=2, line=0, cex.axis=1, cex = 9 )

xx5 <- barplot( ratio5, xaxt = 'n',  xlab = '5-mer motifs (100% m6A)',  ylim = c(0,100),  ylab = "Mut(%)", col="red",  width = 1 )
text(x = xx5, y = ratio5, label = ratio5, pos = 3, cex = 1, col = "black")
axis(1, at=xx5, labels=col_names5, tick=FALSE, las=2, line=0, cex.axis=1, cex = 9 )

dev.off()
 
  






