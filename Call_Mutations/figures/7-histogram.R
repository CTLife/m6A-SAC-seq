
##############################################################################################################################################################################################
suppressPackageStartupMessages( library(optparse) )  ## To run the script in command lines.

getParameters_f <- function() {
	option_list_Local <- list(   ## Options list with associated default value.  
			optparse::make_option(opt_str=c("-i", "--inputFile"),
	                        default="7-getMotifs/Ig1.bed.number",
	                        type="character",   dest="inputFile",
	                        help="Input file. [default: %default]."),
	  
			optparse::make_option(opt_str=c("-o", "--outFile"),
	                        default="7-getMotifs/Ig1.bed.pdf",
	                        type="character",   dest="outFile",
	                        help="Output file. [default: %default].")	
   )
	
	## Now parse the command line to check which option is given and get associated values.
	parser_Local <- optparse::OptionParser(usage="usage: Rscript %prog [options]",
			option_list=option_list_Local, 
			description="May 22, 2021.",                             
			epilogue="For comments, bug reports etc..., please contact Yong Peng <yongp@outlook.com>."
	)
	opt_Local <- optparse::parse_args(parser_Local, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options
	return(opt_Local)
}
##############################################################################################################################################################################################





##############################################################################################################################################################################################
opt_g = getParameters_f()  

inputFile_g    <- as.character(opt_g$inputFile)
outFile_g      <- as.character(opt_g$outFile)

 

#############################################
matrix1 = read.table(file= inputFile_g, header = F, sep = character(0), quote = "\"'", dec = ".")
dim(matrix1) 
 
  
ratio1 = as.numeric( matrix1[,1] )
col_names1 = matrix1[,2]
  
  
pdf( outFile_g ,  height=5,   width=10)

xx1 <- barplot( ratio1, xaxt = 'n',  xlab = 'motifs',   ylab = "number", col="red",  width = 1 )
text(x = xx1, y = ratio1, label = ratio1, pos = 3, cex = 0.9, col = "black")
axis(1, at=xx1, labels=col_names1, tick=FALSE, las=2, line=0, cex.axis=0.6, cex=5 )

dev.off()

  
 
 
