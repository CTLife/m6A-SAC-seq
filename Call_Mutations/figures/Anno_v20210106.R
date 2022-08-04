##############################################################################################################################################################################################
## Author: Yong Peng, yongp@outlook.com.
## Run " Rscript Anno_v20210106.R -h" to get some help.
##############################################################################################################################################################################################





##############################################################################################################################################################################################
suppressPackageStartupMessages( library(optparse) )  ## To run the script in command lines.

getParameters_f <- function() {
	option_list_Local <- list(   ## Options list with associated default value.  
			optparse::make_option(opt_str=c("-i", "--inputDir"),
	                        default="1-BED",
	                        type="character",   dest="inputDir",
	                        help="Path to the directory containing BED files with genomic regions. The 1st to 6th columns of BED files, all of them are required.  [default: %default]."),
	  
			optparse::make_option(opt_str=c("-o", "--outDir"),
	                        default="2-Annotation",
	                        type="character",   dest="outDir",
	                        help="Path to the directory containing all the analysis results. [default: %default]."),
	  
			optparse::make_option(opt_str=c("-r", "--refGenome"),
	                        default="hg38",
	                        type="character",   dest="refGenome",
	                        help="A string that defines the genome assembly or reference genome for your input BED files, such as hg38, hg19, mm10 and mm9. [default: %default]."),
	  
			optparse::make_option(opt_str=c("-s", "--stranded"),
	                        default=FALSE,
	                        type="logical",   dest="stranded",
	                        help="Logical, whether find nearest/overlap gene in the same strand. [default: %default]."),
		  
			optparse::make_option(opt_str=c("-f", "--fast"),
	                        default=FALSE,
	                        type="logical",   dest="fast",
	                        help="Logical, whether use fast mode. If it is TRUE, some steps will be skipped. [default: %default]."),    
	  
			optparse::make_option(opt_str=c("-U", "--upstream"),
			                      default=3000,
			                      type="integer",   dest="upstream",
			                      help="Number of upstream bases of TSS to define promoter. [default: %default]."),
			
			optparse::make_option(opt_str=c("-D", "--downstream"),
			                      default=3000,
			                      type="integer",   dest="downstream",
			                      help="Number of downstream bases of TSS to define promoter. [default: %default].")
			
   )
	
	## Now parse the command line to check which option is given and get associated values.
	parser_Local <- optparse::OptionParser(usage="usage: Rscript %prog [options]",
			option_list=option_list_Local, 
			description="ChIPseeker v6, February 22, 2020.",                             
			epilogue="For comments, bug reports etc..., please contact Yong Peng <yongp@outlook.com>."
	)
	opt_Local <- optparse::parse_args(parser_Local, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options
	return(opt_Local)
}
##############################################################################################################################################################################################





##############################################################################################################################################################################################
opt_g = getParameters_f()  

inputDir_g    <- as.character(opt_g$inputDir)
outDir_g      <- as.character(opt_g$outDir)
refGenome_g   <- as.character(opt_g$refGenome)
stranded_g    <- as.logical(opt_g$stranded)
fast_g        <- as.logical(opt_g$fast)
upstream_g    <- as.integer(opt_g$upstream)
downstream_g  <- as.integer(opt_g$downstream)

rm(getParameters_f)  
rm(opt_g)  

options(digits=10)
continue_on_error_g <- function() {
    print( "NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'. " )
}
options( error=continue_on_error_g )  ## This option is very important.

#  inputDir_g    <- "1-BED"
#  outDir_g      <- "2-Annotation"
#  refGenome_g   <- "hg38"
#  stranded_g    <- TRUE
#  fast_g        <- FALSE
#  upstream_g    <- 1000
#  downstream_g  <- 0

if( ! file.exists(outDir_g)  ) { dir.create(outDir_g,  recursive = TRUE)  }
######################################################################################################################################################





######################################################################################################################################################
continue_on_error <- function() {
  print(" NOTE: THERE WAS AN ERROR HERE. Yong Peng. We are continuing because we have set 'options(error=continue_on_error())'. ")
}
# This is the key option
options(error=continue_on_error) 

suppressPackageStartupMessages( library(GenomicFeatures) ) 
suppressPackageStartupMessages( library(GenomicRanges) ) 
suppressPackageStartupMessages( library(clusterProfiler) ) 
suppressPackageStartupMessages( library(ReactomePA) ) 
suppressPackageStartupMessages( library(ChIPseeker) ) 
suppressPackageStartupMessages( library(DOSE) ) 
suppressPackageStartupMessages( library(ggplot2) ) 
suppressPackageStartupMessages( library(topGO) ) 
suppressPackageStartupMessages( library(KEGG.db) ) 
suppressPackageStartupMessages( library(enrichplot) ) 
suppressPackageStartupMessages( library(ggupset) ) 
suppressPackageStartupMessages( library(ggimage) ) 

my_txdb_g      <- ""
my_orgdb_g     <- ""
my_organism_g  <- ""
my_organism2_g <- ""
if(refGenome_g == "hg38") {
    suppressPackageStartupMessages( library(TxDb.Hsapiens.UCSC.hg38.knownGene) ) 
    suppressPackageStartupMessages( library(org.Hs.eg.db) )  
    my_txdb_g  <- TxDb.Hsapiens.UCSC.hg38.knownGene
    my_orgdb_g <- "org.Hs.eg.db"
    my_organism_g <- "human"     
    my_organism2_g <- "hsa"
}
if(refGenome_g == "mm10") {
    suppressPackageStartupMessages( library(TxDb.Mmusculus.UCSC.mm10.knownGene) ) 
    suppressPackageStartupMessages( library(org.Mm.eg.db) )  
    my_txdb_g  <- TxDb.Mmusculus.UCSC.mm10.knownGene
    my_orgdb_g <- "org.Mm.eg.db"
    my_organism_g <- "mouse"
    my_organism2_g <- "mmu"
}
if(refGenome_g == "danRer11") {
    suppressPackageStartupMessages( library(TxDb.Drerio.UCSC.danRer11.refGene) ) 
    suppressPackageStartupMessages( library(org.Dr.eg.db) )  
    my_txdb_g  <- TxDb.Drerio.UCSC.danRer11.refGene
    my_orgdb_g <- "org.Dr.eg.db"
    my_organism_g <- "zebrafish"
    my_organism2_g <- "dre"
}

sink( paste(outDir_g, "Parameters.log.txt", sep="/")   )
  print( "inputDir_g:"  )
  print(  inputDir_g  )
  cat("\n\n\n") 
  print( "outDir_g: ")  
  print(  outDir_g ) 
  cat("\n\n\n") 
  print( "refGenome_g: ")  
  print(  refGenome_g ) 
  cat("\n\n\n") 
  print( "stranded_g: ")  
  print(  stranded_g ) 
  cat("\n\n\n") 
  print( "fast_g: ")  
  print(  fast_g ) 
  cat("\n\n\n") 
  print( "upstream_g: ")  
  print(  upstream_g ) 
  cat("\n\n\n") 
  print( "downstream_g: ")  
  print(  downstream_g ) 
  cat("\n\n\n") 
  print( "my_organism_g:" )
  print( my_organism_g )
  cat("\n\n\n") 
  print( "my_organism2_g:" )
  print( my_organism2_g )
  cat("\n\n\n") 
sink()   


## read the input files.
peakFiles          <- list.files(path = inputDir_g, pattern = ".bed$",  full.names = TRUE  )
peakFiles_onlyName <- list.files(path = inputDir_g, pattern = ".bed$",  full.names = FALSE )
cat("####################\n\n\n")
print(peakFiles_onlyName)
print(length(peakFiles_onlyName))
cat("####################\n\n\n")

sink( paste(outDir_g, "All_the_input_BED_files.txt", sep="/")   )
print(peakFiles)
print(length(peakFiles))
cat("####################\n\n\n")
print(peakFiles_onlyName)
print(length(peakFiles_onlyName))
sink()   
######################################################################################################################################################
  
 



######################################################################################################################################################
MyTheme_1 <- function(textSize1=14, hjust1=NULL, vjust1=NULL,  angle1=NULL) {    # "hjust=1, vjust=1, angle=30" for some boxplots.
  theme(  
    line  = element_line(colour="black",  size=1.0,   linetype=1,      lineend=NULL),                                                                                        ## all line elements.          局部优先总体,下面3个也是,只对非局部设置有效.   所有线属性.
    rect  = element_rect(colour="black",  size=1.0,   linetype=1,      fill="transparent" ),                                                                                 ## all rectangluar elements.    hjust=1: 靠右对齐.   所有矩形区域属性.
    text  = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),    ## all text elements.           "serif" for a serif font. 所有文本相关属性.
    title = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),    ## all title elements: plot, axes, legends.    hjust:水平对齐的方向.  所有标题属性.
    ## aspect.ratio = 1,   ##高宽比
    
    axis.title    = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## label of axes (element_text; inherits from text).  horizontal: 水平的, 水平线 
    axis.title.x  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## x axis label (element_text; inherits from axis.title)
    axis.title.y  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=90,      lineheight=1.0,  margin = NULL, debug = NULL),       ## y axis label (element_text; inherits from axis.title)
    axis.text     = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## tick labels along axes (element_text; inherits from text). 坐标轴刻度的标签的属性.                                                         
    axis.text.x   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=hjust1, vjust=vjust1, angle=angle1,  lineheight=1.0,  margin = NULL, debug = NULL),       ## x axis tick labels (element_text; inherits from axis.text)
    axis.text.y   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## y axis tick labels (element_text; inherits from axis.text)
    
    axis.ticks        = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## tick marks along axes (element_line; inherits from line). 坐标轴刻度线.
    axis.ticks.x      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## x axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.y      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## y axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.length = grid::unit(2.0,   "mm",   data=NULL),                                      ## length of tick marks (unit), ‘"mm"’ Millimetres.  10 mm = 1 cm.  刻度线长度
    axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	 ## lines along axes (element_line; inherits from line). 坐标轴线
    axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	 ## line along x axis (element_line; inherits from axis.line)
    axis.line.y       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),	   ## line along y axis (element_line; inherits from axis.line)
    
    legend.background    = element_rect(colour="transparent", size=1, linetype=1, fill="transparent" ), 	      ## background of legend (element_rect; inherits from rect)
    legend.spacing       = grid::unit(1, "mm", data=NULL), 	                                                    ## extra space added around legend (unit). linetype=1指的是矩形边框的类型.
    legend.key           = element_rect(colour="transparent", size=2, linetype=1, fill="transparent" ), 	      ## background underneath legend keys. 图例符号. size=1指的是矩形边框的大小.
    legend.key.size      = grid::unit(6,   "mm", data=NULL) , 	                                                ## size of legend keys   (unit; inherits from legend.key.size)
    legend.key.height    = grid::unit(6.5, "mm", data=NULL) , 	                                                ## key background height (unit; inherits from legend.key.size)
    legend.key.width     = grid::unit(8,   "mm", data=NULL) ,                                                   ## key background width  (unit; inherits from legend.key.size)
    legend.text          = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	##legend item labels. 图例文字标签.
    legend.text.align    = 0, 	                    ## alignment of legend labels (number from 0 (left) to 1 (right))
    legend.title         = element_blank(),   	    ## title of legend (element_text; inherits from title)
    legend.title.align   = 0, 	                    ## alignment of legend title (number from 0 (left) to 1 (right))
    legend.position      = "right", 	              ## the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
    legend.direction     = "vertical",        	    ## layout of items in legends  ("horizontal" or "vertical")   图例排列方向
    legend.justification = "center",      	        ## anchor point for positioning legend inside plot ("center" or two-element numeric vector)  图例居中方式
    legend.box           = NULL, 	                  ## arrangement of multiple legends ("horizontal" or "vertical")  多图例的排列方式
    legend.box.just      = NULL, 	                  ## justification of each legend within the overall bounding box, when there are multiple legends ("top", "bottom", "left", or "right")  多图例的居中方式
    
    panel.background   = element_rect(colour="transparent", size=0.0, linetype=1, fill="transparent" ),   	## background of plotting area, drawn underneath plot (element_rect; inherits from rect)
    panel.border       = element_rect(colour="black", size=0.5, linetype=1, fill=NA ), 	                    ## border around plotting area, drawn on top of plot so that it covers tick marks and grid lines. This should be used with fill=NA (element_rect; inherits from rect)
    panel.spacing      = grid::unit(1, "mm", data=NULL) , 	                                                ## margin around facet panels (unit)  分面绘图区之间的边距
    panel.spacing.x    = grid::unit(1, "mm", data=NULL) ,
    panel.spacing.y    = grid::unit(1, "mm", data=NULL) ,
    panel.grid         = element_blank(), 	                                                                ## grid lines (element_line; inherits from line)  绘图区网格线
    panel.grid.major   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## major grid lines (element_line; inherits from panel.grid)  主网格线
    panel.grid.minor   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## minor grid lines (element_line; inherits from panel.grid)  次网格线
    panel.grid.major.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## vertical major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.major.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.minor.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## vertical minor grid lines (element_line; inherits from panel.grid.minor)
    panel.grid.minor.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal minor grid lines (element_line; inherits from panel.grid.minor)
    
    plot.background	= element_rect(colour="transparent", size=NULL, linetype=NULL, fill="transparent" ),                                                ## background of the entire plot (element_rect; inherits from rect)  整个图形的背景
    plot.title      = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=0.5, vjust=0.5,   angle=NULL, lineheight=NULL),     ## plot title (text appearance) (element_text; inherits from title)  图形标题
    plot.margin     = grid::unit(c(5, 5, 5, 5), "mm", data=NULL), 	                                                                                    ## margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
    
    strip.background = element_rect(colour=NULL,    size=NULL, linetype=NULL, fill=NULL ), 	                                                      ## background of facet labels (element_rect; inherits from rect)  分面标签背景
    strip.text       = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels (element_text; inherits from text)
    strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels along horizontal direction (element_text; inherits from strip.text)
    strip.text.y     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL)   	  ## facet labels along vertical direction (element_text; inherits from strip.text) 
  ) 
} 


MySaveGgplot2_1 <- function(path1, fileName1,  height1, width1) {
  SVG1 <- paste(path1,  "/",  "SVG",  sep = "",  collapse = NULL)
  PNG1 <- paste(path1,  "/",  "PNG",  sep = "",  collapse = NULL)
  PDF1 <- paste(path1,  "/",  "PDF",  sep = "",  collapse = NULL)
  EPS1 <- paste(path1,  "/",  "EPS",  sep = "",  collapse = NULL)
  if( ! file.exists(SVG1) ) { dir.create(SVG1) }
  if( ! file.exists(PNG1) ) { dir.create(PNG1) }
  if( ! file.exists(PDF1) ) { dir.create(PDF1) }
  if( ! file.exists(EPS1) ) { dir.create(EPS1) }
  ggsave( filename = paste(SVG1,  "/",  fileName1,  ".svg",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(PNG1,  "/",  fileName1,  ".png",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(PDF1,  "/",  fileName1,  ".pdf",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(EPS1,  "/",  fileName1,  ".eps",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200,   device=cairo_ps)         
}
######################################################################################################################################################





######################################################################################################################################################
#############  For annotating one file or one group.


##
MyPeaksAnno_OneGroup_1_g <- function(myPeak_t1, myPath_t1) { 
  myTempFunction1 <- function() {
      if( ! file.exists(myPath_t1)  ) { dir.create(myPath_t1,  recursive = TRUE)  }
      pdf(file=paste(myPath_t1, "1_GenomicRegions_over_Chromosomes.pdf", sep="/"),  width=20, height=10)
          print( covplot(peak=myPeak_t1, weightCol = "V5", xlab = "Chromosome Size (bp)", ylab = "strength", title = "Genomic Regions over Chromosomes", chrs = NULL, xlim = NULL, lower = 0.01) ) 
      dev.off()
      pdf(file=paste(myPath_t1, "2_GenomicRegions_over_only-chromXY.pdf", sep="/"),  width=20, height=5)
          print( covplot(peak=myPeak_t1, weightCol = "V5", xlab = "Chromosome Size (bp)", ylab = "strength", title = "Genomic Regions over Chromosome X and Y", chrs = c("chrX", "chrY"), xlim = NULL, lower = 0.01) )                                                           
      dev.off()
      pdf(file=paste(myPath_t1, "3_chromLength.pdf", sep="/"),  width=20, height=10) 
          print( covplot(peak=myPeak_t1, weightCol = "V5", xlab = "Chromosome Size (bp)", ylab = "strength", title = "Length of Each Chromosome", chrs = NULL, xlim = NULL, lower = -1) )                                                           
      dev.off()
  }
  tryCatch(
    myTempFunction1(),
    error = function(err){"MyPeaksAnno_OneGroup_1_g:error_1. Yong Peng."}
  )
}


##
MyPeaksSignal_OneGroup_1_g <- function(myPeak_t1, myPath_t1, up1, down1) { 
  myTempFunction2 <- function() {
      if( ! file.exists(myPath_t1)  ) { dir.create(myPath_t1,  recursive = TRUE)  }
      promoter2  <- getPromoters(TxDb=my_txdb_g, upstream=up1, downstream=down1)
      tagMatrix1 <- getTagMatrix(myPeak_t1, windows=promoter2)
      write.table(x=tagMatrix1, file =paste(myPath_t1, "1_peakMatrix_TSS.txt", sep="/"), 
                  append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )

      pdf(file = paste(myPath_t1, "2_Heatmap_of_GenomicRegions_binding_to_TSS_regions.pdf", sep="/"),  width=12, height=10)
          print( tagHeatmap(tagMatrix1, xlim=c(-up1, down1), color="red") )
      dev.off()

      pdf(file=paste(myPath_t1, "3_Average_Profile_of_GenomicRegions_binding_to_TSS_region.pdf", sep="/"),  width=12, height=10)
          print( plotAvgProf(tagMatrix1, xlim=c(-up1, down1),  xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency" )  )
      dev.off()
  }
  tryCatch(
    myTempFunction2(),
    error = function(err){"MyPeaksSignal_OneGroup_1_g:error_2."}
  )
}


##
MyPeaksDistribution_OneGroup_1_g <- function(myPeak_t1, myPath_t1, up1, down1) { 
  myTempFunction3 <- function() {    
      if( ! file.exists(myPath_t1)  ) { dir.create(myPath_t1,  recursive = TRUE)  }
      peakAnno <- annotatePeak( myPeak_t1, tssRegion=c(-up1, down1),  TxDb=my_txdb_g, annoDb=my_orgdb_g, sameStrand = stranded_g , overlap="all" )  
      sink(  file = paste(myPath_t1, "1A_peakAnnotation_summary.txt", sep="/") )
          print( peakAnno )  
      sink()
      write.table(x=as.data.frame(peakAnno), file = paste(myPath_t1, "1B_peakAnnotation_details.txt", sep="/"), append = FALSE,
                   quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
  
      pdf(file=paste(myPath_t1, "2A_GenomicDistribution_by_pieplot.pdf", sep="/"),  width=7, height=5)
          print( plotAnnoPie(peakAnno) ) 
      dev.off()
      svg(filename=paste(myPath_t1, "2B_GenomicDistribution_by_pieplot.svg", sep="/"),  width=7, height=5)
          print( plotAnnoPie(peakAnno) ) 
      dev.off()
  
      pdf(file=paste(myPath_t1, "3A_GenomicDistribution_by_barplot.pdf", sep="/"),  width=7, height=3)
          print( plotAnnoBar(peakAnno) ) 
      dev.off()  
      svg(filename=paste(myPath_t1, "3B_GenomicDistribution_by_barplot.svg", sep="/"),  width=7, height=3)
          print( plotAnnoBar(peakAnno) ) 
      dev.off()
  
      pdf(file=paste(myPath_t1, "4A_GenomicDistribution_by_vennpie.pdf", sep="/"),  width=7, height=5)
          print( vennpie(peakAnno) ) 
      dev.off()
      svg(filename=paste(myPath_t1, "4B_GenomicDistribution_by_vennpie.svg", sep="/"),  width=7, height=5)
          print( vennpie(peakAnno) ) 
      dev.off()
  
      pdf(file=paste(myPath_t1, "5A_GenomicDistribution_by_upsetplot.pdf", sep="/"),  width=8, height=5)
          print( upsetplot(peakAnno) ) 
      dev.off() 
      svg(filename=paste(myPath_t1, "5B_GenomicDistribution_by_upsetplot.svg", sep="/"),  width=8, height=5)
          print( upsetplot(peakAnno) ) 
      dev.off()
  
      pdf(file=paste(myPath_t1, "6A_GenomicDistribution_by_upsetplot_vennpie.pdf", sep="/"),  width=8, height=5)
          print( upsetplot(peakAnno, vennpie=TRUE) ) 
      dev.off()
      svg(filename=paste(myPath_t1, "6B_GenomicDistribution_by_upsetplot_vennpie.svg", sep="/"),  width=8, height=5)
          print( upsetplot(peakAnno, vennpie=TRUE) ) 
      dev.off()
  
      pdf(file=paste(myPath_t1, "7A_Distribution_of_peaks_to_TSS.pdf", sep="/"),  width=8, height=2)
          print( plotDistToTSS(peakAnno,  title="Distribution of peaks relative to TSS") ) 
      dev.off()
      svg(filename=paste(myPath_t1, "7B_Distribution_of_peaks_to_TSS.svg", sep="/"),  width=800, height=200)
          print( plotDistToTSS(peakAnno,  title="Distribution of peaks relative to TSS") ) 
      dev.off()
      
  }
  tryCatch(
    myTempFunction3(),
    error = function(err){"MyPeaksDistribution_OneGroup_1_g:error_3."}
  )
}




#############
MyGeneOntology_BP_OneGroup_1_g <- function(MyMaxDis=900000000,  MyPeaksAnno, MyFolder,  MyFileName ) {
  myTempFunction3 <- function() {    
    peakAnno2 <- as.data.frame(MyPeaksAnno)
    peakAnno2 <- peakAnno2[ abs(peakAnno2$distanceToTSS) < MyMaxDis, ]
    write.table(x=as.data.frame(peakAnno2), file = paste(MyFolder, "Genomic_regions_annotation_within-XXXkb.txt", sep="/"),       
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
    
    
    ## 只有这里需要改动
    enrich_raw <- clusterProfiler::enrichGO(gene=as.data.frame(peakAnno2)$geneId, OrgDb=my_orgdb_g,  keyType = "ENTREZID",  ont = "BP",
                                     pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,     readable = TRUE, pool = FALSE)
    
    
    #
    write.table(x= enrich_raw , file = paste(MyFolder, "/", MyFileName, "1_enrich_raw.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_2_enrich_raw.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_raw, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_raw, showCategory = 10) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw, showCategory = 10) )
      print( enrichplot::dotplot(enrich_raw, showCategory = 5) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_3_enrich_raw.svg", sep=""),  width=8, height=5)
        print( enrichplot::dotplot(enrich_raw, showCategory = 30) )
    dev.off()
    enrich_raw_dataFrame     <- as.data.frame( enrich_raw ) 
    if(nrow(enrich_raw_dataFrame) > 1) {
      enrich_raw_dataFrame$GeneRatio <- sapply(enrich_raw$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_raw_dataFrame$pvalue    <- -log10(enrich_raw_dataFrame$pvalue)
      enrich_raw_dataFrame$p.adjust  <- -log10(enrich_raw_dataFrame$p.adjust)
      enrich_raw_dataFrame$qvalue    <- -log10(enrich_raw_dataFrame$qvalue)
      enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[1:20, -8]
      if( nrow(enrich_raw_dataFrame) <20 ) {enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[, -8]}
      enrich_raw_dataFrame_selected <- transform(enrich_raw_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_raw_dataFrame_selected$GeneRatio) + min(enrich_raw_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_raw_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_raw_dataFrame_selected$ID, labels=rev(enrich_raw_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_3_enrich_raw", sep=""),  height1=6,  width1=5)
    }
    
    # simplify
    enrich_raw_simplify <- clusterProfiler::simplify(enrich_raw, cutoff=0.7, by="p.adjust", select_fun=min)
    write.table(x= enrich_raw_simplify , file = paste(MyFolder, "/", MyFileName, "4_enrich_raw_simplify.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_5_enrich_raw_simplify.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_raw_simplify, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_raw_simplify, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_raw_simplify, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw_simplify, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_raw_simplify, showCategory = 10) )
      print( enrichplot::emapplot(enrich_raw_simplify, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_raw_simplify, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw_simplify, showCategory = 10) )
      print( enrichplot::dotplot(enrich_raw_simplify, showCategory = 5) )
      print( enrichplot::emapplot(enrich_raw_simplify, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_raw_simplify, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw_simplify, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_6_enrich_raw_simplify.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_raw_simplify, showCategory = 30) )
    dev.off()
    enrich_raw_simplify_dataFrame     <- as.data.frame( enrich_raw_simplify ) 
    if(nrow(enrich_raw_simplify_dataFrame) > 1) {
      enrich_raw_simplify_dataFrame$GeneRatio <- sapply(enrich_raw_simplify$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_raw_simplify_dataFrame$pvalue    <- -log10(enrich_raw_simplify_dataFrame$pvalue)
      enrich_raw_simplify_dataFrame$p.adjust  <- -log10(enrich_raw_simplify_dataFrame$p.adjust)
      enrich_raw_simplify_dataFrame$qvalue    <- -log10(enrich_raw_simplify_dataFrame$qvalue)
      enrich_raw_simplify_dataFrame_selected <- enrich_raw_simplify_dataFrame[1:20, -8]
      if( nrow(enrich_raw_simplify_dataFrame) <20 ) {enrich_raw_simplify_dataFrame_selected <- enrich_raw_simplify_dataFrame[, -8]}
      enrich_raw_simplify_dataFrame_selected <- transform(enrich_raw_simplify_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_raw_simplify_dataFrame_selected$GeneRatio) + min(enrich_raw_simplify_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_raw_simplify_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_raw_simplify_dataFrame_selected$ID, labels=rev(enrich_raw_simplify_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_6_enrich_raw_simplify", sep=""),  height1=6,  width1=5)
    }
    
    
    myTempFunction98769 <- function() {
      dev.off()
      dev.off()
      dev.off()
      dev.off()
      dev.off()
    }
    tryCatch(
      myTempFunction98769(),
      error = function(err){"myTempFunction98769_abhfd"}
    )
    
    
    

    #########################
    enrich_drop1 <- dropGO(enrich_raw, level = c(1:4), term = NULL) ## drop specific GO terms or level
    
    write.table(x= enrich_drop1 , file = paste(MyFolder, "/", MyFileName, "7_enrich_drop1.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_8_enrich_drop1.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_drop1, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_drop1, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_drop1, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_drop1, showCategory = 10) )
      print( enrichplot::emapplot(enrich_drop1, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_drop1, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1, showCategory = 10) )
      print( enrichplot::dotplot(enrich_drop1, showCategory = 5) )
      print( enrichplot::emapplot(enrich_drop1, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_drop1, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_9_enrich_drop1.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_drop1, showCategory = 30) )
    dev.off()
    enrich_drop1_dataFrame     <- as.data.frame( enrich_drop1 ) 
    if(nrow(enrich_drop1_dataFrame) > 1) {
      enrich_drop1_dataFrame$GeneRatio <- sapply(enrich_drop1$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_drop1_dataFrame$pvalue    <- -log10(enrich_drop1_dataFrame$pvalue)
      enrich_drop1_dataFrame$p.adjust  <- -log10(enrich_drop1_dataFrame$p.adjust)
      enrich_drop1_dataFrame$qvalue    <- -log10(enrich_drop1_dataFrame$qvalue)
      enrich_drop1_dataFrame_selected <- enrich_drop1_dataFrame[1:20, -8]
      if( nrow(enrich_drop1_dataFrame) <20 ) {enrich_drop1_dataFrame_selected <- enrich_drop1_dataFrame[, -8]}
      enrich_drop1_dataFrame_selected <- transform(enrich_drop1_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_drop1_dataFrame_selected$GeneRatio) + min(enrich_drop1_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_drop1_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_drop1_dataFrame_selected$ID, labels=rev(enrich_drop1_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_9_enrich_drop1", sep=""),  height1=6,  width1=5)
    }
    
    # simplify
    enrich_drop1_simplify <- clusterProfiler::simplify(enrich_drop1, cutoff=0.7, by="p.adjust", select_fun=min)
    write.table(x= enrich_drop1_simplify , file = paste(MyFolder, "/", MyFileName, "10_enrich_drop1_simplify.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_11_enrich_drop1_simplify.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_drop1_simplify, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_drop1_simplify, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_drop1_simplify, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1_simplify, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_drop1_simplify, showCategory = 10) )
      print( enrichplot::emapplot(enrich_drop1_simplify, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_drop1_simplify, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1_simplify, showCategory = 10) )
      print( enrichplot::dotplot(enrich_drop1_simplify, showCategory = 5) )
      print( enrichplot::emapplot(enrich_drop1_simplify, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_drop1_simplify, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1_simplify, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_12_enrich_drop1_simplify.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_drop1_simplify, showCategory = 30) )
    dev.off()
    enrich_drop1_simplify_dataFrame     <- as.data.frame( enrich_drop1_simplify ) 
    if(nrow(enrich_drop1_simplify_dataFrame) > 1) {
      enrich_drop1_simplify_dataFrame$GeneRatio <- sapply(enrich_drop1_simplify$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_drop1_simplify_dataFrame$pvalue    <- -log10(enrich_drop1_simplify_dataFrame$pvalue)
      enrich_drop1_simplify_dataFrame$p.adjust  <- -log10(enrich_drop1_simplify_dataFrame$p.adjust)
      enrich_drop1_simplify_dataFrame$qvalue    <- -log10(enrich_drop1_simplify_dataFrame$qvalue)
      enrich_drop1_simplify_dataFrame_selected <- enrich_drop1_simplify_dataFrame[1:20, -8]
      if( nrow(enrich_drop1_simplify_dataFrame) <20 ) {enrich_drop1_simplify_dataFrame_selected <- enrich_drop1_simplify_dataFrame[, -8]}
      enrich_drop1_simplify_dataFrame_selected <- transform(enrich_drop1_simplify_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_drop1_simplify_dataFrame_selected$GeneRatio) + min(enrich_drop1_simplify_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_drop1_simplify_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_drop1_simplify_dataFrame_selected$ID, labels=rev(enrich_drop1_simplify_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_12_enrich_drop1_simplify", sep=""),  height1=6,  width1=5)
    }
    
    
    myTempFunction98769 <- function() {
      dev.off()
      dev.off()
      dev.off()
      dev.off()
      dev.off()
    }
    tryCatch(
      myTempFunction98769(),
      error = function(err){"myTempFunction98769_abhfd"}
    )
    
    
    
    
    
    #########################
    enrich_drop2 <- dropGO(enrich_raw, level = c(1:5), term = NULL) ## drop specific GO terms or level
    
    write.table(x= enrich_drop2 , file = paste(MyFolder, "/", MyFileName, "13_enrich_drop2.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_14_enrich_drop2.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_drop2, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_drop2, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_drop2, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_drop2, showCategory = 10) )
      print( enrichplot::emapplot(enrich_drop2, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_drop2, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2, showCategory = 10) )
      print( enrichplot::dotplot(enrich_drop2, showCategory = 5) )
      print( enrichplot::emapplot(enrich_drop2, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_drop2, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_15_enrich_drop2.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_drop2, showCategory = 30) )
    dev.off()
    enrich_drop2_dataFrame     <- as.data.frame( enrich_drop2 ) 
    if(nrow(enrich_drop2_dataFrame) > 1) {
      enrich_drop2_dataFrame$GeneRatio <- sapply(enrich_drop2$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_drop2_dataFrame$pvalue    <- -log10(enrich_drop2_dataFrame$pvalue)
      enrich_drop2_dataFrame$p.adjust  <- -log10(enrich_drop2_dataFrame$p.adjust)
      enrich_drop2_dataFrame$qvalue    <- -log10(enrich_drop2_dataFrame$qvalue)
      enrich_drop2_dataFrame_selected <- enrich_drop2_dataFrame[1:20, -8]
      if( nrow(enrich_drop2_dataFrame) <20 ) {enrich_drop2_dataFrame_selected <- enrich_drop2_dataFrame[, -8]}
      enrich_drop2_dataFrame_selected <- transform(enrich_drop2_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_drop2_dataFrame_selected$GeneRatio) + min(enrich_drop2_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_drop2_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_drop2_dataFrame_selected$ID, labels=rev(enrich_drop2_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_15_enrich_drop2", sep=""),  height1=6,  width1=5)
    }
    
    # simplify
    enrich_drop2_simplify <- clusterProfiler::simplify(enrich_drop2, cutoff=0.7, by="p.adjust", select_fun=min)
    write.table(x= enrich_drop2_simplify , file = paste(MyFolder, "/", MyFileName, "16_enrich_drop2_simplify.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_17_enrich_drop2_simplify.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_drop2_simplify, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_drop2_simplify, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_drop2_simplify, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2_simplify, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_drop2_simplify, showCategory = 10) )
      print( enrichplot::emapplot(enrich_drop2_simplify, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_drop2_simplify, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2_simplify, showCategory = 10) )
      print( enrichplot::dotplot(enrich_drop2_simplify, showCategory = 5) )
      print( enrichplot::emapplot(enrich_drop2_simplify, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_drop2_simplify, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2_simplify, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_18_enrich_drop2_simplify.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_drop2_simplify, showCategory = 30) )
    dev.off()
    enrich_drop2_simplify_dataFrame     <- as.data.frame( enrich_drop2_simplify ) 
    if(nrow(enrich_drop2_simplify_dataFrame) > 1) {
      enrich_drop2_simplify_dataFrame$GeneRatio <- sapply(enrich_drop2_simplify$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_drop2_simplify_dataFrame$pvalue    <- -log10(enrich_drop2_simplify_dataFrame$pvalue)
      enrich_drop2_simplify_dataFrame$p.adjust  <- -log10(enrich_drop2_simplify_dataFrame$p.adjust)
      enrich_drop2_simplify_dataFrame$qvalue    <- -log10(enrich_drop2_simplify_dataFrame$qvalue)
      enrich_drop2_simplify_dataFrame_selected <- enrich_drop2_simplify_dataFrame[1:20, -8]
      if( nrow(enrich_drop2_simplify_dataFrame) <20 ) {enrich_drop2_simplify_dataFrame_selected <- enrich_drop2_simplify_dataFrame[, -8]}
      enrich_drop2_simplify_dataFrame_selected <- transform(enrich_drop2_simplify_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_drop2_simplify_dataFrame_selected$GeneRatio) + min(enrich_drop2_simplify_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_drop2_simplify_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_drop2_simplify_dataFrame_selected$ID, labels=rev(enrich_drop2_simplify_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_18_enrich_drop2_simplify", sep=""),  height1=6,  width1=5)
    }
    
    
    myTempFunction98769 <- function() {
      dev.off()
      dev.off()
      dev.off()
      dev.off()
      dev.off()
    }
    tryCatch(
      myTempFunction98769(),
      error = function(err){"myTempFunction98769_abhfd"}
    )
  
      
  }
  tryCatch(
    myTempFunction3(),
    error = function(err){"MyPeaksDistribution_OneGroup_1_g:error_3."}
  )  
}




#############
MyGeneOntology_MF_OneGroup_1_g <- function(MyMaxDis=900000000,  MyPeaksAnno, MyFolder,  MyFileName ) {
  myTempFunction3 <- function() {    
    peakAnno2 <- as.data.frame(MyPeaksAnno)
    peakAnno2 <- peakAnno2[ abs(peakAnno2$distanceToTSS) < MyMaxDis, ]
    write.table(x=as.data.frame(peakAnno2), file = paste(MyFolder, "Genomic_regions_annotation_within-XXXkb.txt", sep="/"),       
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
    
    
    ## 只有这里需要改动
    enrich_raw <- clusterProfiler::enrichGO(gene=as.data.frame(peakAnno2)$geneId, OrgDb=my_orgdb_g,  keyType = "ENTREZID",  ont = "MF",
                                     pvalueCutoff = 0.05, pAdjustMethod = "BH",   qvalueCutoff = 0.2,  readable = TRUE, pool = FALSE)
    
    
    #
    write.table(x= enrich_raw , file = paste(MyFolder, "/", MyFileName, "1_enrich_raw.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_2_enrich_raw.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_raw, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_raw, showCategory = 10) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw, showCategory = 10) )
      print( enrichplot::dotplot(enrich_raw, showCategory = 5) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_3_enrich_raw.svg", sep=""),  width=8, height=5)
        print( enrichplot::dotplot(enrich_raw, showCategory = 30) )
    dev.off()
    enrich_raw_dataFrame     <- as.data.frame( enrich_raw ) 
    if(nrow(enrich_raw_dataFrame) > 1) {
      enrich_raw_dataFrame$GeneRatio <- sapply(enrich_raw$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_raw_dataFrame$pvalue    <- -log10(enrich_raw_dataFrame$pvalue)
      enrich_raw_dataFrame$p.adjust  <- -log10(enrich_raw_dataFrame$p.adjust)
      enrich_raw_dataFrame$qvalue    <- -log10(enrich_raw_dataFrame$qvalue)
      enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[1:20, -8]
      if( nrow(enrich_raw_dataFrame) <20 ) {enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[, -8]}
      enrich_raw_dataFrame_selected <- transform(enrich_raw_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_raw_dataFrame_selected$GeneRatio) + min(enrich_raw_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_raw_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_raw_dataFrame_selected$ID, labels=rev(enrich_raw_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_3_enrich_raw", sep=""),  height1=6,  width1=5)
    }
    
    # simplify
    enrich_raw_simplify <- clusterProfiler::simplify(enrich_raw, cutoff=0.7, by="p.adjust", select_fun=min)
    write.table(x= enrich_raw_simplify , file = paste(MyFolder, "/", MyFileName, "4_enrich_raw_simplify.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_5_enrich_raw_simplify.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_raw_simplify, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_raw_simplify, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_raw_simplify, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw_simplify, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_raw_simplify, showCategory = 10) )
      print( enrichplot::emapplot(enrich_raw_simplify, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_raw_simplify, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw_simplify, showCategory = 10) )
      print( enrichplot::dotplot(enrich_raw_simplify, showCategory = 5) )
      print( enrichplot::emapplot(enrich_raw_simplify, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_raw_simplify, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw_simplify, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_6_enrich_raw_simplify.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_raw_simplify, showCategory = 30) )
    dev.off()
    enrich_raw_simplify_dataFrame     <- as.data.frame( enrich_raw_simplify ) 
    if(nrow(enrich_raw_simplify_dataFrame) > 1) {
      enrich_raw_simplify_dataFrame$GeneRatio <- sapply(enrich_raw_simplify$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_raw_simplify_dataFrame$pvalue    <- -log10(enrich_raw_simplify_dataFrame$pvalue)
      enrich_raw_simplify_dataFrame$p.adjust  <- -log10(enrich_raw_simplify_dataFrame$p.adjust)
      enrich_raw_simplify_dataFrame$qvalue    <- -log10(enrich_raw_simplify_dataFrame$qvalue)
      enrich_raw_simplify_dataFrame_selected <- enrich_raw_simplify_dataFrame[1:20, -8]
      if( nrow(enrich_raw_simplify_dataFrame) <20 ) {enrich_raw_simplify_dataFrame_selected <- enrich_raw_simplify_dataFrame[, -8]}
      enrich_raw_simplify_dataFrame_selected <- transform(enrich_raw_simplify_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_raw_simplify_dataFrame_selected$GeneRatio) + min(enrich_raw_simplify_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_raw_simplify_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_raw_simplify_dataFrame_selected$ID, labels=rev(enrich_raw_simplify_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_6_enrich_raw_simplify", sep=""),  height1=6,  width1=5)
    }
    
    
    myTempFunction98769 <- function() {
      dev.off()
      dev.off()
      dev.off()
      dev.off()
      dev.off()
    }
    tryCatch(
      myTempFunction98769(),
      error = function(err){"myTempFunction98769_abhfd"}
    )
    
    
    

    #########################
    enrich_drop1 <- dropGO(enrich_raw, level = c(1:4), term = NULL) ## drop specific GO terms or level
    
    write.table(x= enrich_drop1 , file = paste(MyFolder, "/", MyFileName, "7_enrich_drop1.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_8_enrich_drop1.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_drop1, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_drop1, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_drop1, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_drop1, showCategory = 10) )
      print( enrichplot::emapplot(enrich_drop1, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_drop1, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1, showCategory = 10) )
      print( enrichplot::dotplot(enrich_drop1, showCategory = 5) )
      print( enrichplot::emapplot(enrich_drop1, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_drop1, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_9_enrich_drop1.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_drop1, showCategory = 30) )
    dev.off()
    enrich_drop1_dataFrame     <- as.data.frame( enrich_drop1 ) 
    if(nrow(enrich_drop1_dataFrame) > 1) {
      enrich_drop1_dataFrame$GeneRatio <- sapply(enrich_drop1$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_drop1_dataFrame$pvalue    <- -log10(enrich_drop1_dataFrame$pvalue)
      enrich_drop1_dataFrame$p.adjust  <- -log10(enrich_drop1_dataFrame$p.adjust)
      enrich_drop1_dataFrame$qvalue    <- -log10(enrich_drop1_dataFrame$qvalue)
      enrich_drop1_dataFrame_selected <- enrich_drop1_dataFrame[1:20, -8]
      if( nrow(enrich_drop1_dataFrame) <20 ) {enrich_drop1_dataFrame_selected <- enrich_drop1_dataFrame[, -8]}
      enrich_drop1_dataFrame_selected <- transform(enrich_drop1_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_drop1_dataFrame_selected$GeneRatio) + min(enrich_drop1_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_drop1_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_drop1_dataFrame_selected$ID, labels=rev(enrich_drop1_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_9_enrich_drop1", sep=""),  height1=6,  width1=5)
    }
    
    # simplify
    enrich_drop1_simplify <- clusterProfiler::simplify(enrich_drop1, cutoff=0.7, by="p.adjust", select_fun=min)
    write.table(x= enrich_drop1_simplify , file = paste(MyFolder, "/", MyFileName, "10_enrich_drop1_simplify.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_11_enrich_drop1_simplify.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_drop1_simplify, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_drop1_simplify, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_drop1_simplify, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1_simplify, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_drop1_simplify, showCategory = 10) )
      print( enrichplot::emapplot(enrich_drop1_simplify, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_drop1_simplify, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1_simplify, showCategory = 10) )
      print( enrichplot::dotplot(enrich_drop1_simplify, showCategory = 5) )
      print( enrichplot::emapplot(enrich_drop1_simplify, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_drop1_simplify, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1_simplify, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_12_enrich_drop1_simplify.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_drop1_simplify, showCategory = 30) )
    dev.off()
    enrich_drop1_simplify_dataFrame     <- as.data.frame( enrich_drop1_simplify ) 
    if(nrow(enrich_drop1_simplify_dataFrame) > 1) {
      enrich_drop1_simplify_dataFrame$GeneRatio <- sapply(enrich_drop1_simplify$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_drop1_simplify_dataFrame$pvalue    <- -log10(enrich_drop1_simplify_dataFrame$pvalue)
      enrich_drop1_simplify_dataFrame$p.adjust  <- -log10(enrich_drop1_simplify_dataFrame$p.adjust)
      enrich_drop1_simplify_dataFrame$qvalue    <- -log10(enrich_drop1_simplify_dataFrame$qvalue)
      enrich_drop1_simplify_dataFrame_selected <- enrich_drop1_simplify_dataFrame[1:20, -8]
      if( nrow(enrich_drop1_simplify_dataFrame) <20 ) {enrich_drop1_simplify_dataFrame_selected <- enrich_drop1_simplify_dataFrame[, -8]}
      enrich_drop1_simplify_dataFrame_selected <- transform(enrich_drop1_simplify_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_drop1_simplify_dataFrame_selected$GeneRatio) + min(enrich_drop1_simplify_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_drop1_simplify_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_drop1_simplify_dataFrame_selected$ID, labels=rev(enrich_drop1_simplify_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_12_enrich_drop1_simplify", sep=""),  height1=6,  width1=5)
    }
    
    
    myTempFunction98769 <- function() {
      dev.off()
      dev.off()
      dev.off()
      dev.off()
      dev.off()
    }
    tryCatch(
      myTempFunction98769(),
      error = function(err){"myTempFunction98769_abhfd"}
    )
    
    
    
    
    
    #########################
    enrich_drop2 <- dropGO(enrich_raw, level = c(1:5), term = NULL) ## drop specific GO terms or level
    
    write.table(x= enrich_drop2 , file = paste(MyFolder, "/", MyFileName, "13_enrich_drop2.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_14_enrich_drop2.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_drop2, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_drop2, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_drop2, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_drop2, showCategory = 10) )
      print( enrichplot::emapplot(enrich_drop2, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_drop2, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2, showCategory = 10) )
      print( enrichplot::dotplot(enrich_drop2, showCategory = 5) )
      print( enrichplot::emapplot(enrich_drop2, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_drop2, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_15_enrich_drop2.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_drop2, showCategory = 30) )
    dev.off()
    enrich_drop2_dataFrame     <- as.data.frame( enrich_drop2 ) 
    if(nrow(enrich_drop2_dataFrame) > 1) {
      enrich_drop2_dataFrame$GeneRatio <- sapply(enrich_drop2$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_drop2_dataFrame$pvalue    <- -log10(enrich_drop2_dataFrame$pvalue)
      enrich_drop2_dataFrame$p.adjust  <- -log10(enrich_drop2_dataFrame$p.adjust)
      enrich_drop2_dataFrame$qvalue    <- -log10(enrich_drop2_dataFrame$qvalue)
      enrich_drop2_dataFrame_selected <- enrich_drop2_dataFrame[1:20, -8]
      if( nrow(enrich_drop2_dataFrame) <20 ) {enrich_drop2_dataFrame_selected <- enrich_drop2_dataFrame[, -8]}
      enrich_drop2_dataFrame_selected <- transform(enrich_drop2_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_drop2_dataFrame_selected$GeneRatio) + min(enrich_drop2_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_drop2_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_drop2_dataFrame_selected$ID, labels=rev(enrich_drop2_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_15_enrich_drop2", sep=""),  height1=6,  width1=5)
    }
    
    # simplify
    enrich_drop2_simplify <- clusterProfiler::simplify(enrich_drop2, cutoff=0.7, by="p.adjust", select_fun=min)
    write.table(x= enrich_drop2_simplify , file = paste(MyFolder, "/", MyFileName, "16_enrich_drop2_simplify.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_17_enrich_drop2_simplify.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_drop2_simplify, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_drop2_simplify, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_drop2_simplify, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2_simplify, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_drop2_simplify, showCategory = 10) )
      print( enrichplot::emapplot(enrich_drop2_simplify, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_drop2_simplify, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2_simplify, showCategory = 10) )
      print( enrichplot::dotplot(enrich_drop2_simplify, showCategory = 5) )
      print( enrichplot::emapplot(enrich_drop2_simplify, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_drop2_simplify, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2_simplify, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_18_enrich_drop2_simplify.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_drop2_simplify, showCategory = 30) )
    dev.off()
    enrich_drop2_simplify_dataFrame     <- as.data.frame( enrich_drop2_simplify ) 
    if(nrow(enrich_drop2_simplify_dataFrame) > 1) {
      enrich_drop2_simplify_dataFrame$GeneRatio <- sapply(enrich_drop2_simplify$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_drop2_simplify_dataFrame$pvalue    <- -log10(enrich_drop2_simplify_dataFrame$pvalue)
      enrich_drop2_simplify_dataFrame$p.adjust  <- -log10(enrich_drop2_simplify_dataFrame$p.adjust)
      enrich_drop2_simplify_dataFrame$qvalue    <- -log10(enrich_drop2_simplify_dataFrame$qvalue)
      enrich_drop2_simplify_dataFrame_selected <- enrich_drop2_simplify_dataFrame[1:20, -8]
      if( nrow(enrich_drop2_simplify_dataFrame) <20 ) {enrich_drop2_simplify_dataFrame_selected <- enrich_drop2_simplify_dataFrame[, -8]}
      enrich_drop2_simplify_dataFrame_selected <- transform(enrich_drop2_simplify_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_drop2_simplify_dataFrame_selected$GeneRatio) + min(enrich_drop2_simplify_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_drop2_simplify_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_drop2_simplify_dataFrame_selected$ID, labels=rev(enrich_drop2_simplify_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_18_enrich_drop2_simplify", sep=""),  height1=6,  width1=5)
    }
    
    
    myTempFunction98769 <- function() {
      dev.off()
      dev.off()
      dev.off()
      dev.off()
      dev.off()
    }
    tryCatch(
      myTempFunction98769(),
      error = function(err){"myTempFunction98769_abhfd"}
    )
    
      
  }
  tryCatch(
    myTempFunction3(),
    error = function(err){"MyPeaksDistribution_OneGroup_1_g:error_3."}
  )
}




#############
MyGeneOntology_CC_OneGroup_1_g <- function(MyMaxDis=900000000,  MyPeaksAnno, MyFolder,  MyFileName ) {
  myTempFunction3 <- function() {    
    peakAnno2 <- as.data.frame(MyPeaksAnno)
    peakAnno2 <- peakAnno2[ abs(peakAnno2$distanceToTSS) < MyMaxDis, ]
    write.table(x=as.data.frame(peakAnno2), file = paste(MyFolder, "Genomic_regions_annotation_within-XXXkb.txt", sep="/"),       
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
    
    
    ## 只有这里需要改动
    enrich_raw <- clusterProfiler::enrichGO(gene=as.data.frame(peakAnno2)$geneId, OrgDb=my_orgdb_g,  keyType = "ENTREZID",  ont = "CC",
                                     pvalueCutoff = 0.05, pAdjustMethod = "BH",   qvalueCutoff = 0.2,   readable = TRUE, pool = FALSE)
    
    
    #
    write.table(x= enrich_raw , file = paste(MyFolder, "/", MyFileName, "1_enrich_raw.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_2_enrich_raw.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_raw, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_raw, showCategory = 10) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw, showCategory = 10) )
      print( enrichplot::dotplot(enrich_raw, showCategory = 5) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_3_enrich_raw.svg", sep=""),  width=8, height=5)
        print( enrichplot::dotplot(enrich_raw, showCategory = 30) )
    dev.off()
    enrich_raw_dataFrame     <- as.data.frame( enrich_raw ) 
    if(nrow(enrich_raw_dataFrame) > 1) {
      enrich_raw_dataFrame$GeneRatio <- sapply(enrich_raw$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_raw_dataFrame$pvalue    <- -log10(enrich_raw_dataFrame$pvalue)
      enrich_raw_dataFrame$p.adjust  <- -log10(enrich_raw_dataFrame$p.adjust)
      enrich_raw_dataFrame$qvalue    <- -log10(enrich_raw_dataFrame$qvalue)
      enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[1:20, -8]
      if( nrow(enrich_raw_dataFrame) <20 ) {enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[, -8]}
      enrich_raw_dataFrame_selected <- transform(enrich_raw_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_raw_dataFrame_selected$GeneRatio) + min(enrich_raw_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_raw_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_raw_dataFrame_selected$ID, labels=rev(enrich_raw_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_3_enrich_raw", sep=""),  height1=6,  width1=5)
    }
    
    # simplify
    enrich_raw_simplify <- clusterProfiler::simplify(enrich_raw, cutoff=0.7, by="p.adjust", select_fun=min)
    write.table(x= enrich_raw_simplify , file = paste(MyFolder, "/", MyFileName, "4_enrich_raw_simplify.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_5_enrich_raw_simplify.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_raw_simplify, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_raw_simplify, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_raw_simplify, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw_simplify, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_raw_simplify, showCategory = 10) )
      print( enrichplot::emapplot(enrich_raw_simplify, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_raw_simplify, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw_simplify, showCategory = 10) )
      print( enrichplot::dotplot(enrich_raw_simplify, showCategory = 5) )
      print( enrichplot::emapplot(enrich_raw_simplify, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_raw_simplify, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_raw_simplify, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_6_enrich_raw_simplify.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_raw_simplify, showCategory = 30) )
    dev.off()
    enrich_raw_simplify_dataFrame     <- as.data.frame( enrich_raw_simplify ) 
    if(nrow(enrich_raw_simplify_dataFrame) > 1) {
      enrich_raw_simplify_dataFrame$GeneRatio <- sapply(enrich_raw_simplify$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_raw_simplify_dataFrame$pvalue    <- -log10(enrich_raw_simplify_dataFrame$pvalue)
      enrich_raw_simplify_dataFrame$p.adjust  <- -log10(enrich_raw_simplify_dataFrame$p.adjust)
      enrich_raw_simplify_dataFrame$qvalue    <- -log10(enrich_raw_simplify_dataFrame$qvalue)
      enrich_raw_simplify_dataFrame_selected <- enrich_raw_simplify_dataFrame[1:20, -8]
      if( nrow(enrich_raw_simplify_dataFrame) <20 ) {enrich_raw_simplify_dataFrame_selected <- enrich_raw_simplify_dataFrame[, -8]}
      enrich_raw_simplify_dataFrame_selected <- transform(enrich_raw_simplify_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_raw_simplify_dataFrame_selected$GeneRatio) + min(enrich_raw_simplify_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_raw_simplify_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_raw_simplify_dataFrame_selected$ID, labels=rev(enrich_raw_simplify_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_6_enrich_raw_simplify", sep=""),  height1=6,  width1=5)
    }
    
    
    myTempFunction98769 <- function() {
      dev.off()
      dev.off()
      dev.off()
      dev.off()
      dev.off()
    }
    tryCatch(
      myTempFunction98769(),
      error = function(err){"myTempFunction98769_abhfd"}
    )
    
    
    

    #########################
    enrich_drop1 <- dropGO(enrich_raw, level = c(1:4), term = NULL) ## drop specific GO terms or level
    
    write.table(x= enrich_drop1 , file = paste(MyFolder, "/", MyFileName, "7_enrich_drop1.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_8_enrich_drop1.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_drop1, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_drop1, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_drop1, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_drop1, showCategory = 10) )
      print( enrichplot::emapplot(enrich_drop1, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_drop1, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1, showCategory = 10) )
      print( enrichplot::dotplot(enrich_drop1, showCategory = 5) )
      print( enrichplot::emapplot(enrich_drop1, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_drop1, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_9_enrich_drop1.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_drop1, showCategory = 30) )
    dev.off()
    enrich_drop1_dataFrame     <- as.data.frame( enrich_drop1 ) 
    if(nrow(enrich_drop1_dataFrame) > 1) {
      enrich_drop1_dataFrame$GeneRatio <- sapply(enrich_drop1$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_drop1_dataFrame$pvalue    <- -log10(enrich_drop1_dataFrame$pvalue)
      enrich_drop1_dataFrame$p.adjust  <- -log10(enrich_drop1_dataFrame$p.adjust)
      enrich_drop1_dataFrame$qvalue    <- -log10(enrich_drop1_dataFrame$qvalue)
      enrich_drop1_dataFrame_selected <- enrich_drop1_dataFrame[1:20, -8]
      if( nrow(enrich_drop1_dataFrame) <20 ) {enrich_drop1_dataFrame_selected <- enrich_drop1_dataFrame[, -8]}
      enrich_drop1_dataFrame_selected <- transform(enrich_drop1_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_drop1_dataFrame_selected$GeneRatio) + min(enrich_drop1_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_drop1_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_drop1_dataFrame_selected$ID, labels=rev(enrich_drop1_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_9_enrich_drop1", sep=""),  height1=6,  width1=5)
    }
    
    # simplify
    enrich_drop1_simplify <- clusterProfiler::simplify(enrich_drop1, cutoff=0.7, by="p.adjust", select_fun=min)
    write.table(x= enrich_drop1_simplify , file = paste(MyFolder, "/", MyFileName, "10_enrich_drop1_simplify.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_11_enrich_drop1_simplify.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_drop1_simplify, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_drop1_simplify, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_drop1_simplify, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1_simplify, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_drop1_simplify, showCategory = 10) )
      print( enrichplot::emapplot(enrich_drop1_simplify, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_drop1_simplify, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1_simplify, showCategory = 10) )
      print( enrichplot::dotplot(enrich_drop1_simplify, showCategory = 5) )
      print( enrichplot::emapplot(enrich_drop1_simplify, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_drop1_simplify, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop1_simplify, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_12_enrich_drop1_simplify.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_drop1_simplify, showCategory = 30) )
    dev.off()
    enrich_drop1_simplify_dataFrame     <- as.data.frame( enrich_drop1_simplify ) 
    if(nrow(enrich_drop1_simplify_dataFrame) > 1) {
      enrich_drop1_simplify_dataFrame$GeneRatio <- sapply(enrich_drop1_simplify$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_drop1_simplify_dataFrame$pvalue    <- -log10(enrich_drop1_simplify_dataFrame$pvalue)
      enrich_drop1_simplify_dataFrame$p.adjust  <- -log10(enrich_drop1_simplify_dataFrame$p.adjust)
      enrich_drop1_simplify_dataFrame$qvalue    <- -log10(enrich_drop1_simplify_dataFrame$qvalue)
      enrich_drop1_simplify_dataFrame_selected <- enrich_drop1_simplify_dataFrame[1:20, -8]
      if( nrow(enrich_drop1_simplify_dataFrame) <20 ) {enrich_drop1_simplify_dataFrame_selected <- enrich_drop1_simplify_dataFrame[, -8]}
      enrich_drop1_simplify_dataFrame_selected <- transform(enrich_drop1_simplify_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_drop1_simplify_dataFrame_selected$GeneRatio) + min(enrich_drop1_simplify_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_drop1_simplify_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_drop1_simplify_dataFrame_selected$ID, labels=rev(enrich_drop1_simplify_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_12_enrich_drop1_simplify", sep=""),  height1=6,  width1=5)
    }
    
    
    myTempFunction98769 <- function() {
      dev.off()
      dev.off()
      dev.off()
      dev.off()
      dev.off()
    }
    tryCatch(
      myTempFunction98769(),
      error = function(err){"myTempFunction98769_abhfd"}
    )
    
    
    
    
    
    #########################
    enrich_drop2 <- dropGO(enrich_raw, level = c(1:5), term = NULL) ## drop specific GO terms or level
    
    write.table(x= enrich_drop2 , file = paste(MyFolder, "/", MyFileName, "13_enrich_drop2.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_14_enrich_drop2.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_drop2, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_drop2, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_drop2, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_drop2, showCategory = 10) )
      print( enrichplot::emapplot(enrich_drop2, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_drop2, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2, showCategory = 10) )
      print( enrichplot::dotplot(enrich_drop2, showCategory = 5) )
      print( enrichplot::emapplot(enrich_drop2, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_drop2, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_15_enrich_drop2.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_drop2, showCategory = 30) )
    dev.off()
    enrich_drop2_dataFrame     <- as.data.frame( enrich_drop2 ) 
    if(nrow(enrich_drop2_dataFrame) > 1) {
      enrich_drop2_dataFrame$GeneRatio <- sapply(enrich_drop2$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_drop2_dataFrame$pvalue    <- -log10(enrich_drop2_dataFrame$pvalue)
      enrich_drop2_dataFrame$p.adjust  <- -log10(enrich_drop2_dataFrame$p.adjust)
      enrich_drop2_dataFrame$qvalue    <- -log10(enrich_drop2_dataFrame$qvalue)
      enrich_drop2_dataFrame_selected <- enrich_drop2_dataFrame[1:20, -8]
      if( nrow(enrich_drop2_dataFrame) <20 ) {enrich_drop2_dataFrame_selected <- enrich_drop2_dataFrame[, -8]}
      enrich_drop2_dataFrame_selected <- transform(enrich_drop2_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_drop2_dataFrame_selected$GeneRatio) + min(enrich_drop2_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_drop2_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_drop2_dataFrame_selected$ID, labels=rev(enrich_drop2_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_15_enrich_drop2", sep=""),  height1=6,  width1=5)
    }
    
    # simplify
    enrich_drop2_simplify <- clusterProfiler::simplify(enrich_drop2, cutoff=0.7, by="p.adjust", select_fun=min)
    write.table(x= enrich_drop2_simplify , file = paste(MyFolder, "/", MyFileName, "16_enrich_drop2_simplify.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_17_enrich_drop2_simplify.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_drop2_simplify, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_drop2_simplify, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_drop2_simplify, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2_simplify, showCategory = 20 ) )
      print( enrichplot::dotplot(enrich_drop2_simplify, showCategory = 10) )
      print( enrichplot::emapplot(enrich_drop2_simplify, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_drop2_simplify, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2_simplify, showCategory = 10) )
      print( enrichplot::dotplot(enrich_drop2_simplify, showCategory = 5) )
      print( enrichplot::emapplot(enrich_drop2_simplify, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_drop2_simplify, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::goplot(enrich_drop2_simplify, showCategory = 5) )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_18_enrich_drop2_simplify.svg", sep=""),  width=8, height=5)
    print( enrichplot::dotplot(enrich_drop2_simplify, showCategory = 30) )
    dev.off()
    enrich_drop2_simplify_dataFrame     <- as.data.frame( enrich_drop2_simplify ) 
    if(nrow(enrich_drop2_simplify_dataFrame) > 1) {
      enrich_drop2_simplify_dataFrame$GeneRatio <- sapply(enrich_drop2_simplify$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_drop2_simplify_dataFrame$pvalue    <- -log10(enrich_drop2_simplify_dataFrame$pvalue)
      enrich_drop2_simplify_dataFrame$p.adjust  <- -log10(enrich_drop2_simplify_dataFrame$p.adjust)
      enrich_drop2_simplify_dataFrame$qvalue    <- -log10(enrich_drop2_simplify_dataFrame$qvalue)
      enrich_drop2_simplify_dataFrame_selected <- enrich_drop2_simplify_dataFrame[1:20, -8]
      if( nrow(enrich_drop2_simplify_dataFrame) <20 ) {enrich_drop2_simplify_dataFrame_selected <- enrich_drop2_simplify_dataFrame[, -8]}
      enrich_drop2_simplify_dataFrame_selected <- transform(enrich_drop2_simplify_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_drop2_simplify_dataFrame_selected$GeneRatio) + min(enrich_drop2_simplify_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_drop2_simplify_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_drop2_simplify_dataFrame_selected$ID, labels=rev(enrich_drop2_simplify_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_18_enrich_drop2_simplify", sep=""),  height1=6,  width1=5)
    }
    
    
    myTempFunction98769 <- function() {
      dev.off()
      dev.off()
      dev.off()
      dev.off()
      dev.off()
    }
    tryCatch(
      myTempFunction98769(),
      error = function(err){"myTempFunction98769_abhfd"}
    )
   
      
  }
  tryCatch(
    myTempFunction3(),
    error = function(err){"MyPeaksDistribution_OneGroup_1_g:error_3."}
  ) 
}





#############
MyReactome_OneGroup_1_g <- function(MyMaxDis=900000000,  MyPeaksAnno, MyFolder,  MyFileName ) {
  myTempFunction3 <- function() {    
    peakAnno2 <- as.data.frame(MyPeaksAnno)
    peakAnno2 <- peakAnno2[ abs(peakAnno2$distanceToTSS) < MyMaxDis, ]
    write.table(x=as.data.frame(peakAnno2), file = paste(MyFolder, "Genomic_regions_annotation_within-XXXkb.txt", sep="/"),       
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
    
    
    ## 只有这里需要改动
    enrich_raw <- enrichPathway(gene=as.data.frame(peakAnno2)$geneId, organism = my_organism_g, pvalueCutoff = 0.05,
                       pAdjustMethod = "BH", qvalueCutoff = 0.2,     readable = TRUE)
    
    #
    write.table(x= enrich_raw , file = paste(MyFolder, "/", MyFileName, "1_enrich_raw.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_2_enrich_raw.pdf", sep=""),  width=8, height=5)
      print( enrichplot::dotplot(enrich_raw, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::dotplot(enrich_raw, showCategory = 10) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::dotplot(enrich_raw, showCategory = 5) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_3_enrich_raw.svg", sep=""),  width=8, height=5)
        print( enrichplot::dotplot(enrich_raw, showCategory = 30) )
    dev.off()
    enrich_raw_dataFrame     <- as.data.frame( enrich_raw ) 
    if(nrow(enrich_raw_dataFrame) > 1) {
      enrich_raw_dataFrame$GeneRatio <- sapply(enrich_raw$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_raw_dataFrame$pvalue    <- -log10(enrich_raw_dataFrame$pvalue)
      enrich_raw_dataFrame$p.adjust  <- -log10(enrich_raw_dataFrame$p.adjust)
      enrich_raw_dataFrame$qvalue    <- -log10(enrich_raw_dataFrame$qvalue)
      enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[1:20, -8]
      if( nrow(enrich_raw_dataFrame) <20 ) {enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[, -8]}
      enrich_raw_dataFrame_selected <- transform(enrich_raw_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_raw_dataFrame_selected$GeneRatio) + min(enrich_raw_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_raw_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_raw_dataFrame_selected$ID, labels=rev(enrich_raw_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_3_enrich_raw", sep=""),  height1=6,  width1=5)
    }
 
    myTempFunction98769 <- function() {
      dev.off()
      dev.off()
      dev.off()
      dev.off()
      dev.off()
    }
    tryCatch(
      myTempFunction98769(),
      error = function(err){"myTempFunction98769_abhfd"}
    )
     
      
  }
  tryCatch(
    myTempFunction3(),
    error = function(err){"MyPeaksDistribution_OneGroup_1_g:error_3."}
  )
}




#############
MyKEGG_OneGroup_1_g <- function(MyMaxDis=900000000,  MyPeaksAnno, MyFolder,  MyFileName ) {  
  myTempFunction3 <- function() {    
    peakAnno2 <- as.data.frame(MyPeaksAnno)
    peakAnno2 <- peakAnno2[ abs(peakAnno2$distanceToTSS) < MyMaxDis, ]
    write.table(x=as.data.frame(peakAnno2), file = paste(MyFolder, "Genomic_regions_annotation_within-XXXkb.txt", sep="/"),       
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
    
    
    ## 只有这里需要改动
    enrich_raw <- enrichKEGG(gene=as.data.frame(peakAnno2)$geneId, organism = my_organism2_g, keyType = "kegg",  pvalueCutoff = 0.05,
                       pAdjustMethod = "BH", qvalueCutoff = 0.2,      use_internal_data = TRUE )

  
    #
    write.table(x= enrich_raw , file = paste(MyFolder, "/", MyFileName, "1_enrich_raw.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_2_enrich_raw.pdf", sep=""),  width=8, height=5)
      print( barplot(enrich_raw, drop=TRUE, showCategory=20, x = "Count",     color = "p.adjust", font.size = 12, title = "")  )
      print( barplot(enrich_raw, drop=TRUE, showCategory=20, x = "GeneRatio", color = "p.adjust", font.size = 12, title = "")  )
      print( enrichplot::dotplot(enrich_raw, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::dotplot(enrich_raw, showCategory = 10) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::dotplot(enrich_raw, showCategory = 5) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_3_enrich_raw.svg", sep=""),  width=8, height=5)
        print( enrichplot::dotplot(enrich_raw, showCategory = 30) )
    dev.off()
    enrich_raw_dataFrame     <- as.data.frame( enrich_raw ) 
    if(nrow(enrich_raw_dataFrame) > 1) {
      enrich_raw_dataFrame$GeneRatio <- sapply(enrich_raw$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_raw_dataFrame$pvalue    <- -log10(enrich_raw_dataFrame$pvalue)
      enrich_raw_dataFrame$p.adjust  <- -log10(enrich_raw_dataFrame$p.adjust)
      enrich_raw_dataFrame$qvalue    <- -log10(enrich_raw_dataFrame$qvalue)
      enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[1:20, -8]
      if( nrow(enrich_raw_dataFrame) <20 ) {enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[, -8]}
      enrich_raw_dataFrame_selected <- transform(enrich_raw_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_raw_dataFrame_selected$GeneRatio) + min(enrich_raw_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_raw_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_raw_dataFrame_selected$ID, labels=rev(enrich_raw_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_3_enrich_raw", sep=""),  height1=6,  width1=5)
    }
  
    myTempFunction98769 <- function() {
      dev.off()
      dev.off()
      dev.off()
      dev.off()
      dev.off()
    }
    tryCatch(
      myTempFunction98769(),
      error = function(err){"myTempFunction98769_abhfd"}
    )
 
      
  }
  tryCatch(
    myTempFunction3(),
    error = function(err){"MyPeaksDistribution_OneGroup_1_g:error_3."}
  )
}





#############
MyModuleKEGG_OneGroup_1_g <- function(MyMaxDis=900000000,  MyPeaksAnno, MyFolder,  MyFileName ) {
  myTempFunction3 <- function() {    
    peakAnno2 <- as.data.frame(MyPeaksAnno)
    peakAnno2 <- peakAnno2[ abs(peakAnno2$distanceToTSS) < MyMaxDis, ]
    write.table(x=as.data.frame(peakAnno2), file = paste(MyFolder, "Genomic_regions_annotation_within-XXXkb.txt", sep="/"),       
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
    
    
    ## 只有这里需要改动
    enrich_raw <- enrichMKEGG(gene=as.data.frame(peakAnno2)$geneId, organism = my_organism2_g, keyType = "kegg",  pvalueCutoff = 0.05,
                       pAdjustMethod = "BH", qvalueCutoff = 0.2     )


    #
    write.table(x= enrich_raw , file = paste(MyFolder, "/", MyFileName, "1_enrich_raw.txt", sep=""), 
                append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
    myTempFunction1234gty <- function() {
      pdf(file=paste(MyFolder, "/", MyFileName, "_2_enrich_raw.pdf", sep=""),  width=8, height=5)
      print( barplot(enrich_raw, drop=TRUE, showCategory=20, x = "Count",     color = "p.adjust", font.size = 12, title = "")  )
      print( barplot(enrich_raw, drop=TRUE, showCategory=20, x = "GeneRatio", color = "p.adjust", font.size = 12, title = "")  )
      print( enrichplot::dotplot(enrich_raw, showCategory = 20 ) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 20) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::dotplot(enrich_raw, showCategory = 10) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 10) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
      print( enrichplot::dotplot(enrich_raw, showCategory = 5) )
      print( enrichplot::emapplot(enrich_raw, showCategory = 5) )
      print( enrichplot::cnetplot(enrich_raw, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
      dev.off()
    }
    tryCatch(
      myTempFunction1234gty(),
      error = function(err){"myTempFunction1234gty_444888"}
    )
    svg(filename=paste(MyFolder, "/", MyFileName, "_3_enrich_raw.svg", sep=""),  width=8, height=5)
        print( enrichplot::dotplot(enrich_raw, showCategory = 30) )
    dev.off()
    enrich_raw_dataFrame     <- as.data.frame( enrich_raw ) 
    if(nrow(enrich_raw_dataFrame) > 1) {
      enrich_raw_dataFrame$GeneRatio <- sapply(enrich_raw$GeneRatio, function(x) eval(parse(text=x))) * 100
      enrich_raw_dataFrame$pvalue    <- -log10(enrich_raw_dataFrame$pvalue)
      enrich_raw_dataFrame$p.adjust  <- -log10(enrich_raw_dataFrame$p.adjust)
      enrich_raw_dataFrame$qvalue    <- -log10(enrich_raw_dataFrame$qvalue)
      enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[1:20, -8]
      if( nrow(enrich_raw_dataFrame) <20 ) {enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[, -8]}
      enrich_raw_dataFrame_selected <- transform(enrich_raw_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
      myMidValue <- ( max(enrich_raw_dataFrame_selected$GeneRatio) + min(enrich_raw_dataFrame_selected$GeneRatio)  )/2
      FigureTemp1 <- ggplot( data = enrich_raw_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
        geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
        scale_y_discrete( breaks=enrich_raw_dataFrame_selected$ID, labels=rev(enrich_raw_dataFrame_selected$ID) ) + 
        xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
      MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_3_enrich_raw", sep=""),  height1=6,  width1=5)
    }
    
    myTempFunction98769 <- function() {
      dev.off()
      dev.off()
      dev.off()
      dev.off()
      dev.off()
    }
    tryCatch(
      myTempFunction98769(),
      error = function(err){"myTempFunction98769_abhfd"}
    )

      
  }
  tryCatch(
    myTempFunction3(),
    error = function(err){"MyPeaksDistribution_OneGroup_1_g:error_3."}
  )
}



 
#############
MyEnrichDGN_OneGroup_1_g <- function(MyMaxDis=900000000,  MyPeaksAnno, MyFolder,  MyFileName ) {  ## Enrichment analysis based on the DisGeNET (http://www.disgenet.org/)
  myTempFunction3 <- function() {    
  peakAnno2 <- as.data.frame(MyPeaksAnno)
  peakAnno2 <- peakAnno2[ abs(peakAnno2$distanceToTSS) < MyMaxDis, ]
  write.table(x=as.data.frame(peakAnno2), file = paste(MyFolder, "Genomic_regions_annotation_within-XXXkb.txt", sep="/"),       
              append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
  
  
  ## 只有这里需要改动
  enrich_raw <- enrichDGN(gene=as.data.frame(peakAnno2)$geneId,   pvalueCutoff = 0.05,  pAdjustMethod = "BH", 
                          qvalueCutoff = 0.2,  readable = TRUE  )
  
  #
  write.table(x= enrich_raw , file = paste(MyFolder, "/", MyFileName, "1_enrich_raw.txt", sep=""), 
              append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
  myTempFunction1234gty <- function() {
    pdf(file=paste(MyFolder, "/", MyFileName, "_2_enrich_raw.pdf", sep=""),  width=8, height=5)
    print( barplot(enrich_raw, drop=TRUE, showCategory=20, x = "Count",     color = "p.adjust", font.size = 12, title = "")  )
    print( barplot(enrich_raw, drop=TRUE, showCategory=20, x = "GeneRatio", color = "p.adjust", font.size = 12, title = "")  )
    print( enrichplot::dotplot(enrich_raw, showCategory = 20 ) )
    print( enrichplot::Map(enrich_raw, showCategory = 20) )
    print( enrichplot::cnetplot(enrich_raw, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
    print( enrichplot::dotplot(enrich_raw, showCategory = 10) )
    print( enrichplot::Map(enrich_raw, showCategory = 10) )
    print( enrichplot::cnetplot(enrich_raw, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
    print( enrichplot::dotplot(enrich_raw, showCategory = 5) )
    print( enrichplot::Map(enrich_raw, showCategory = 5) )
    print( enrichplot::cnetplot(enrich_raw, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
    dev.off()
  }
  tryCatch(
    myTempFunction1234gty(),
    error = function(err){"myTempFunction1234gty_444888"}
  )
  svg(filename=paste(MyFolder, "/", MyFileName, "_3_enrich_raw.svg", sep=""),  width=8, height=5)
  print( enrichplot::dotplot(enrich_raw, showCategory = 30) )
  dev.off()
  enrich_raw_dataFrame     <- as.data.frame( enrich_raw ) 
  if(nrow(enrich_raw_dataFrame) > 1) {
    enrich_raw_dataFrame$GeneRatio <- sapply(enrich_raw$GeneRatio, function(x) eval(parse(text=x))) * 100
    enrich_raw_dataFrame$pvalue    <- -log10(enrich_raw_dataFrame$pvalue)
    enrich_raw_dataFrame$p.adjust  <- -log10(enrich_raw_dataFrame$p.adjust)
    enrich_raw_dataFrame$qvalue    <- -log10(enrich_raw_dataFrame$qvalue)
    enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[1:20, -8]
    if( nrow(enrich_raw_dataFrame) <20 ) {enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[, -8]}
    enrich_raw_dataFrame_selected <- transform(enrich_raw_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
    myMidValue <- ( max(enrich_raw_dataFrame_selected$GeneRatio) + min(enrich_raw_dataFrame_selected$GeneRatio)  )/2
    FigureTemp1 <- ggplot( data = enrich_raw_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
      geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
      scale_y_discrete( breaks=enrich_raw_dataFrame_selected$ID, labels=rev(enrich_raw_dataFrame_selected$ID) ) + 
      xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
    MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_3_enrich_raw", sep=""),  height1=6,  width1=5)
  }
  
  myTempFunction98769 <- function() {
    dev.off()
    dev.off()
    dev.off()
    dev.off()
    dev.off()
  }
  tryCatch(
    myTempFunction98769(),
    error = function(err){"myTempFunction98769_abhfd"}
  )
  
  
      
  }
  tryCatch(
    myTempFunction3(),
    error = function(err){"MyPeaksDistribution_OneGroup_1_g:error_3."}
  )
}




 
#############
MyEnrichDO_OneGroup_1_g <- function(MyMaxDis=900000000,  MyPeaksAnno, MyFolder,  MyFileName ) {  ## DO Enrichment Analysis
  myTempFunction3 <- function() {    
  peakAnno2 <- as.data.frame(MyPeaksAnno)
  peakAnno2 <- peakAnno2[ abs(peakAnno2$distanceToTSS) < MyMaxDis, ]
  write.table(x=as.data.frame(peakAnno2), file = paste(MyFolder, "Genomic_regions_annotation_within-XXXkb.txt", sep="/"),       
              append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
  
  
  ## 只有这里需要改动
  enrich_raw <- enrichDO(gene=as.data.frame(peakAnno2)$geneId,   ont = "DO", pvalueCutoff = 0.05,  pAdjustMethod = "BH", 
                         qvalueCutoff = 0.2,  readable = TRUE  )
  
  #
  write.table(x= enrich_raw , file = paste(MyFolder, "/", MyFileName, "1_enrich_raw.txt", sep=""), 
              append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
  myTempFunction1234gty <- function() {
    pdf(file=paste(MyFolder, "/", MyFileName, "_2_enrich_raw.pdf", sep=""),  width=8, height=5)
    print( barplot(enrich_raw, drop=TRUE, showCategory=20, x = "Count",     color = "p.adjust", font.size = 12, title = "")  )
    print( barplot(enrich_raw, drop=TRUE, showCategory=20, x = "GeneRatio", color = "p.adjust", font.size = 12, title = "")  )
    print( enrichplot::dotplot(enrich_raw, showCategory = 20 ) )
    print( enrichplot::emapplot(enrich_raw, showCategory = 20) )
    print( enrichplot::cnetplot(enrich_raw, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
    print( enrichplot::dotplot(enrich_raw, showCategory = 10) )
    print( enrichplot::emapplot(enrich_raw, showCategory = 10) )
    print( enrichplot::cnetplot(enrich_raw, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
    print( enrichplot::dotplot(enrich_raw, showCategory = 5) )
    print( enrichplot::emapplot(enrich_raw, showCategory = 5) )
    print( enrichplot::cnetplot(enrich_raw, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
    dev.off()
  }
  tryCatch(
    myTempFunction1234gty(),
    error = function(err){"myTempFunction1234gty_444888"}
  )
  svg(filename=paste(MyFolder, "/", MyFileName, "_3_enrich_raw.svg", sep=""),  width=8, height=5)
  print( enrichplot::dotplot(enrich_raw, showCategory = 30) )
  dev.off()
  enrich_raw_dataFrame     <- as.data.frame( enrich_raw ) 
  if(nrow(enrich_raw_dataFrame) > 1) {
    enrich_raw_dataFrame$GeneRatio <- sapply(enrich_raw$GeneRatio, function(x) eval(parse(text=x))) * 100
    enrich_raw_dataFrame$pvalue    <- -log10(enrich_raw_dataFrame$pvalue)
    enrich_raw_dataFrame$p.adjust  <- -log10(enrich_raw_dataFrame$p.adjust)
    enrich_raw_dataFrame$qvalue    <- -log10(enrich_raw_dataFrame$qvalue)
    enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[1:20, -8]
    if( nrow(enrich_raw_dataFrame) <20 ) {enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[, -8]}
    enrich_raw_dataFrame_selected <- transform(enrich_raw_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
    myMidValue <- ( max(enrich_raw_dataFrame_selected$GeneRatio) + min(enrich_raw_dataFrame_selected$GeneRatio)  )/2
    FigureTemp1 <- ggplot( data = enrich_raw_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
      geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
      scale_y_discrete( breaks=enrich_raw_dataFrame_selected$ID, labels=rev(enrich_raw_dataFrame_selected$ID) ) + 
      xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
    MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_3_enrich_raw", sep=""),  height1=6,  width1=5)
  }
  
  myTempFunction98769 <- function() {
    dev.off()
    dev.off()
    dev.off()
    dev.off()
    dev.off()
  }
  tryCatch(
    myTempFunction98769(),
    error = function(err){"myTempFunction98769_abhfd"}
  )
  
       
  }
  tryCatch(
    myTempFunction3(),
    error = function(err){"MyPeaksDistribution_OneGroup_1_g:error_3."}
  ) 
}




#############
MyEnrichNCG_OneGroup_1_g <- function(MyMaxDis=900000000,  MyPeaksAnno, MyFolder,  MyFileName ) {  ## Enrichment analysis based on the Network of Cancer Genes database (http://ncg.kcl.ac.uk/)
  myTempFunction3 <- function() {    
  peakAnno2 <- as.data.frame(MyPeaksAnno)
  peakAnno2 <- peakAnno2[ abs(peakAnno2$distanceToTSS) < MyMaxDis, ]
  write.table(x=as.data.frame(peakAnno2), file = paste(MyFolder, "Genomic_regions_annotation_within-XXXkb.txt", sep="/"),       
              append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
  
  
  ## 只有这里需要改动
  enrich_raw <- enrichNCG(gene=as.data.frame(peakAnno2)$geneId,   pvalueCutoff = 0.05,  pAdjustMethod = "BH", 
                          qvalueCutoff = 0.2,    readable = TRUE  )
  
  #
  write.table(x= enrich_raw , file = paste(MyFolder, "/", MyFileName, "1_enrich_raw.txt", sep=""), 
              append = FALSE,  quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE )
  myTempFunction1234gty <- function() {
    pdf(file=paste(MyFolder, "/", MyFileName, "_2_enrich_raw.pdf", sep=""),  width=8, height=5)
    print( barplot(enrich_raw, drop=TRUE, showCategory=20, x = "Count",     color = "p.adjust", font.size = 12, title = "")  )
    print( barplot(enrich_raw, drop=TRUE, showCategory=20, x = "GeneRatio", color = "p.adjust", font.size = 12, title = "")  )
    print( enrichplot::dotplot(enrich_raw, showCategory = 20 ) )
    print( enrichplot::Map(enrich_raw, showCategory = 20) )
    print( enrichplot::cnetplot(enrich_raw, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
    print( enrichplot::dotplot(enrich_raw, showCategory = 10) )
    print( enrichplot::Map(enrich_raw, showCategory = 10) )
    print( enrichplot::cnetplot(enrich_raw, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
    print( enrichplot::dotplot(enrich_raw, showCategory = 5) )
    print( enrichplot::Map(enrich_raw, showCategory = 5) )
    print( enrichplot::cnetplot(enrich_raw, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
    print( enrichMap(enrich_raw, n = 10, fixed = TRUE, vertex.label.font = 1) )
    print( enrichMap(enrich_raw, n = 20, fixed = TRUE, vertex.label.font = 1) )
    dev.off()
  }
  tryCatch(
    myTempFunction1234gty(),
    error = function(err){"myTempFunction1234gty_444888"}
  )
  svg(filename=paste(MyFolder, "/", MyFileName, "_3_enrich_raw.svg", sep=""),  width=8, height=5)
  print( enrichplot::dotplot(enrich_raw, showCategory = 30) )
  dev.off()
  enrich_raw_dataFrame     <- as.data.frame( enrich_raw ) 
  if(nrow(enrich_raw_dataFrame) > 1) {
    enrich_raw_dataFrame$GeneRatio <- sapply(enrich_raw$GeneRatio, function(x) eval(parse(text=x))) * 100
    enrich_raw_dataFrame$pvalue    <- -log10(enrich_raw_dataFrame$pvalue)
    enrich_raw_dataFrame$p.adjust  <- -log10(enrich_raw_dataFrame$p.adjust)
    enrich_raw_dataFrame$qvalue    <- -log10(enrich_raw_dataFrame$qvalue)
    enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[1:20, -8]
    if( nrow(enrich_raw_dataFrame) <20 ) {enrich_raw_dataFrame_selected <- enrich_raw_dataFrame[, -8]}
    enrich_raw_dataFrame_selected <- transform(enrich_raw_dataFrame_selected, ID= rev( factor(ID, levels=unique(ID))) )    
    myMidValue <- ( max(enrich_raw_dataFrame_selected$GeneRatio) + min(enrich_raw_dataFrame_selected$GeneRatio)  )/2
    FigureTemp1 <- ggplot( data = enrich_raw_dataFrame_selected, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
      geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
      scale_y_discrete( breaks=enrich_raw_dataFrame_selected$ID, labels=rev(enrich_raw_dataFrame_selected$ID) ) + 
      xlab("-log10(ajusted p-value)") +   ylab("Ontology terms") +   ggtitle("Ontology terms") +  theme_bw()
    MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_3_enrich_raw", sep=""),  height1=6,  width1=5)
  }
  
  myTempFunction98769 <- function() {
    dev.off()
    dev.off()
    dev.off()
    dev.off()
    dev.off()
  }
  tryCatch(
    myTempFunction98769(),
    error = function(err){"myTempFunction98769_abhfd"}
  )
  
        
  }
  tryCatch(
    myTempFunction3(),
    error = function(err){"MyPeaksDistribution_OneGroup_1_g:error_3."}
  )
}

 
######################################################################################################################################################





######################################################################################################################################################
MyAnnotation_OneGroup_1_g <- function(myFile_t1, myPath_t1) {
  if( ! file.exists( myPath_t1 )  ) { dir.create(myPath_t1,  recursive = TRUE)  }
  sink(file = paste(myPath_t1, "1_fileName.txt", sep="/") )
      print( myFile_t1 )
  sink()
  myPeak_t1 = readPeakFile( myFile_t1 )
  
  myPeak_t1_data.frame = as.data.frame(myPeak_t1)
  sink(file = paste(myPath_t1, "2_dimension-of-inputFile.txt", sep="/") )
      print( dim(myPeak_t1_data.frame) )
  sink()
  
  write.table(x=myPeak_t1_data.frame, file = paste(myPath_t1, "3_Genomic-regions.txt", sep="/"),       
              append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".",  row.names = FALSE, col.names = TRUE )

  #MyPeaksAnno_OneGroup_1_g(myPeak_t1 = myPeak_t1, myPath_t1 = paste(myPath_t1,  "1_GenomicRegions-on-chromosomes",  sep="/") )
  #MyPeaksSignal_OneGroup_1_g(myPeak_t1 = myPeak_t1, myPath_t1 = paste(myPath_t1,  "2_GenomicRegions-on-TSSs-10kb",  sep="/"),  up1=upstream_g, down1=downstream_g)
  MyPeaksDistribution_OneGroup_1_g(myPeak_t1 = myPeak_t1, myPath_t1 = paste(myPath_t1,  "3_GenomicRegions-distribution",  sep="/"),  up1=upstream_g, down1=downstream_g)
  #MyPeaksDistribution_OneGroup_1_g(myPeak_t1 = myPeak_t1, myPath_t1 = paste(myPath_t1,  "4_GenomicRegions-distribution_x2",  sep="/"),  up1=upstream_g+1000, down1=downstream_g+1000)
  
  peakAnno <- annotatePeak( myPeak_t1, tssRegion=c(-upstream_g, downstream_g),  TxDb=my_txdb_g, annoDb=my_orgdb_g, sameStrand = stranded_g , overlap="all" )  
  write.table(x=as.data.frame(peakAnno), file = paste(myPath_t1, "4_Genomic_regions_annotation.txt", sep="/"),      
              append = FALSE, quote = FALSE, sep = "\t",   eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE )
  
  
  if(fast_g == FALSE) {  
  my_folder_5A <- paste(myPath_t1,  "5A_GeneOntology-BiologicalProcess-All",  sep="/") 
  if( ! file.exists(my_folder_5A)  ) { dir.create(my_folder_5A,  recursive = TRUE)  }
  MyGeneOntology_BP_OneGroup_1_g(MyMaxDis=900000000, MyPeaksAnno=peakAnno, MyFolder=my_folder_5A,  MyFileName="5A_all-NearestGenes_")
  
  my_folder_5B <- paste(myPath_t1,  "5B_GeneOntology-BiologicalProcess-within100kb",  sep="/") 
  if( ! file.exists(my_folder_5B)  ) { dir.create(my_folder_5B,  recursive = TRUE)  }
  MyGeneOntology_BP_OneGroup_1_g(MyMaxDis=100000, MyPeaksAnno=peakAnno, MyFolder=my_folder_5B,  MyFileName="5B_100kb-NearestGenes_")
  
  my_folder_5C <- paste(myPath_t1,  "5C_GeneOntology-BiologicalProcess-within50kb",  sep="/") 
  if( ! file.exists(my_folder_5C)  ) { dir.create(my_folder_5C,  recursive = TRUE)  }
  MyGeneOntology_BP_OneGroup_1_g(MyMaxDis=50000, MyPeaksAnno=peakAnno, MyFolder=my_folder_5C,  MyFileName="5C_50kb-NearestGenes_")
  
  my_folder_5D <- paste(myPath_t1,  "5D_GeneOntology-BiologicalProcess-within20kb",  sep="/") 
  if( ! file.exists(my_folder_5D)  ) { dir.create(my_folder_5D,  recursive = TRUE)  }
  MyGeneOntology_BP_OneGroup_1_g(MyMaxDis=20000, MyPeaksAnno=peakAnno, MyFolder=my_folder_5D,  MyFileName="5D_20kb-NearestGenes_")
  
  my_folder_5E <- paste(myPath_t1,  "5E_GeneOntology-BiologicalProcess-within5kb",  sep="/") 
  if( ! file.exists(my_folder_5E)  ) { dir.create(my_folder_5E,  recursive = TRUE)  }
  MyGeneOntology_BP_OneGroup_1_g(MyMaxDis=5000, MyPeaksAnno=peakAnno, MyFolder=my_folder_5E,  MyFileName="5E_5kb-NearestGenes_")
  
  
  
  my_folder_6A <- paste(myPath_t1,  "6A_GeneOntology-MolecularFunction-All",  sep="/") 
  if( ! file.exists(my_folder_6A)  ) { dir.create(my_folder_6A,  recursive = TRUE)  }
  MyGeneOntology_MF_OneGroup_1_g(MyMaxDis=900000000, MyPeaksAnno=peakAnno, MyFolder=my_folder_6A,  MyFileName="6A_all-NearestGenes_")
  
  my_folder_6B <- paste(myPath_t1,  "6B_GeneOntology-MolecularFunction-within100kb",  sep="/") 
  if( ! file.exists(my_folder_6B)  ) { dir.create(my_folder_6B,  recursive = TRUE)  }
  MyGeneOntology_MF_OneGroup_1_g(MyMaxDis=100000, MyPeaksAnno=peakAnno, MyFolder=my_folder_6B,  MyFileName="6B_100kb-NearestGenes_")
  
  my_folder_6C <- paste(myPath_t1,  "6C_GeneOntology-MolecularFunction-within50kb",  sep="/") 
  if( ! file.exists(my_folder_6C)  ) { dir.create(my_folder_6C,  recursive = TRUE)  }
  MyGeneOntology_MF_OneGroup_1_g(MyMaxDis=50000, MyPeaksAnno=peakAnno, MyFolder=my_folder_6C,  MyFileName="6C_50kb-NearestGenes_")
  
  my_folder_6D <- paste(myPath_t1,  "6D_GeneOntology-MolecularFunction-within20kb",  sep="/") 
  if( ! file.exists(my_folder_6D)  ) { dir.create(my_folder_6D,  recursive = TRUE)  }
  MyGeneOntology_MF_OneGroup_1_g(MyMaxDis=20000, MyPeaksAnno=peakAnno, MyFolder=my_folder_6D,  MyFileName="6D_20kb-NearestGenes_")
  
  my_folder_6E <- paste(myPath_t1,  "6E_GeneOntology-MolecularFunction-within5kb",  sep="/") 
  if( ! file.exists(my_folder_6E)  ) { dir.create(my_folder_6E,  recursive = TRUE)  }
  MyGeneOntology_MF_OneGroup_1_g(MyMaxDis=5000, MyPeaksAnno=peakAnno, MyFolder=my_folder_6E,  MyFileName="6E_5kb-NearestGenes_")
  
  
  
  
  my_folder_7A <- paste(myPath_t1,  "7A_GeneOntology-CellularComponent-All",  sep="/") 
  if( ! file.exists(my_folder_7A)  ) { dir.create(my_folder_7A,  recursive = TRUE)  }
  MyGeneOntology_CC_OneGroup_1_g(MyMaxDis=900000000, MyPeaksAnno=peakAnno, MyFolder=my_folder_7A,  MyFileName="7A_all-NearestGenes_")
  
  my_folder_7B <- paste(myPath_t1,  "7B_GeneOntology-CellularComponent-within100kb",  sep="/") 
  if( ! file.exists(my_folder_7B)  ) { dir.create(my_folder_7B,  recursive = TRUE)  }
  MyGeneOntology_CC_OneGroup_1_g(MyMaxDis=100000, MyPeaksAnno=peakAnno, MyFolder=my_folder_7B,  MyFileName="7B_100kb-NearestGenes_")
  
  my_folder_7C <- paste(myPath_t1,  "7C_GeneOntology-CellularComponent-within50kb",  sep="/") 
  if( ! file.exists(my_folder_7C)  ) { dir.create(my_folder_7C,  recursive = TRUE)  }
  MyGeneOntology_CC_OneGroup_1_g(MyMaxDis=50000, MyPeaksAnno=peakAnno, MyFolder=my_folder_7C,  MyFileName="7C_50kb-NearestGenes_")
  
  my_folder_7D <- paste(myPath_t1,  "7D_GeneOntology-CellularComponent-within20kb",  sep="/") 
  if( ! file.exists(my_folder_7D)  ) { dir.create(my_folder_7D,  recursive = TRUE)  }
  MyGeneOntology_CC_OneGroup_1_g(MyMaxDis=20000, MyPeaksAnno=peakAnno, MyFolder=my_folder_7D,  MyFileName="7D_20kb-NearestGenes_")
  
  my_folder_7E <- paste(myPath_t1,  "7E_GeneOntology-CellularComponent-within5kb",  sep="/") 
  if( ! file.exists(my_folder_7E)  ) { dir.create(my_folder_7E,  recursive = TRUE)  }
  MyGeneOntology_CC_OneGroup_1_g(MyMaxDis=5000, MyPeaksAnno=peakAnno, MyFolder=my_folder_7E,  MyFileName="7E_5kb-NearestGenes_")
  
  
  
  
  my_folder_8A <- paste(myPath_t1,  "8A_Reactome-All",  sep="/") 
  if( ! file.exists(my_folder_8A)  ) { dir.create(my_folder_8A,  recursive = TRUE)  }
  MyReactome_OneGroup_1_g(MyMaxDis=900000000, MyPeaksAnno=peakAnno, MyFolder=my_folder_8A,  MyFileName="8A_all-NearestGenes_")
  
  my_folder_8B <- paste(myPath_t1,  "8B_Reactome-within100kb",  sep="/") 
  if( ! file.exists(my_folder_8B)  ) { dir.create(my_folder_8B,  recursive = TRUE)  }
  MyReactome_OneGroup_1_g(MyMaxDis=100000, MyPeaksAnno=peakAnno, MyFolder=my_folder_8B,  MyFileName="8B_100kb-NearestGenes_")
  
  my_folder_8C <- paste(myPath_t1,  "8C_Reactome-within50kb",  sep="/") 
  if( ! file.exists(my_folder_8C)  ) { dir.create(my_folder_8C,  recursive = TRUE)  }
  MyReactome_OneGroup_1_g(MyMaxDis=50000, MyPeaksAnno=peakAnno, MyFolder=my_folder_8C,  MyFileName="8C_50kb-NearestGenes_")
  
  my_folder_8D <- paste(myPath_t1,  "8D_Reactome-within20kb",  sep="/") 
  if( ! file.exists(my_folder_8D)  ) { dir.create(my_folder_8D,  recursive = TRUE)  }
  MyReactome_OneGroup_1_g(MyMaxDis=20000, MyPeaksAnno=peakAnno, MyFolder=my_folder_8D,  MyFileName="8D_20kb-NearestGenes_")
  
  my_folder_8E <- paste(myPath_t1,  "8E_Reactome-within5kb",  sep="/") 
  if( ! file.exists(my_folder_8E)  ) { dir.create(my_folder_8E,  recursive = TRUE)  }
  MyReactome_OneGroup_1_g(MyMaxDis=5000, MyPeaksAnno=peakAnno, MyFolder=my_folder_8E,  MyFileName="8E_5kb-NearestGenes_")
  
  
  
  
  my_folder_9A <- paste(myPath_t1,  "9A_KEGG-All",  sep="/") 
  if( ! file.exists(my_folder_9A)  ) { dir.create(my_folder_9A,  recursive = TRUE)  }
  MyKEGG_OneGroup_1_g(MyMaxDis=900000000, MyPeaksAnno=peakAnno, MyFolder=my_folder_9A,  MyFileName="9A_all-NearestGenes_")
  
  my_folder_9B <- paste(myPath_t1,  "9B_KEGG-within100kb",  sep="/") 
  if( ! file.exists(my_folder_9B)  ) { dir.create(my_folder_9B,  recursive = TRUE)  }
  MyKEGG_OneGroup_1_g(MyMaxDis=100000, MyPeaksAnno=peakAnno, MyFolder=my_folder_9B,  MyFileName="9B_100kb-NearestGenes_")
  
  my_folder_9C <- paste(myPath_t1,  "9C_KEGG-within50kb",  sep="/") 
  if( ! file.exists(my_folder_9C)  ) { dir.create(my_folder_9C,  recursive = TRUE)  }
  MyKEGG_OneGroup_1_g(MyMaxDis=50000, MyPeaksAnno=peakAnno, MyFolder=my_folder_9C,  MyFileName="9C_50kb-NearestGenes_")
  
  my_folder_9D <- paste(myPath_t1,  "9D_KEGG-within20kb",  sep="/") 
  if( ! file.exists(my_folder_9D)  ) { dir.create(my_folder_9D,  recursive = TRUE)  }
  MyKEGG_OneGroup_1_g(MyMaxDis=20000, MyPeaksAnno=peakAnno, MyFolder=my_folder_9D,  MyFileName="9D_20kb-NearestGenes_")
  
  my_folder_9E <- paste(myPath_t1,  "9E_KEGG-within5kb",  sep="/") 
  if( ! file.exists(my_folder_9E)  ) { dir.create(my_folder_9E,  recursive = TRUE)  }
  MyKEGG_OneGroup_1_g(MyMaxDis=5000, MyPeaksAnno=peakAnno, MyFolder=my_folder_9E,  MyFileName="9E_5kb-NearestGenes_")
  
  
  
  my_folder_10A <- paste(myPath_t1,  "10A_ModuleKEGG-All",  sep="/") 
  if( ! file.exists(my_folder_10A)  ) { dir.create(my_folder_10A,  recursive = TRUE)  }
  MyModuleKEGG_OneGroup_1_g(MyMaxDis=900000000, MyPeaksAnno=peakAnno, MyFolder=my_folder_10A,  MyFileName="10A_all-NearestGenes_")
  
  my_folder_10B <- paste(myPath_t1,  "10B_ModuleKEGG-within100kb",  sep="/") 
  if( ! file.exists(my_folder_10B)  ) { dir.create(my_folder_10B,  recursive = TRUE)  }
  MyModuleKEGG_OneGroup_1_g(MyMaxDis=100000, MyPeaksAnno=peakAnno, MyFolder=my_folder_10B,  MyFileName="10B_100kb-NearestGenes_")
  
  my_folder_10C <- paste(myPath_t1,  "10C_ModuleKEGG-within50kb",  sep="/") 
  if( ! file.exists(my_folder_10C)  ) { dir.create(my_folder_10C,  recursive = TRUE)  }
  MyModuleKEGG_OneGroup_1_g(MyMaxDis=50000, MyPeaksAnno=peakAnno, MyFolder=my_folder_10C,  MyFileName="10C_50kb-NearestGenes_")
  
  my_folder_10D <- paste(myPath_t1,  "10D_ModuleKEGG-within20kb",  sep="/") 
  if( ! file.exists(my_folder_10D)  ) { dir.create(my_folder_10D,  recursive = TRUE)  }
  MyModuleKEGG_OneGroup_1_g(MyMaxDis=20000, MyPeaksAnno=peakAnno, MyFolder=my_folder_10D,  MyFileName="10D_20kb-NearestGenes_")
  
  my_folder_10E <- paste(myPath_t1,  "10E_ModuleKEGG-within5kb",  sep="/") 
  if( ! file.exists(my_folder_10E)  ) { dir.create(my_folder_10E,  recursive = TRUE)  }
  MyModuleKEGG_OneGroup_1_g(MyMaxDis=5000, MyPeaksAnno=peakAnno, MyFolder=my_folder_10E,  MyFileName="10E_5kb-NearestGenes_")
  
  
  
  my_folder_11A <- paste(myPath_t1,  "11A_DiseaseEnrichmentOnDisGeNET-All",  sep="/") 
  if( ! file.exists(my_folder_11A)  ) { dir.create(my_folder_11A,  recursive = TRUE)  }
  MyEnrichDGN_OneGroup_1_g(MyMaxDis=900000000, MyPeaksAnno=peakAnno, MyFolder=my_folder_11A,  MyFileName="11A_all-NearestGenes_")
  
  my_folder_11B <- paste(myPath_t1,  "11B_DiseaseEnrichmentOnDisGeNET-within100kb",  sep="/") 
  if( ! file.exists(my_folder_11B)  ) { dir.create(my_folder_11B,  recursive = TRUE)  }
  MyEnrichDGN_OneGroup_1_g(MyMaxDis=100000, MyPeaksAnno=peakAnno, MyFolder=my_folder_11B,  MyFileName="11B_100kb-NearestGenes_")
  
  my_folder_11C <- paste(myPath_t1,  "11C_DiseaseEnrichmentOnDisGeNET-within50kb",  sep="/") 
  if( ! file.exists(my_folder_11C)  ) { dir.create(my_folder_11C,  recursive = TRUE)  }
  MyEnrichDGN_OneGroup_1_g(MyMaxDis=50000, MyPeaksAnno=peakAnno, MyFolder=my_folder_11C,  MyFileName="11C_50kb-NearestGenes_")
  
  my_folder_11D <- paste(myPath_t1,  "11D_DiseaseEnrichmentOnDisGeNET-within20kb",  sep="/") 
  if( ! file.exists(my_folder_11D)  ) { dir.create(my_folder_11D,  recursive = TRUE)  }
  MyEnrichDGN_OneGroup_1_g(MyMaxDis=20000, MyPeaksAnno=peakAnno, MyFolder=my_folder_11D,  MyFileName="11D_20kb-NearestGenes_")
  
  my_folder_11E <- paste(myPath_t1,  "11E_DiseaseEnrichmentOnDisGeNET-within5kb",  sep="/") 
  if( ! file.exists(my_folder_11E)  ) { dir.create(my_folder_11E,  recursive = TRUE)  }
  MyEnrichDGN_OneGroup_1_g(MyMaxDis=5000, MyPeaksAnno=peakAnno, MyFolder=my_folder_11E,  MyFileName="11E_5kb-NearestGenes_")
  
  
  
  
  my_folder_12A <- paste(myPath_t1,  "12A_DiseaseEnrichmentOnDO-All",  sep="/") 
  if( ! file.exists(my_folder_12A)  ) { dir.create(my_folder_12A,  recursive = TRUE)  }
  MyEnrichDO_OneGroup_1_g(MyMaxDis=900000000, MyPeaksAnno=peakAnno, MyFolder=my_folder_12A,  MyFileName="12A_all-NearestGenes_")
  
  my_folder_12B <- paste(myPath_t1,  "12B_DiseaseEnrichmentOnDO-within100kb",  sep="/") 
  if( ! file.exists(my_folder_12B)  ) { dir.create(my_folder_12B,  recursive = TRUE)  }
  MyEnrichDO_OneGroup_1_g(MyMaxDis=100000, MyPeaksAnno=peakAnno, MyFolder=my_folder_12B,  MyFileName="12B_100kb-NearestGenes_")
  
  my_folder_12C <- paste(myPath_t1,  "12C_DiseaseEnrichmentOnDO-within50kb",  sep="/") 
  if( ! file.exists(my_folder_12C)  ) { dir.create(my_folder_12C,  recursive = TRUE)  }
  MyEnrichDO_OneGroup_1_g(MyMaxDis=50000, MyPeaksAnno=peakAnno, MyFolder=my_folder_12C,  MyFileName="12C_50kb-NearestGenes_")
  
  my_folder_12D <- paste(myPath_t1,  "12D_DiseaseEnrichmentOnDO-within20kb",  sep="/") 
  if( ! file.exists(my_folder_12D)  ) { dir.create(my_folder_12D,  recursive = TRUE)  }
  MyEnrichDO_OneGroup_1_g(MyMaxDis=20000, MyPeaksAnno=peakAnno, MyFolder=my_folder_12D,  MyFileName="12D_20kb-NearestGenes_")
  
  my_folder_12E <- paste(myPath_t1,  "12E_DiseaseEnrichmentOnDO-within5kb",  sep="/") 
  if( ! file.exists(my_folder_12E)  ) { dir.create(my_folder_12E,  recursive = TRUE)  }
  MyEnrichDO_OneGroup_1_g(MyMaxDis=5000, MyPeaksAnno=peakAnno, MyFolder=my_folder_12E,  MyFileName="12E_5kb-NearestGenes_")
  
  
  
  my_folder_13A <- paste(myPath_t1,  "13A_DiseaseEnrichmentOnTheNetworkOfCancerGenesDatabase-All",  sep="/") 
  if( ! file.exists(my_folder_13A)  ) { dir.create(my_folder_13A,  recursive = TRUE)  }
  MyEnrichNCG_OneGroup_1_g(MyMaxDis=900000000, MyPeaksAnno=peakAnno, MyFolder=my_folder_13A,  MyFileName="13A_all-NearestGenes_")
  
  my_folder_13B <- paste(myPath_t1,  "13B_DiseaseEnrichmentOnTheNetworkOfCancerGenesDatabase-within100kb",  sep="/") 
  if( ! file.exists(my_folder_13B)  ) { dir.create(my_folder_13B,  recursive = TRUE)  }
  MyEnrichNCG_OneGroup_1_g(MyMaxDis=100000, MyPeaksAnno=peakAnno, MyFolder=my_folder_13B,  MyFileName="13B_100kb-NearestGenes_")
  
  my_folder_13C <- paste(myPath_t1,  "13C_DiseaseEnrichmentOnTheNetworkOfCancerGenesDatabase-within50kb",  sep="/") 
  if( ! file.exists(my_folder_13C)  ) { dir.create(my_folder_13C,  recursive = TRUE)  }
  MyEnrichNCG_OneGroup_1_g(MyMaxDis=50000, MyPeaksAnno=peakAnno, MyFolder=my_folder_13C,  MyFileName="13C_50kb-NearestGenes_")
  
  my_folder_13D <- paste(myPath_t1,  "13D_DiseaseEnrichmentOnTheNetworkOfCancerGenesDatabase-within20kb",  sep="/") 
  if( ! file.exists(my_folder_13D)  ) { dir.create(my_folder_13D,  recursive = TRUE)  }
  MyEnrichNCG_OneGroup_1_g(MyMaxDis=20000, MyPeaksAnno=peakAnno, MyFolder=my_folder_13D,  MyFileName="13D_20kb-NearestGenes_")
  
  my_folder_13E <- paste(myPath_t1,  "13E_DiseaseEnrichmentOnTheNetworkOfCancerGenesDatabase-within5kb",  sep="/") 
  if( ! file.exists(my_folder_13E)  ) { dir.create(my_folder_13E,  recursive = TRUE)  }
  MyEnrichNCG_OneGroup_1_g(MyMaxDis=5000, MyPeaksAnno=peakAnno, MyFolder=my_folder_13E,  MyFileName="13E_5kb-NearestGenes_")
  }
   
  myTempFun_1 <- function() {
    dev.off()
    dev.off()
    dev.off()
    dev.off()
    dev.off()
  }
  tryCatch(
    myTempFun_1(),
    error = function(err){"myTempFun_1_YongPeng."}
  )
  
}
######################################################################################################################################################




 
######################################################################################################################################################
for(i in c(1:length(peakFiles)) ) { 
  myTempFunction3 <- function() {  
     ##MyAnnotation_OneGroup_1_g( myFile_t1 = peakFiles[i],  myPath_t1 = paste(outDir_g, "/File.", i, "._", peakFiles_onlyName[i],  sep="") )      
  }
  tryCatch(
    myTempFunction3(),
    error = function(err){"My_OneGroup_1_g:error_3."}
  )
}
######################################################################################################################################################










########################################### 
myPath_multi_1 = paste(outDir_g, "Multiple_Files",sep="/") 
if( ! file.exists( myPath_multi_1 )  ) { dir.create(myPath_multi_1,  recursive = TRUE)  }


peakAnnoList_1 <- lapply(peakFiles, annotatePeak,  TxDb=my_txdb_g, annoDb=my_orgdb_g, tssRegion=c(-upstream_g, downstream_g), verbose=FALSE, sameStrand = stranded_g , overlap="all")
names( peakAnnoList_1 ) <- peakFiles_onlyName

pdf( file=paste(myPath_multi_1, "Genomic_Locations.pdf", sep="/"),  width=10,  height=length(peakFiles_onlyName)/4 + 1 )
    plotAnnoBar(peakAnnoList_1)
    plotDistToTSS(peakAnnoList_1)
dev.off() 
 
sink( file = paste(myPath_multi_1, "Genomic_Locations.txt", sep="/")  )
print(peakAnnoList_1)  
sink()

########################################### 

 

 





