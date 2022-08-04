
## To run the script in command lines.
#################################################
suppressPackageStartupMessages( library(optparse) )      

getParameters_f <- function() {
  option_list_Local <- list(   # Options list with associated default value.   
      optparse::make_option(opt_str=c("-m", "--motif"),
			default="AAACA",
			type="character",   dest="motif",
			help="One 5-mer. [default: %default]."),
			
      optparse::make_option(opt_str=c("-r", "--rcmotif"),
			default="TGNTT",
			type="character",   dest="rcmotif",
			help="Reverse complememt of the 5-mer. [default: %default]."),

      optparse::make_option(opt_str=c("-o", "--outDir"),
			default="fiveLevels",
			type="character",   dest="outDir",
			help="Path or dir name for output files. [default: %default].")
  )

  # now parse the command line to check which option is given and get associated values
  parser_Local <- optparse::OptionParser(usage="usage: %prog [options]",
		option_list=option_list_Local, 
		description="fit_motif.R: linear regerssion for the correlation between mutation ratio and m6A fraction. August 2, 2020.",                             
		epilogue="For comments, bug reports etc..., please contact Yong Peng <yongp@outlook.com>"
  )

  opt_Local <- optparse::parse_args(parser_Local, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options
  return(opt_Local)   
}


opt_g = getParameters_f()
motif_g      <- opt_g$motif
rcmotif_g    <- opt_g$rcmotif
outDir_g     <- opt_g$outDir
# motif_g     <- "AAACA"
# rcmotif_g   <- "TGNTT"
# outDir_g    <- "fiveLevels"

print(motif_g)
print(rcmotif_g )
print(outDir_g)

outDir_g = paste(outDir_g, motif_g, sep="/")
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

cat("##########################\n")
col_names1 == col_names2
col_names1 == col_names3
col_names1 == col_names4
col_names1 == col_names5
cat("##########################\n")

## TGACA


DRACH_index = grepl(pattern = motif_g, x=col_names1, ignore.case = FALSE, perl = TRUE,  fixed = FALSE, useBytes = FALSE)  
col_names1[DRACH_index]
cat("##########################\n\n\n\n\n\n")

# motif_g     
# rcmotif_g
# outDir_g     

matrix1 = matrix1[, DRACH_index]
matrix2 = matrix2[, DRACH_index]
matrix3 = matrix3[, DRACH_index]
matrix4 = matrix4[, DRACH_index]
matrix5 = matrix5[, DRACH_index]

numbers1 = as.vector(matrix1)
numbers2 = as.vector(matrix2)
numbers3 = as.vector(matrix3)
numbers4 = as.vector(matrix4)
numbers5 = as.vector(matrix5)


ratio1 = numbers1/sum(numbers1)
ratio2 = numbers2/sum(numbers2)
ratio3 = numbers3/sum(numbers3)
ratio4 = numbers4/sum(numbers4)
ratio5 = numbers5/sum(numbers5)



allMutations = c(1-ratio1[1],   1-ratio2[1],   1-ratio3[1],    1-ratio5[1])
rawLevel = c(0.00, 0.25, 0.50,   1.00)





##################################################################################################################
library(ggplot2) 


MyTheme_1_g <- function(textSize1=14, hjust1=NULL, vjust1=NULL,  angle1=NULL) {    # "hjust=1, vjust=1, angle=30" for some boxplots.
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
    axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),    ## lines along axes (element_line; inherits from line). 坐标轴线
    axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),    ## line along x axis (element_line; inherits from axis.line)
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
    strip.text       = element_text(family="serif", face="bold", colour=NULL, size=rel(1.2), hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels (element_text; inherits from text)
    strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels along horizontal direction (element_text; inherits from strip.text)
    strip.text.y     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL)   	  ## facet labels along vertical direction (element_text; inherits from strip.text) 
  ) 
} 


MySaveGgplot2_1_g <- function(ggplot2Figure1,  path1, fileName1,  height1, width1) {
  SVG1 <- paste(path1,  "/",  "SVG",  sep = "",  collapse = NULL)
  PNG1 <- paste(path1,  "/",  "PNG",  sep = "",  collapse = NULL)
  PDF1 <- paste(path1,  "/",  "PDF",  sep = "",  collapse = NULL)
  EPS1 <- paste(path1,  "/",  "EPS",  sep = "",  collapse = NULL)
  if( ! file.exists(SVG1) ) { dir.create(SVG1) }
  #if( ! file.exists(PNG1) ) { dir.create(PNG1) }
  if( ! file.exists(PDF1) ) { dir.create(PDF1) }
  if( ! file.exists(EPS1) ) { dir.create(EPS1) }
  ggsave( filename = paste(SVG1,  "/",  fileName1,  ".svg",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 , limitsize = FALSE)
  #ggsave( filename = paste(PNG1,  "/",  fileName1,  ".png",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 , limitsize = FALSE )
  ggsave( filename = paste(PDF1,  "/",  fileName1,  ".pdf",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 , limitsize = FALSE )
  ggsave( filename = paste(EPS1,  "/",  fileName1,  ".eps",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200,   device=cairo_ps , limitsize = FALSE)         
}


MyRegression_1 <- function( xAxis2,  yAxis2,  file2 ) {  ## for any number of points
  regression_A <- lm(yAxis2 ~ xAxis2, model=TRUE, x=TRUE, y=TRUE, qr=TRUE)    
  print("####################################################################")
  print(regression_A)
  print(summary(regression_A))
  print("####################################################################")
  r2 = summary(regression_A)$adj.r.squared     ## names(summary(regression_A))
  pv = summary(regression_A)$coefficients[2,4]    
  b = coef(regression_A)[2]   ## b is slope.
  write.table(x=r2, file = paste(file2,  "_Adjusted-R-squared.txt", sep = ""), append = TRUE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
  write.table(x=pv, file = paste(file2,  "_Pvalue.txt",      sep = ""),        append = TRUE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
  write.table(x=b,  file = paste(file2,  "_slope.txt",      sep = ""),           append = TRUE, quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")                                                        
  return(b)
}


MyLmEquation_1 <- function(m) {  #### m is the output of lm function.
  l <- list( a = format(coef(m)[1], digits = 3, nsmall=3),   b = format(abs(coef(m)[2]), digits = 3, nsmall=4),   r2 = format(summary(m)$r.squared, digits = 3, nsmall=3), pv = format(as.numeric(unlist(summary(m)$coefficients[,4][2] )), digits = 3, nsmall=5) )                                                          
  if (coef(m)[2] >= 0)  {
    eq <- substitute(y == a~+~b* "x," ~~r^2~"="~r2* ","  ~~p~"="~pv, l)  ## substitute函数中需要替换的变量用列表参数方式给出。~ is space.
  } else {
    eq <- substitute(y == a~-~b* "x," ~~r^2~"="~r2* ","  ~~p~"="~pv, l)    
  }
  as.character(as.expression(eq));                 
}


my_scatter_equation <- function( xAxis2,  yAxis2,   path2,  fileName2, titel2, xlab2, ylab2,   yMin2=0, yMax2=1,   xMin2=0, xMax2=1 ) { 
  #xAxis2[xAxis2<xMin2] = xMin2
  #xAxis2[xAxis2>xMax2] = xMax2
  #yAxis2[yAxis2<yMin2] = yMin2
  #yAxis2[yAxis2>yMax2] = yMax2
  dataFrame_NTR_point <- data.frame( xAxis=xAxis2,  yAxis=yAxis2 ) 
  regressionA <- lm(yAxis2 ~ xAxis2)
  FigureTemp1 <- ggplot(data = dataFrame_NTR_point, aes(x = xAxis, y = yAxis)) + 
    geom_point(size=1, alpha=1, color="red3") + geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +  
    geom_text( aes(x = xMax2-0.5, y = yMax2-0.05,  label = MyLmEquation_1( regressionA )), parse = TRUE,  colour="cyan4", size=2.5 ) +
    xlab(xlab2) +  ylab(ylab2) + xlim(xMin2, xMax2) +  scale_y_continuous(breaks=c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5 ), limits=c(yMin2, yMax2) ) +
    ggtitle(titel2) + MyTheme_1_g( hjust1=1, vjust1=1, textSize=12 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-scatter-equation", sep="",  collapse=NULL),  height1=4, width1=4)
}



sink( paste(outDir_g, "Levels.txt", sep="/") )
print( rawLevel )
print( allMutations)
sink()


sink( paste(outDir_g, "Regression.log.txt", sep="/") ) 
MyRegression_1( xAxis2 = rawLevel,  yAxis2 = allMutations,  file2 = paste(outDir_g, "Regression", sep="/") )
sink()

my_scatter_equation( xAxis2 = rawLevel,  yAxis2 = allMutations, 
                     path2=outDir_g,  fileName2=motif_g, titel2="", xlab2="m6A fraction", ylab2="Mutation frequency" ,
                     yMin2=0.00, yMax2=1.00,   xMin2=0, xMax2=1 )



 






