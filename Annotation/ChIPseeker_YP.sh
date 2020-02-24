sleep  5h

Rscript ChIPseeker_YP.R    \
        --inputDir=1-BED    \
        --outDir=2-Annotation   \
        --refGenome=hg38  \
        --stranded=TRUE    \
        --fast=TRUE        \
        --upstream=1000      \
        --downstream=0  

