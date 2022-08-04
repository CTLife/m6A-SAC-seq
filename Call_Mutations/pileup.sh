inputDir=1-Split-strand/minus_strand
outDir=2-Pileup/minus_strand

mkdir -p $outDir


 
bam1=input
samtools mpileup  --fasta-ref hg38.genome.fa   --output $outDir/$bam1.pileup   $inputDir/$bam1.bam   
 
bam1=scA0-5op_S15
samtools mpileup  --fasta-ref hg38.genome.fa   --output $outDir/$bam1.pileup   $inputDir/$bam1.bam   
  
bam1=scA5op_S14
samtools mpileup  --fasta-ref hg38.genome.fa   --output $outDir/$bam1.pileup   $inputDir/$bam1.bam   
  
bam1=scA50op_S13
samtools mpileup  --fasta-ref hg38.genome.fa   --output $outDir/$bam1.pileup   $inputDir/$bam1.bam   
  
 
 
 
 

 inputDir=1-Split-strand/plus_strand
outDir=2-Pileup/plus_strand

mkdir -p $outDir


 
bam1=input
samtools mpileup  --fasta-ref hg38.genome.fa   --output $outDir/$bam1.pileup   $inputDir/$bam1.bam   
 
bam1=scA0-5op_S15
samtools mpileup  --fasta-ref hg38.genome.fa   --output $outDir/$bam1.pileup   $inputDir/$bam1.bam   
  
bam1=scA5op_S14
samtools mpileup  --fasta-ref hg38.genome.fa   --output $outDir/$bam1.pileup   $inputDir/$bam1.bam   
  
bam1=scA50op_S13
samtools mpileup  --fasta-ref hg38.genome.fa   --output $outDir/$bam1.pileup   $inputDir/$bam1.bam   
  
 
 
 
 

 
