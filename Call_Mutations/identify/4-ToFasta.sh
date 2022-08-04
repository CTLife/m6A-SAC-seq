out="4-fasta"
mkdir -p $out

bedtools getfasta  -fo $out/T.fasta    -name  -s    -fi hg38.genome.fa   -fullHeader   -bed 3-BEDregion/minusStrand/T.bed                
bedtools getfasta  -fo $out/A.fasta    -name  -s    -fi hg38.genome.fa   -fullHeader   -bed 3-BEDregion/plusStrand/A.bed                

cat   $out/T.fasta  $out/A.fasta   >  $out/all.fasta

