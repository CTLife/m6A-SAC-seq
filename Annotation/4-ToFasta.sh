out="4-fasta"
mkdir -p $out

bedtools getfasta  -fo $out/T.fasta    -name  -s    -fi hg38.fa   -fullHeader   -bed 3-BEDregion/minusStrand/T.bed                
bedtools getfasta  -fo $out/A.fasta    -name  -s    -fi hg38.fa   -fullHeader   -bed 3-BEDregion/plusStrand/A.bed                


