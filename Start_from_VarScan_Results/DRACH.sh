## conserved motif DRACH (D=G/A/U, R=G/A, H=A/U/C) 

cat FP.fasta | grep -iP "^[AGT][AG]AC[ACT]" | wc -l  >> number.DRACH.txt
cat TP.fasta | grep -iP "^[AGT][AG]AC[ACT]" | wc -l  >> number.DRACH.txt


