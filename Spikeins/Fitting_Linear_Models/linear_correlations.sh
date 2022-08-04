perl  2-probesFrequency.pl    -in 1-seq   -out 2-probesFrequency 
bash  3-fit_motif.sh
bash  3-fit_motif_5levels.sh
bash  4-DRACH.select.sh
bash  4-DRACH.select-5levels.sh
Rscript  5-histogram.R
perl  6-format.pl

