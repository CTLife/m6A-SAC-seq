A1_1="2-getFraction/sc0.5ng.bed";

A1_3="2-getFraction/sc5ng.bed";

A1_5="2-getFraction/sc50ng.bed";

 
A1_1a="2-getFraction/HeLa.m6A.Peaks.bed";
A1_2a="2-getFraction/Lulu_Hela.bed";
A1_3a="2-getFraction/Lulu_Hela_RIP.bed";
A1_4a="2-getFraction/A1_11_common.bed";
A1_5a="2-getFraction/A2_11_common.bed";
A1_6a="2-getFraction/As1_11_common.bed";
A1_7a="2-getFraction/As2_11_common.bed";
A1_8a="2-getFraction/Ass1_11_common.bed";
A1_9a="2-getFraction/Ass2_11_common.bed";
  

##################################
 


outpath1="6-overlap/pairwise.fraction"
mkdir -p $outpath1

intervene  pairwise  \
--input  \
$A1_1	\
$A1_3	\
$A1_5	\
$A1_1a	\
$A1_2a	\
$A1_3a	\
$A1_4a	\
$A1_5a	\
$A1_6a	\
$A1_7a	\
$A1_8a	\
$A1_9a	\
--compute  frac   --htype color   --output $outpath1          \
    --dpi 1200    --figtype svg     --figsize 10 10   --fontsize 10   









outpath2="6-overlap/pairwise.count"
mkdir -p $outpath2

intervene  pairwise  \
--input  \
$A1_1	\
$A1_3	\
$A1_5	\
$A1_1a	\
$A1_2a	\
$A1_3a	\
$A1_4a	\
$A1_5a	\
$A1_6a	\
$A1_7a	\
$A1_8a	\
$A1_9a	\
--compute  count   --htype color   --output $outpath2          \
       --dpi 1200    --figtype svg     --figsize 10 10   --fontsize 10   







