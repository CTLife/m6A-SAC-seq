minus_FTO_vs_Input_A="FTO_vs_Input/minus/A/sets/11_A_B.bed"
minus_FTO_vs_Input_C="FTO_vs_Input/minus/C/sets/11_A_B.bed"
minus_FTO_vs_Input_G="FTO_vs_Input/minus/G/sets/11_A_B.bed"
minus_FTO_vs_Input_T="FTO_vs_Input/minus/T/sets/11_A_B.bed"
plus_FTO_vs_Input_A="FTO_vs_Input/plus/A/sets/11_A_B.bed"
plus_FTO_vs_Input_C="FTO_vs_Input/plus/C/sets/11_A_B.bed"
plus_FTO_vs_Input_G="FTO_vs_Input/plus/G/sets/11_A_B.bed"
plus_FTO_vs_Input_T="FTO_vs_Input/plus/T/sets/11_A_B.bed"

minus_WT_vs_FTO_A="WT_vs_FTO/minus/A/sets/11_A_B.bed"
minus_WT_vs_FTO_C="WT_vs_FTO/minus/C/sets/11_A_B.bed"
minus_WT_vs_FTO_G="WT_vs_FTO/minus/G/sets/11_A_B.bed"
minus_WT_vs_FTO_T="WT_vs_FTO/minus/T/sets/11_A_B.bed"
plus_WT_vs_FTO_A="WT_vs_FTO/plus/A/sets/11_A_B.bed"
plus_WT_vs_FTO_C="WT_vs_FTO/plus/C/sets/11_A_B.bed"
plus_WT_vs_FTO_G="WT_vs_FTO/plus/G/sets/11_A_B.bed"
plus_WT_vs_FTO_T="WT_vs_FTO/plus/T/sets/11_A_B.bed"

minus_WT_vs_Input_A="WT_vs_Input/minus/A/sets/11_A_B.bed"
minus_WT_vs_Input_C="WT_vs_Input/minus/C/sets/11_A_B.bed"
minus_WT_vs_Input_G="WT_vs_Input/minus/G/sets/11_A_B.bed"
minus_WT_vs_Input_T="WT_vs_Input/minus/T/sets/11_A_B.bed"
plus_WT_vs_Input_A="WT_vs_Input/plus/A/sets/11_A_B.bed"
plus_WT_vs_Input_C="WT_vs_Input/plus/C/sets/11_A_B.bed"
plus_WT_vs_Input_G="WT_vs_Input/plus/G/sets/11_A_B.bed"
plus_WT_vs_Input_T="WT_vs_Input/plus/T/sets/11_A_B.bed"




outPath1="3groups_overlap/minus"
mkdir -p $outPath1/A
mkdir -p $outPath1/C
mkdir -p $outPath1/G
mkdir -p $outPath1/T

intervene  venn  --input   $minus_WT_vs_Input_A     $minus_WT_vs_FTO_A    $minus_FTO_vs_Input_A    --names=A,B,C             \
--output $outPath1/A          --save-overlaps      --dpi 1200    --figtype svg     

intervene  venn  --input   $minus_WT_vs_Input_C     $minus_WT_vs_FTO_C    $minus_FTO_vs_Input_C    --names=A,B,C             \
--output $outPath1/C          --save-overlaps      --dpi 1200    --figtype svg     

intervene  venn  --input   $minus_WT_vs_Input_G     $minus_WT_vs_FTO_G    $minus_FTO_vs_Input_G    --names=A,B,C             \
--output $outPath1/G          --save-overlaps      --dpi 1200    --figtype svg     

intervene  venn  --input   $minus_WT_vs_Input_T     $minus_WT_vs_FTO_T    $minus_FTO_vs_Input_T    --names=A,B,C             \
--output $outPath1/T          --save-overlaps      --dpi 1200    --figtype svg     


outPath1="3groups_overlap/plus"
mkdir -p $outPath1/A
mkdir -p $outPath1/C
mkdir -p $outPath1/G
mkdir -p $outPath1/T

intervene  venn  --input   $plus_WT_vs_Input_A     $plus_WT_vs_FTO_A    $plus_FTO_vs_Input_A    --names=A,B,C             \
--output $outPath1/A          --save-overlaps      --dpi 1200    --figtype svg     

intervene  venn  --input   $plus_WT_vs_Input_C     $plus_WT_vs_FTO_C    $plus_FTO_vs_Input_C    --names=A,B,C             \
--output $outPath1/C          --save-overlaps      --dpi 1200    --figtype svg     

intervene  venn  --input   $plus_WT_vs_Input_G     $plus_WT_vs_FTO_G    $plus_FTO_vs_Input_G    --names=A,B,C             \
--output $outPath1/G          --save-overlaps      --dpi 1200    --figtype svg     

intervene  venn  --input   $plus_WT_vs_Input_T     $plus_WT_vs_FTO_T    $plus_FTO_vs_Input_T    --names=A,B,C             \
--output $outPath1/T          --save-overlaps      --dpi 1200    --figtype svg     







