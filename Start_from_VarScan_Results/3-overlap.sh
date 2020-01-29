minus_FTO1_vs_Input_A="2-furtherRemove/FTO1_vs_Input/minus/select2.tumor_A.bed"
minus_FTO1_vs_Input_C="2-furtherRemove/FTO1_vs_Input/minus/select2.tumor_C.bed"
minus_FTO1_vs_Input_G="2-furtherRemove/FTO1_vs_Input/minus/select2.tumor_G.bed"
minus_FTO1_vs_Input_T="2-furtherRemove/FTO1_vs_Input/minus/select2.tumor_T.bed"
plus_FTO1_vs_Input_A="2-furtherRemove/FTO1_vs_Input/plus/select2.tumor_A.bed"
plus_FTO1_vs_Input_C="2-furtherRemove/FTO1_vs_Input/plus/select2.tumor_C.bed"
plus_FTO1_vs_Input_G="2-furtherRemove/FTO1_vs_Input/plus/select2.tumor_G.bed"
plus_FTO1_vs_Input_T="2-furtherRemove/FTO1_vs_Input/plus/select2.tumor_T.bed"

minus_FTO2_vs_Input_A="2-furtherRemove/FTO2_vs_Input/minus/select2.tumor_A.bed"
minus_FTO2_vs_Input_C="2-furtherRemove/FTO2_vs_Input/minus/select2.tumor_C.bed"
minus_FTO2_vs_Input_G="2-furtherRemove/FTO2_vs_Input/minus/select2.tumor_G.bed"
minus_FTO2_vs_Input_T="2-furtherRemove/FTO2_vs_Input/minus/select2.tumor_T.bed"
plus_FTO2_vs_Input_A="2-furtherRemove/FTO2_vs_Input/plus/select2.tumor_A.bed"
plus_FTO2_vs_Input_C="2-furtherRemove/FTO2_vs_Input/plus/select2.tumor_C.bed"
plus_FTO2_vs_Input_G="2-furtherRemove/FTO2_vs_Input/plus/select2.tumor_G.bed"
plus_FTO2_vs_Input_T="2-furtherRemove/FTO2_vs_Input/plus/select2.tumor_T.bed"



minus_WT1_vs_FTO_A="2-furtherRemove/WT1_vs_FTO/minus/select2.tumor_A.bed"
minus_WT1_vs_FTO_C="2-furtherRemove/WT1_vs_FTO/minus/select2.tumor_C.bed"
minus_WT1_vs_FTO_G="2-furtherRemove/WT1_vs_FTO/minus/select2.tumor_G.bed"
minus_WT1_vs_FTO_T="2-furtherRemove/WT1_vs_FTO/minus/select2.tumor_T.bed"
plus_WT1_vs_FTO_A="2-furtherRemove/WT1_vs_FTO/plus/select2.tumor_A.bed"
plus_WT1_vs_FTO_C="2-furtherRemove/WT1_vs_FTO/plus/select2.tumor_C.bed"
plus_WT1_vs_FTO_G="2-furtherRemove/WT1_vs_FTO/plus/select2.tumor_G.bed"
plus_WT1_vs_FTO_T="2-furtherRemove/WT1_vs_FTO/plus/select2.tumor_T.bed"

minus_WT2_vs_FTO_A="2-furtherRemove/WT2_vs_FTO/minus/select2.tumor_A.bed"
minus_WT2_vs_FTO_C="2-furtherRemove/WT2_vs_FTO/minus/select2.tumor_C.bed"
minus_WT2_vs_FTO_G="2-furtherRemove/WT2_vs_FTO/minus/select2.tumor_G.bed"
minus_WT2_vs_FTO_T="2-furtherRemove/WT2_vs_FTO/minus/select2.tumor_T.bed"
plus_WT2_vs_FTO_A="2-furtherRemove/WT2_vs_FTO/plus/select2.tumor_A.bed"
plus_WT2_vs_FTO_C="2-furtherRemove/WT2_vs_FTO/plus/select2.tumor_C.bed"
plus_WT2_vs_FTO_G="2-furtherRemove/WT2_vs_FTO/plus/select2.tumor_G.bed"
plus_WT2_vs_FTO_T="2-furtherRemove/WT2_vs_FTO/plus/select2.tumor_T.bed"



minus_WT1_vs_Input_A="2-furtherRemove/WT1_vs_Input/minus/select2.tumor_A.bed"
minus_WT1_vs_Input_C="2-furtherRemove/WT1_vs_Input/minus/select2.tumor_C.bed"
minus_WT1_vs_Input_G="2-furtherRemove/WT1_vs_Input/minus/select2.tumor_G.bed"
minus_WT1_vs_Input_T="2-furtherRemove/WT1_vs_Input/minus/select2.tumor_T.bed"
plus_WT1_vs_Input_A="2-furtherRemove/WT1_vs_Input/plus/select2.tumor_A.bed"
plus_WT1_vs_Input_C="2-furtherRemove/WT1_vs_Input/plus/select2.tumor_C.bed"
plus_WT1_vs_Input_G="2-furtherRemove/WT1_vs_Input/plus/select2.tumor_G.bed"
plus_WT1_vs_Input_T="2-furtherRemove/WT1_vs_Input/plus/select2.tumor_T.bed"

minus_WT2_vs_Input_A="2-furtherRemove/WT2_vs_Input/minus/select2.tumor_A.bed"
minus_WT2_vs_Input_C="2-furtherRemove/WT2_vs_Input/minus/select2.tumor_C.bed"
minus_WT2_vs_Input_G="2-furtherRemove/WT2_vs_Input/minus/select2.tumor_G.bed"
minus_WT2_vs_Input_T="2-furtherRemove/WT2_vs_Input/minus/select2.tumor_T.bed"
plus_WT2_vs_Input_A="2-furtherRemove/WT2_vs_Input/plus/select2.tumor_A.bed"
plus_WT2_vs_Input_C="2-furtherRemove/WT2_vs_Input/plus/select2.tumor_C.bed"
plus_WT2_vs_Input_G="2-furtherRemove/WT2_vs_Input/plus/select2.tumor_G.bed"
plus_WT2_vs_Input_T="2-furtherRemove/WT2_vs_Input/plus/select2.tumor_T.bed"





outPath1="3-overlap/FTO_vs_Input/minus"
outPath2="3-overlap/FTO_vs_Input/plus"
mkdir -p $outPath1/A
mkdir -p $outPath1/C
mkdir -p $outPath1/G
mkdir -p $outPath1/T
mkdir -p $outPath2/A
mkdir -p $outPath2/C
mkdir -p $outPath2/G
mkdir -p $outPath2/T

intervene  venn  --input   $minus_FTO1_vs_Input_A     $minus_FTO2_vs_Input_A       --names=A,B             \
--output $outPath1/A          --save-overlaps      --dpi 1200    --figtype svg     

intervene  venn  --input   $minus_FTO1_vs_Input_C     $minus_FTO2_vs_Input_C       --names=A,B             \
--output $outPath1/C          --save-overlaps      --dpi 1200    --figtype svg     
 

intervene  venn  --input   $minus_FTO1_vs_Input_G     $minus_FTO2_vs_Input_G       --names=A,B             \
--output $outPath1/G          --save-overlaps      --dpi 1200    --figtype svg     
 

intervene  venn  --input   $minus_FTO1_vs_Input_T     $minus_FTO2_vs_Input_T       --names=A,B             \
--output $outPath1/T          --save-overlaps      --dpi 1200    --figtype svg     
 


intervene  venn  --input   $plus_FTO1_vs_Input_A     $plus_FTO2_vs_Input_A       --names=A,B             \
--output $outPath2/A          --save-overlaps      --dpi 1200    --figtype svg     

intervene  venn  --input   $plus_FTO1_vs_Input_C     $plus_FTO2_vs_Input_C       --names=A,B             \
--output $outPath2/C          --save-overlaps      --dpi 1200    --figtype svg     
 

intervene  venn  --input   $plus_FTO1_vs_Input_G     $plus_FTO2_vs_Input_G       --names=A,B             \
--output $outPath2/G          --save-overlaps      --dpi 1200    --figtype svg     
 

intervene  venn  --input   $plus_FTO1_vs_Input_T     $plus_FTO2_vs_Input_T       --names=A,B             \
--output $outPath2/T          --save-overlaps      --dpi 1200    --figtype svg     
 




outPath1="3-overlap/WT_vs_Input/minus"
outPath2="3-overlap/WT_vs_Input/plus"
mkdir -p $outPath1/A
mkdir -p $outPath1/C
mkdir -p $outPath1/G
mkdir -p $outPath1/T
mkdir -p $outPath2/A
mkdir -p $outPath2/C
mkdir -p $outPath2/G
mkdir -p $outPath2/T

intervene  venn  --input   $minus_WT1_vs_Input_A     $minus_WT2_vs_Input_A       --names=A,B             \
--output $outPath1/A          --save-overlaps      --dpi 1200    --figtype svg     

intervene  venn  --input   $minus_WT1_vs_Input_C     $minus_WT2_vs_Input_C       --names=A,B             \
--output $outPath1/C          --save-overlaps      --dpi 1200    --figtype svg     
 

intervene  venn  --input   $minus_WT1_vs_Input_G     $minus_WT2_vs_Input_G       --names=A,B             \
--output $outPath1/G          --save-overlaps      --dpi 1200    --figtype svg     
 

intervene  venn  --input   $minus_WT1_vs_Input_T     $minus_WT2_vs_Input_T       --names=A,B             \
--output $outPath1/T          --save-overlaps      --dpi 1200    --figtype svg     
 


intervene  venn  --input   $plus_WT1_vs_Input_A     $plus_WT2_vs_Input_A       --names=A,B             \
--output $outPath2/A          --save-overlaps      --dpi 1200    --figtype svg     

intervene  venn  --input   $plus_WT1_vs_Input_C     $plus_WT2_vs_Input_C       --names=A,B             \
--output $outPath2/C          --save-overlaps      --dpi 1200    --figtype svg     
 

intervene  venn  --input   $plus_WT1_vs_Input_G     $plus_WT2_vs_Input_G       --names=A,B             \
--output $outPath2/G          --save-overlaps      --dpi 1200    --figtype svg     
 

intervene  venn  --input   $plus_WT1_vs_Input_T     $plus_WT2_vs_Input_T       --names=A,B             \
--output $outPath2/T          --save-overlaps      --dpi 1200    --figtype svg     
 







outPath1="3-overlap/WT_vs_FTO/minus"
outPath2="3-overlap/WT_vs_FTO/plus"
mkdir -p $outPath1/A
mkdir -p $outPath1/C
mkdir -p $outPath1/G
mkdir -p $outPath1/T
mkdir -p $outPath2/A
mkdir -p $outPath2/C
mkdir -p $outPath2/G
mkdir -p $outPath2/T

intervene  venn  --input   $minus_WT1_vs_FTO_A     $minus_WT2_vs_FTO_A       --names=A,B             \
--output $outPath1/A          --save-overlaps      --dpi 1200    --figtype svg     

intervene  venn  --input   $minus_WT1_vs_FTO_C     $minus_WT2_vs_FTO_C       --names=A,B             \
--output $outPath1/C          --save-overlaps      --dpi 1200    --figtype svg     
 

intervene  venn  --input   $minus_WT1_vs_FTO_G     $minus_WT2_vs_FTO_G       --names=A,B             \
--output $outPath1/G          --save-overlaps      --dpi 1200    --figtype svg     
 

intervene  venn  --input   $minus_WT1_vs_FTO_T     $minus_WT2_vs_FTO_T       --names=A,B             \
--output $outPath1/T          --save-overlaps      --dpi 1200    --figtype svg     
 


intervene  venn  --input   $plus_WT1_vs_FTO_A     $plus_WT2_vs_FTO_A       --names=A,B             \
--output $outPath2/A          --save-overlaps      --dpi 1200    --figtype svg     

intervene  venn  --input   $plus_WT1_vs_FTO_C     $plus_WT2_vs_FTO_C       --names=A,B             \
--output $outPath2/C          --save-overlaps      --dpi 1200    --figtype svg     
 

intervene  venn  --input   $plus_WT1_vs_FTO_G     $plus_WT2_vs_FTO_G       --names=A,B             \
--output $outPath2/G          --save-overlaps      --dpi 1200    --figtype svg     
 

intervene  venn  --input   $plus_WT1_vs_FTO_T     $plus_WT2_vs_FTO_T       --names=A,B             \
--output $outPath2/T          --save-overlaps      --dpi 1200    --figtype svg     
 
    
 


