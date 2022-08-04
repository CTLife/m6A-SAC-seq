prefix1=input        ## normal 
prefix2=scA0-5op_S15    ## turmor
pileup1="2-Pileup/plus_strand/$prefix1.pileup"    ## normal 
pileup2="2-Pileup/plus_strand/$prefix2.pileup"    ## turmor
out="3-VarScan/$prefix1.vs.$prefix2/plus"
mkdir -p   $out
commands="java -jar  /home/yp/.MyProgramFiles/7_Downstream/VarScan/VarScan.jar"
##Tumor-normal Comparison:
$commands  somatic     $pileup1  $pileup2  $out    --validation 1   > $out/runLog.txt 2>&1
  

prefix1=input        ## normal 
prefix2=scA5op_S14    ## turmor
pileup1="2-Pileup/plus_strand/$prefix1.pileup"    ## normal 
pileup2="2-Pileup/plus_strand/$prefix2.pileup"    ## turmor
out="3-VarScan/$prefix1.vs.$prefix2/plus"
mkdir -p   $out
commands="java -jar  /home/yp/.MyProgramFiles/7_Downstream/VarScan/VarScan.jar"
##Tumor-normal Comparison:
$commands  somatic     $pileup1  $pileup2  $out    --validation 1   > $out/runLog.txt 2>&1
  

prefix1=input        ## normal 
prefix2=scA50op_S13    ## turmor
pileup1="2-Pileup/plus_strand/$prefix1.pileup"    ## normal 
pileup2="2-Pileup/plus_strand/$prefix2.pileup"    ## turmor
out="3-VarScan/$prefix1.vs.$prefix2/plus"
mkdir -p   $out
commands="java -jar  /home/yp/.MyProgramFiles/7_Downstream/VarScan/VarScan.jar"
##Tumor-normal Comparison:
$commands  somatic     $pileup1  $pileup2  $out    --validation 1   > $out/runLog.txt 2>&1
  




  


prefix1=input        ## normal 
prefix2=scA0-5op_S15    ## turmor
pileup1="2-Pileup/minus_strand/$prefix1.pileup"    ## normal 
pileup2="2-Pileup/minus_strand/$prefix2.pileup"    ## turmor
out="3-VarScan/$prefix1.vs.$prefix2/minus"
mkdir -p   $out
commands="java -jar  /home/yp/.MyProgramFiles/7_Downstream/VarScan/VarScan.jar"
##Tumor-normal Comparison:
$commands  somatic     $pileup1  $pileup2  $out    --validation 1   > $out/runLog.txt 2>&1
  

prefix1=input        ## normal 
prefix2=scA5op_S14    ## turmor
pileup1="2-Pileup/minus_strand/$prefix1.pileup"    ## normal 
pileup2="2-Pileup/minus_strand/$prefix2.pileup"    ## turmor
out="3-VarScan/$prefix1.vs.$prefix2/minus"
mkdir -p   $out
commands="java -jar  /home/yp/.MyProgramFiles/7_Downstream/VarScan/VarScan.jar"
##Tumor-normal Comparison:
$commands  somatic     $pileup1  $pileup2  $out    --validation 1   > $out/runLog.txt 2>&1
  

prefix1=input        ## normal 
prefix2=scA50op_S13    ## turmor
pileup1="2-Pileup/minus_strand/$prefix1.pileup"    ## normal 
pileup2="2-Pileup/minus_strand/$prefix2.pileup"    ## turmor
out="3-VarScan/$prefix1.vs.$prefix2/minus"
mkdir -p   $out
commands="java -jar  /home/yp/.MyProgramFiles/7_Downstream/VarScan/VarScan.jar"
##Tumor-normal Comparison:
$commands  somatic     $pileup1  $pileup2  $out    --validation 1   > $out/runLog.txt 2>&1
  




  

