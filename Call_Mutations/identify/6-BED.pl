#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "5-5mers/DRACH.txt", global variable.
my $output_g = '';  ## such as "6-BED/DRACH", global variable.
{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Classify A sites into 256 groups based on their 5-mer contexts.

        Usage:
               perl  6-BED.pl    [-version]    [-help]     [-in inputDir]    [-out outDir]
        For instance:
               perl  6-BED.pl    -in 5-5mers/DRACH.txt   -out 6-BED/DRACH    > 6-BED.runLog.txt

        Optional arguments:
        -version        Show version number of this program and exit.
        -help           Show this help message and exit.

        Required arguments:
        -in inputDir        "inputDir" is the name of your input file.  (no default)
        -out outDir         "outDir" is the name of output path that contains your running results of this script.  (no default)
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    Version 0.2,  2020-08-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '5-5mers/DRACH.txt';             ## This is only an initialization value or suggesting value, not default value.
$output_g = '6-BED/DRACH';                   ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  6-BED.pl  -help' \n";
    exit 0;
}

## Get Arguments
if ( exists $args{'-version' }   )     { say  "\n$version\n";    exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP\n";       exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in'      }; }else{say   "\n -in     is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'     }; }else{say   "\n -out    is required.\n";   say  "\n$HELP\n";    exit 0; }

## Conditions
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$output_g =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Input   File:  $input_g
                Output  Path:  $output_g
        ##########################################################
\n";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Running ......";

sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die;       }
}
&myMakeDir($output_g);

open(my $FH1_g, $input_g)  ||  die;
my @lines1 = <$FH1_g> ; 

open(fileA,  ">", "$output_g/1mer.bed" )  or  die;
open(fileB,  ">", "$output_g/5mer.bed" )  or  die;
open(fileC,  ">", "$output_g/numb.txt" )  or  die;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $n1 = 0;
my $n2 = 0;
my $n3 = 0;
my $n4 = 0;
my $n5 = 0;
my $n6 = 0;
my $n7 = 0;
my $n8 = 0;
my $n9 = 0;
my $n10 = 0;
my $n11 = 0;
## >chr1-14362-T-G---234-1-0.43%-T---164-1-0.61%-T---Reference-1.0-0.6554511278195432---0-164-0-1---0-234-0-1......0.18......0.6554511278195432......165::chr1:14359-14364(-)......GGATG
## >chr1-14517-T-A---36-0-0%-T---32-1-3.03%-T---Reference-1.0-0.47826086956525227---0-32-0-1---0-36-0-0......3.03......0.47826086956525227......33::chr1:14514-14519(-)......GGACT
for(my $i=0; $i<=$#lines1; $i=$i+1) {
        my $temp = $lines1[$i];
        $temp =~ m/^>(chr\w+)\-(\d+)\-\w\-\S+\-\-\-\d+\-\d+\-(\d+\.?\d*\%)\-\S{1,50}\-\-\-\d+\-\d+\-(\d+\.?\d*\%)\-\S{1,50}\-\-\-[a-z]+\-\d\S+\.{6}([\d\.]+)\.{6}([-\d\.E]+)\.{6}(\d+)\:\:(chr\w+)\:(\d+)\-(\d+)\(([-+])\)\.{6}([AGCTN]{2}A[AGCTN]{2})/i  or  die "##$temp##";                                
        my $chrom    = $1;	
        my $position = $2;		
        my $score    = $5;
        my $pvalue   = $6;
        my $support  = $7;
        my $chrom2   = $8;
        my $start2   = $9;
        my $end2     = $10;
        my $strand   = $11;
        
        ($chrom eq $chrom2) or die;
        ($position-3 == $start2)  or  die;
        ($position+2 == $end2)    or  die;
                  
        my $name = $temp;
        $name =~ s/\n$// or die;
        $name =~ s/^>// or die;
        my $start3 = $position-1;
        my $end3 = $position;
        
        if( ($score <= 2)  or ($pvalue >= 0.1) )  {next;}
        if( $support <= 5 )  {next;}
        
        print   fileA    "$chrom\t$start3\t$end3\t$name\t$score\t$strand\n";
        print   fileB    "$chrom\t$start2\t$end2\t$name\t$score\t$strand\n";
        
        my $percentage = $score ;
        my $bool = 0;
        if(  $percentage==0  ) { $bool++; } 
        if( ($percentage>0)   and ($percentage<5)  ) {$n1++; $bool++; }
        if( ($percentage>=5)  and ($percentage<10) ) {$n2++; $bool++; }
        if( ($percentage>=10) and ($percentage<20) ) {$n3++; $bool++; }
        if( ($percentage>=20) and ($percentage<30) ) {$n4++; $bool++; }
        if( ($percentage>=30) and ($percentage<40) ) {$n5++; $bool++; }
        if( ($percentage>=40) and ($percentage<50) ) {$n6++; $bool++; }
        if( ($percentage>=50) and ($percentage<60) ) {$n7++; $bool++; }
        if( ($percentage>=60) and ($percentage<70) ) {$n8++; $bool++; }
        if( ($percentage>=70) and ($percentage<80) ) {$n9++; $bool++; }
        if( ($percentage>=80) and ($percentage<90) ) {$n10++; $bool++; }
        if( ($percentage>=90) and ($percentage<101)) {$n11++; $bool++; }
        ($bool == 1) or die;
}
 
my $total_n = $#lines1 + 1;

my $ratio_n1 = 100*$n1/$total_n;
my $ratio_n2 = 100*$n2/$total_n;
my $ratio_n3 = 100*$n3/$total_n;
my $ratio_n4 = 100*$n4/$total_n;
my $ratio_n5 = 100*$n5/$total_n;
my $ratio_n6 = 100*$n6/$total_n;
my $ratio_n7 = 100*$n7/$total_n;
my $ratio_n8 = 100*$n8/$total_n;
my $ratio_n9 = 100*$n9/$total_n;
my $ratio_n10 = 100*$n10/$total_n;
my $ratio_n11 = 100*$n11/$total_n;

print   fileC  "(0, 5):,   $n1, $ratio_n1\n";
print   fileC  "[5,  10):, $n2, $ratio_n2\n";
print   fileC  "[10, 20):, $n3, $ratio_n3\n";
print   fileC  "[20, 30):, $n4, $ratio_n4\n";
print   fileC  "[30, 40):, $n5, $ratio_n5\n";
print   fileC  "[40, 50):, $n6, $ratio_n6\n";
print   fileC  "[50, 60):, $n7, $ratio_n7\n";
print   fileC  "[60, 70):, $n8, $ratio_n8\n";
print   fileC  "[70, 80):, $n9, $ratio_n9\n";
print   fileC  "[80, 90):, $n10, $ratio_n10\n";
print   fileC  "[90, 100):, $n11, $ratio_n11\n";
###################################################################################################################################################################################################






say   "Done ......";






#####
