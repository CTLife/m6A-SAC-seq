#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "5-5mers/DRACH.txt", global variable.
my $output_g = '';  ## such as "6-BED", global variable.
{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Classify A sites into 256 groups based on their 5-mer contexts.

        Usage:
               perl  6-BED.pl    [-version]    [-help]     [-in inputDir]    [-out outDir]
        For instance:
               perl  6-BED.pl    -in 5-5mers/DRACH.txt   -out 6-BED    > 6-BED.runLog.txt

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
my $version = "    Version 0.1,  2020-02-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '5-5mers/DRACH.txt';             ## This is only an initialization value or suggesting value, not default value.
$output_g = '6-BED';     ## This is only an initialization value or suggesting value, not default value.

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
open(my $FH_input1_g, $input_g)  ||  die;
my $prefix_g = $input_g;

$prefix_g =~ s/^\S+\/// or die;

my @lines1 = <$FH_input1_g> ; 
open(fileA, ">", "$output_g/$prefix_g" )  or  die;
open(fileA1, ">", "$output_g/$prefix_g.number.txt" )  or  die;
open(fileB, ">", "$output_g/$prefix_g.10percent.bed" )  or  die;
open(fileC, ">", "$output_g/$prefix_g.5percent.bed" )  or  die;
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
## >chr1-14644-T-G-62-0-0%-T-27-1-3.57%-T-Reference-1.0-0.31111111111111067-0-27-0-1-0-62-0-0::chr1:14641-14646(-)
##say "11111";
for(my $i=0; $i<=$#lines1; $i=$i+1) {
        my $temp = $lines1[$i];
        $temp =~ m/^>(chr[^-]+)\-(\d+)\-\S\-\S+\-\d+\-\d+\-([.%\d]+)\-\S{1,30}\-\d+\-\d+\-([.%\d]+)\-\S{1,30}\-[a-z]+\-.+\(([-+])\)\n$/i  or  die "##$temp##"; 
        my $chrom = $1;	
        my $position = $2;		
        my $normal_var_freq = $3;		
        my $tumor_var_freq = $4;
        my $strand_g = $5;

	my $start  = $position - 1; 
        my $end    = $position ; 
        $normal_var_freq =~ s/%// or die;
        $tumor_var_freq  =~ s/%// or die;
        my $score  = $tumor_var_freq - $normal_var_freq;  
        my $name = $temp;
        $name =~ s/\n$// or die;
        $name =~ s/^>// or die;

        print   fileA    "$chrom\t$start\t$end\t$name\t$score\t$strand_g\n";

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

        if( $score >= 10 ) {print   fileB    "$chrom\t$start\t$end\t$name\t$score\t$strand_g\n";}
        if( $score >= 5  ) {print   fileC    "$chrom\t$start\t$end\t$name\t$score\t$strand_g\n";}

}
##say "22222";
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

print   fileA1  "(0, 5):,   $n1, $ratio_n1\n";
print   fileA1  "[5,  10):, $n2, $ratio_n2\n";
print   fileA1  "[10, 20):, $n3, $ratio_n3\n";
print   fileA1  "[20, 30):, $n4, $ratio_n4\n";
print   fileA1  "[30, 40):, $n5, $ratio_n5\n";
print   fileA1  "[40, 50):, $n6, $ratio_n6\n";
print   fileA1  "[50, 60):, $n7, $ratio_n7\n";
print   fileA1  "[60, 70):, $n8, $ratio_n8\n";
print   fileA1  "[70, 80):, $n9, $ratio_n9\n";
print   fileA1  "[80, 90):, $n10, $ratio_n10\n";
print   fileA1  "[90, 100):, $n11, $ratio_n11\n";
###################################################################################################################################################################################################






say   "Done ......";






#####
