#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "2-sepAGCT/minusStrand/T.txt", global variable.
my $output_g = '';  ## such as "3-BEDregion/minusStrand/T", global variable.
{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Separate each line based on the 3rd column "ref" in the varscan2 output files.

        Usage:
               perl  3-BEDregion.pl    [-version]    [-help]     [-in inputDir]    [-out outDir]
        For instance:
               perl  3-BEDregion.pl    -in 2-sepAGCT/minusStrand/T.txt   -out 3-BEDregion/minusStrand/T    > 3-BEDregion.runLog.txt

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
$input_g  = '2-sepAGCT/minusStrand/T.txt';         ## This is only an initialization value or suggesting value, not default value.
$output_g = '3-BEDregion/minusStrand/T';            ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  3-BEDregion.pl  -help' \n";
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
my $output2_g = $output_g;
$output2_g =~ s/\/[AGCT]$// or die;
&myMakeDir($output2_g);

my $strand_g = "wrong";
my $bool_temp = 0;
if($input_g =~ m/minus/i) {$strand_g = "-"; $bool_temp++; }
if($input_g =~ m/plus/i)  {$strand_g = "+"; $bool_temp++; }
($bool_temp == 1)  or die;

open(FHin,  "<",   $input_g)   or  die;
open(FHout1, ">",  "$output_g.bed")  or  die;
open(refX_tumorA, ">", "$output_g.to.A.bed")  or  die;
open(refX_tumorC, ">", "$output_g.to.C.bed")  or  die;
open(refX_tumorG, ">", "$output_g.to.G.bed")  or  die;
open(refX_tumorT, ">", "$output_g.to.T.bed")  or  die;
open(refX_tumorN, ">", "$output_g.to.N.bed")  or  die;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my @lines = <FHin>;

for(my $i=1; $i<=$#lines; $i++) {
    $lines[$i] =~ m/^(\S+)\t(\d+)\t([A-Z]+)\t(\S*)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\n$/ or die;

    my $chrom = $1;	
    my $position = $2;	
    my $ref = $3;	
    my $var = $4;
	
    my $normal_reads1 = $5;	
    my $normal_reads2 = $6;	
    my $normal_var_freq = $7;	
    my $normal_gt = $8;
	
    my $tumor_reads1 = $9;	
    my $tumor_reads2 = $10;	
    my $tumor_var_freq = $11;	
    my $tumor_gt = $12;	

    my $somatic_status  = $13;	
    my $variant_p_value = $14;	
    my $somatic_p_value = $15;	

    my $tumor_reads1_plus  = $16;	
    my $tumor_reads1_minus = $17;	
    my $tumor_reads2_plus  = $18;	
    my $tumor_reads2_minus = $19;	

    my $normal_reads1_plus  = $20;	
    my $normal_reads1_minus = $21;	
    my $normal_reads2_plus  = $22;	
    my $normal_reads2_minus = $23;

    ##if($var eq "") {$var = "*";}
    my $start  = $position - 3; 
    my $end    = $position + 2; 
    my $name1  = "$chrom-$position-$ref-$var---$normal_reads1-$normal_reads2-$normal_var_freq-$normal_gt---$tumor_reads1-$tumor_reads2-$tumor_var_freq-$tumor_gt";
    my $name2  = "$somatic_status-$variant_p_value-$somatic_p_value---$tumor_reads1_plus-$tumor_reads1_minus-$tumor_reads2_plus-$tumor_reads2_minus---$normal_reads1_plus-$normal_reads1_minus-$normal_reads2_plus-$normal_reads2_minus";                               
    my $name   = "$name1---$name2"; 
    $normal_var_freq =~ s/%// or die;
    $tumor_var_freq =~ s/%// or die;
    my $score  = $tumor_var_freq - $normal_var_freq;
    
    my $pvalue = $somatic_p_value;
    if($pvalue > $variant_p_value) {$pvalue = $variant_p_value; }
    my $tumor_reads = $tumor_reads1 + $tumor_reads2 ; 
    $name   = "$name......$score......$pvalue......$tumor_reads"; 

    if(  $tumor_var_freq > $normal_var_freq  ) {  ## So, unmutated sites are remvoed.
            ($score > 0)  or die "$tumor_var_freq - $normal_var_freq";        
            print FHout1 "$chrom\t$start\t$end\t$name\t$score\t$strand_g\n";
            if($var eq "C") {
                print   refX_tumorC    "$chrom\t$start\t$end\t$name\t$score\t$strand_g\n"; 
            }elsif($var eq "G") {
                print   refX_tumorG    "$chrom\t$start\t$end\t$name\t$score\t$strand_g\n"; 
            }elsif($var eq "T") {
                print   refX_tumorT    "$chrom\t$start\t$end\t$name\t$score\t$strand_g\n"; 
            }elsif($var eq "A") {
                print   refX_tumorA    "$chrom\t$start\t$end\t$name\t$score\t$strand_g\n"; 
            }else{
                print   refX_tumorN    "$chrom\t$start\t$end\t$name\t$score\t$strand_g\n";  
            }

    }

}




#####
print  "Done!!!\n\n\n\n\n";















#####
