#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "1-raw/lu1_r1.minus.varscan.validation", global variable.
my $output_g = '';  ## such as "2-sepAGCT/lu1_r1.minus.varscan.validation", global variable.
{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Separate each line based on the 3rd column "ref" in the varscan2 output files.

        Usage:
               perl  sepAGCT_varscan.pl    [-version]    [-help]     [-in inputDir]    [-out outDir]
        For instance:
               perl  sepAGCT_varscan.pl    -in 1-raw/lu1_r1.minus.varscan.validation   -out 2-sepAGCT/lu1_r1.minus.varscan.validation    > sepAGCT_varscan.runLog.txt

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
$input_g  = '1-raw/lu1_r1.minus.varscan.validation';         ## This is only an initialization value or suggesting value, not default value.
$output_g = '2-sepAGCT/lu1_r1.minus.varscan.validation';     ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  sepAGCT_varscan.pl  -help' \n";
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

&myMakeDir("$output_g");
open(my $FH_input_g, $input_g)  ||  die;
my $strand_g = "+";
if($input_g =~ m/\.minus\./) {$strand_g = "-"; }

open(fileA, ">", "$output_g/A.txt")  or  die;
open(fileG, ">", "$output_g/G.txt")  or  die;
open(fileC, ">", "$output_g/C.txt")  or  die;
open(fileT, ">", "$output_g/T.txt")  or  die;
open(fileN, ">", "$output_g/N.txt")  or  die;
open(file_non, ">", "$output_g/inconsistency.txt")  or  die;

open(normal_bg_fileA, ">", "$output_g/normal_A.bg")  or  die;
open(normal_bg_fileG, ">", "$output_g/normal_G.bg")  or  die;
open(normal_bg_fileC, ">", "$output_g/normal_C.bg")  or  die;
open(normal_bg_fileT, ">", "$output_g/normal_T.bg")  or  die;
open(normal_bg_fileN, ">", "$output_g/normal_N.bg")  or  die;

open(normal_bed_fileA, ">", "$output_g/normal_A.bed")  or  die;
open(normal_bed_fileG, ">", "$output_g/normal_G.bed")  or  die;
open(normal_bed_fileC, ">", "$output_g/normal_C.bed")  or  die;
open(normal_bed_fileT, ">", "$output_g/normal_T.bed")  or  die;
open(normal_bed_fileN, ">", "$output_g/normal_N.bed")  or  die;

open(tumor_bg_fileA, ">", "$output_g/tumor_A.bg")  or  die;
open(tumor_bg_fileG, ">", "$output_g/tumor_G.bg")  or  die;
open(tumor_bg_fileC, ">", "$output_g/tumor_C.bg")  or  die;
open(tumor_bg_fileT, ">", "$output_g/tumor_T.bg")  or  die;
open(tumor_bg_fileN, ">", "$output_g/tumor_N.bg")  or  die;

open(tumor_bed_fileA, ">", "$output_g/tumor_A.bed")  or  die;
open(tumor_bed_fileG, ">", "$output_g/tumor_G.bed")  or  die;
open(tumor_bed_fileC, ">", "$output_g/tumor_C.bed")  or  die;
open(tumor_bed_fileT, ">", "$output_g/tumor_T.bed")  or  die;
open(tumor_bed_fileN, ">", "$output_g/tumor_N.bed")  or  die;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $line1 = <$FH_input_g>;
print fileA $line1;
print fileG $line1;
print fileC $line1;
print fileT $line1;
print fileN $line1;
print file_non $line1;

while(my $line=<$FH_input_g>) {
        #$line =~ m/^chr/ or die "##$line##\n" ;
        $line =~ m/^(\S+)\t(\d+)\t([A-Z])\t(.*)\t(\d+)\t(\d+)\t(\S+)\t(\S*)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t/ or die "##$line##";
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

        if($var eq "") {$var = "*";}
        my $name  = "$chrom-$position-$ref-$var--$normal_reads1-$normal_reads2-$normal_var_freq-$normal_gt--$tumor_reads1-$tumor_reads2-$tumor_var_freq-$tumor_gt";
        my $start = $position - 1;
        my $end   = $position;
        $normal_var_freq =~ s/\%$// or die;
        $normal_var_freq = $normal_var_freq/100;
        $tumor_var_freq  =~ s/\%$// or die;
        $tumor_var_freq  = $tumor_var_freq/100;
        
        given($ref) {
        	when("A") {print fileA $line; 
                           print normal_bg_fileA  "$chrom\t$start\t$end\t$normal_var_freq\n";  
                           print normal_bed_fileA "$chrom\t$start\t$end\t$name\t$normal_var_freq\t$strand_g\n"; 
                           print tumor_bg_fileA   "$chrom\t$start\t$end\t$tumor_var_freq\n";  
                           print tumor_bed_fileA  "$chrom\t$start\t$end\t$name\t$tumor_var_freq\t$strand_g\n"; 
                          }
		when("G") {print fileG $line; 
                           print normal_bg_fileG  "$chrom\t$start\t$end\t$normal_var_freq\n";  
                           print normal_bed_fileG "$chrom\t$start\t$end\t$name\t$normal_var_freq\t$strand_g\n"; 
                           print tumor_bg_fileG   "$chrom\t$start\t$end\t$tumor_var_freq\n";  
                           print tumor_bed_fileG  "$chrom\t$start\t$end\t$name\t$tumor_var_freq\t$strand_g\n"; 
                          }
		when("C") {print fileC $line; 
                           print normal_bg_fileC  "$chrom\t$start\t$end\t$normal_var_freq\n";  
                           print normal_bed_fileC "$chrom\t$start\t$end\t$name\t$normal_var_freq\t$strand_g\n"; 
                           print tumor_bg_fileC   "$chrom\t$start\t$end\t$tumor_var_freq\n";  
                           print tumor_bed_fileC  "$chrom\t$start\t$end\t$name\t$tumor_var_freq\t$strand_g\n"; 
                          }
		when("T") {print fileT $line; 
                           print normal_bg_fileT  "$chrom\t$start\t$end\t$normal_var_freq\n";  
                           print normal_bed_fileT "$chrom\t$start\t$end\t$name\t$normal_var_freq\t$strand_g\n"; 
                           print tumor_bg_fileT   "$chrom\t$start\t$end\t$tumor_var_freq\n";  
                           print tumor_bed_fileT  "$chrom\t$start\t$end\t$name\t$tumor_var_freq\t$strand_g\n"; 
                          }
		default   {print fileN $line; 
                           print normal_bg_fileN  "$chrom\t$start\t$end\t$normal_var_freq\n";  
                           print normal_bed_fileN "$chrom\t$start\t$end\t$name\t$normal_var_freq\t$strand_g\n"; 
                           print tumor_bg_fileN   "$chrom\t$start\t$end\t$tumor_var_freq\n";  
                           print tumor_bed_fileN  "$chrom\t$start\t$end\t$name\t$tumor_var_freq\t$strand_g\n"; 
                          }
        }
        if( ($ref ne $normal_gt) or ($ref ne $tumor_gt) ){print file_non $line; }
}


my $chromSize_g = "/home/yp/.MyProgramFiles/RefGenome/hg38/Standard_genome_sequence_files/hg38.chrom.sizes";


system("cat  $output_g/normal_A.bg   |  sort  -k1,1  -k2,2n     >  $output_g/normal_A.sorted.bg");
system("bedGraphToBigWig  $output_g/normal_A.sorted.bg   $chromSize_g      $output_g/normal_A.bw ");

system("cat  $output_g/normal_G.bg   |  sort  -k1,1  -k2,2n     >  $output_g/normal_G.sorted.bg");
system("bedGraphToBigWig  $output_g/normal_G.sorted.bg   $chromSize_g      $output_g/normal_G.bw ");

system("cat  $output_g/normal_C.bg   |  sort  -k1,1  -k2,2n     >  $output_g/normal_C.sorted.bg");
system("bedGraphToBigWig  $output_g/normal_C.sorted.bg   $chromSize_g      $output_g/normal_C.bw ");

system("cat  $output_g/normal_T.bg   |  sort  -k1,1  -k2,2n     >  $output_g/normal_T.sorted.bg");
system("bedGraphToBigWig  $output_g/normal_T.sorted.bg   $chromSize_g      $output_g/normal_T.bw ");

system("cat  $output_g/normal_N.bg   |  sort  -k1,1  -k2,2n     >  $output_g/normal_N.sorted.bg");
system("bedGraphToBigWig  $output_g/normal_N.sorted.bg   $chromSize_g      $output_g/normal_N.bw ");


system("cat  $output_g/tumor_A.bg   |  sort  -k1,1  -k2,2n     >  $output_g/tumor_A.sorted.bg");
system("bedGraphToBigWig  $output_g/tumor_A.sorted.bg   $chromSize_g      $output_g/tumor_A.bw ");

system("cat  $output_g/tumor_G.bg   |  sort  -k1,1  -k2,2n     >  $output_g/tumor_G.sorted.bg");
system("bedGraphToBigWig  $output_g/tumor_G.sorted.bg   $chromSize_g      $output_g/tumor_G.bw ");

system("cat  $output_g/tumor_C.bg   |  sort  -k1,1  -k2,2n     >  $output_g/tumor_C.sorted.bg");
system("bedGraphToBigWig  $output_g/tumor_C.sorted.bg   $chromSize_g      $output_g/tumor_C.bw ");

system("cat  $output_g/tumor_T.bg   |  sort  -k1,1  -k2,2n     >  $output_g/tumor_T.sorted.bg");
system("bedGraphToBigWig  $output_g/tumor_T.sorted.bg   $chromSize_g      $output_g/tumor_T.bw ");

system("cat  $output_g/tumor_N.bg   |  sort  -k1,1  -k2,2n     >  $output_g/tumor_N.sorted.bg");
system("bedGraphToBigWig  $output_g/tumor_N.sorted.bg   $chromSize_g      $output_g/tumor_N.bw ");





#####
