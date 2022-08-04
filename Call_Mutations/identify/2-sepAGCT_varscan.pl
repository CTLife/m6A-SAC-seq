#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "1-raw/minus.T.txt", global variable.
my $output_g = '';  ## such as "2-sepAGCT/minus.T.txt", global variable.
{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Separate each line based on the 3rd column "ref" in the varscan2 output files.

        Usage:
               perl  2-sepAGCT_varscan.pl    [-version]    [-help]     [-in inputFile]    [-out outFile]
        For instance:
               perl  2-sepAGCT_varscan.pl    -in 1-raw/minus.T.txt   -out 2-sepAGCT/minus.T.txt    > 2-sepAGCT_varscan.runLog.txt

        Optional arguments:
        -version        Show version number of this program and exit.
        -help           Show this help message and exit.

        Required arguments:
        -in inputFile        "inputFile" is the name of your input file which is the output of VarScan.  (no default)
        -out outFile         "outFile" is the name of the running results of this script.  (no default)
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
$input_g  = '1-raw/minus.T.txt';         ## This is only an initialization value or suggesting value, not default value.
$output_g = '2-sepAGCT/minus.T.txt';     ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  2-sepAGCT_varscan.pl  -help' \n";
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

my $strand_g = "wrong";
my $bool_temp = 0;
if($input_g =~ m/minus\./) {$strand_g = "-"; $bool_temp++; }
if($input_g =~ m/plus\./)  {$strand_g = "+"; $bool_temp++; }
($bool_temp == 1)  or die;

open(fileA, ">", "$output_g/A.txt")  or  die;
open(fileG, ">", "$output_g/G.txt")  or  die;
open(fileC, ">", "$output_g/C.txt")  or  die;
open(fileT, ">", "$output_g/T.txt")  or  die;
open(fileN, ">", "$output_g/N.txt")  or  die;
open(nonMatch, ">", "$output_g/nonMatch.txt")  or  die;

open(normal_bed_fileA, ">", "$output_g/Znormal_A.bed")  or  die;
open(normal_bed_fileG, ">", "$output_g/Znormal_G.bed")  or  die;
open(normal_bed_fileC, ">", "$output_g/Znormal_C.bed")  or  die;
open(normal_bed_fileT, ">", "$output_g/Znormal_T.bed")  or  die;
open(normal_bed_fileN, ">", "$output_g/Znormal_N.bed")  or  die;

open(tumor_bed_fileA, ">", "$output_g/Ztumor_A.bed")  or  die;
open(tumor_bed_fileG, ">", "$output_g/Ztumor_G.bed")  or  die;
open(tumor_bed_fileC, ">", "$output_g/Ztumor_C.bed")  or  die;
open(tumor_bed_fileT, ">", "$output_g/Ztumor_T.bed")  or  die;
open(tumor_bed_fileN, ">", "$output_g/Ztumor_N.bed")  or  die;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
open(my $FH_input_g, $input_g)  ||  die;
my $line1 = <$FH_input_g>;
print fileA $line1;
print fileG $line1;
print fileC $line1;
print fileT $line1;
print fileN $line1;

my $n  = 0;
my $n1 = 0;
my $n2 = 0;
my $n3 = 0;
my $n4 = 0;
my $n5 = 0;

while(my $line=<$FH_input_g>) {
    if($line =~ m/^(\S+)\t(\d+)\t([A-Z])\t(\S*)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\n$/) {
    $line =~ m/^(\S+)\t(\d+)\t([A-Z])\t(\S*)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\n$/ or die "##$line##\n";
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

    if($strand_g eq "+") {
        $tumor_reads1_minus==0  or  die; $tumor_reads2_minus==0  or  die; 
        $normal_reads1_minus==0 or  die; $normal_reads2_minus==0 or  die; 
    }

    if($strand_g eq "-") {
        $tumor_reads1_plus==0  or  die; $tumor_reads2_plus==0  or  die; 
        $normal_reads1_plus==0 or  die; $normal_reads2_plus==0 or  die; 
    }

    if($var eq "") {$var = "*";}
    my $start  = $position - 1; 
    my $end    = $position; 
    my $name1  = "$chrom-$position-$ref-$var-$normal_reads1-$normal_reads2-$normal_var_freq-$normal_gt-$tumor_reads1-$tumor_reads2-$tumor_var_freq-$tumor_gt";
    my $name2  = "$somatic_status-$variant_p_value-$somatic_p_value-$tumor_reads1_plus-$tumor_reads1_minus-$tumor_reads2_plus-$tumor_reads2_minus-$normal_reads1_plus-$normal_reads1_minus-$normal_reads2_plus-$normal_reads2_minus";                               
    my $name   = "$name1-$name2"; 
    $normal_var_freq =~ s/%// or die;
    $tumor_var_freq =~ s/%// or die;

    $n++;
    given($ref) {
        	when("A") {  print fileA $line; 
                            print normal_bed_fileA "$chrom\t$start\t$end\t$name\t$normal_var_freq\t$strand_g\n"; 
                            print tumor_bed_fileA  "$chrom\t$start\t$end\t$name\t$tumor_var_freq\t$strand_g\n"; 
                            $n1++;
                          }
		when("G") { print fileG $line; 
                           print normal_bed_fileG "$chrom\t$start\t$end\t$name\t$normal_var_freq\t$strand_g\n"; 
                           print tumor_bed_fileG  "$chrom\t$start\t$end\t$name\t$tumor_var_freq\t$strand_g\n"; 
                           $n2++;
                          }
		when("C") { print fileC $line; 
                           print normal_bed_fileC "$chrom\t$start\t$end\t$name\t$normal_var_freq\t$strand_g\n"; 
                           print tumor_bed_fileC  "$chrom\t$start\t$end\t$name\t$tumor_var_freq\t$strand_g\n"; 
                           $n3++;
                          }
		when("T") { print fileT $line; 
                           print normal_bed_fileT "$chrom\t$start\t$end\t$name\t$normal_var_freq\t$strand_g\n"; 
                           print tumor_bed_fileT  "$chrom\t$start\t$end\t$name\t$tumor_var_freq\t$strand_g\n"; 
                           $n4++;
                          }
		default   { print fileN $line; 
                           print normal_bed_fileN "$chrom\t$start\t$end\t$name\t$normal_var_freq\t$strand_g\n"; 
                           print tumor_bed_fileN  "$chrom\t$start\t$end\t$name\t$tumor_var_freq\t$strand_g\n"; 
                           $n5++;
                          }
        }
    }else{
        print nonMatch   $line; 
    }
}



my $m = $n1+$n2+$n3+$n4+$n5;
if($n == $m){
         print "OK:\n";
         print "$n == $m; $n1,$n2,$n3,$n4,$n5\n\n";
}else{
         print "Wrong:\n";
         print "$n != $m; $n1,$n2,$n3,$n4,$n5\n\n";     
}


#####
print  "Done!!!\n\n\n\n\n";













