#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "2-twoToOne/DRACH.txt", global variable.
my $output_g = '';  ## such as "3-full-infor/DRACH.txt", global variable.
{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Classify A sites into 256 groups based on their 5-mer contexts.

        Usage:
               perl  3-full-infor.pl    [-version]    [-help]     [-in inputDir]    [-out outDir]
        For instance:
               perl  3-full-infor.pl    -in 2-twoToOne/DRACH.txt   -out 3-full-infor/DRACH.txt    > 3-full-infor.runLog.txt

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
$input_g  = '2-twoToOne/DRACH.txt';       ## This is only an initialization value or suggesting value, not default value.
$output_g = '3-full-infor/DRACH.txt';     ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  3-full-infor.pl  -help' \n";
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

## &myMakeDir($output_g);

open(my $FH_input_g, $input_g)  ||  die;
my @lines = <$FH_input_g> ; 
open(fileA, ">", $output_g)  or  die;
###################################################################################################################################################################################################




print fileA "5mer_chrom\t5mer_start\t5mer_end\tname\tvar_freq\tstrand\tA_position\t5mer\n";
###################################################################################################################################################################################################
for(my $i=0; $i<=$#lines; $i=$i+1) {
        my $temp = $lines[$i];
        $temp =~ m/^>chr/ or die;
        $temp =~ m/^>(chr\S+)\-(\d+)\-(\S)\-(\S+)\-\-(\d+)\-(\d+)\-(\S+)\-(\S)\-\-(\d+)\-(\d+)\-(\S+)\-(\S)::(\S+):(\d+)\-(\d+)\(([-+])\)([ACGTN]{5})\n$/i  or  die  "$temp";
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

        my $chr2 = $13;
        my $start = $14;
        my $end = $15;
        my $strand = $16;
        my $mer5 = $17;
        my $name = "$chrom-$position-$ref-$var-$normal_reads1-$normal_reads2-$normal_var_freq-$normal_gt-$tumor_reads1-$tumor_reads2-$tumor_var_freq-$tumor_gt";
        $tumor_var_freq =~ s/\%$// or die;
        $tumor_var_freq = $tumor_var_freq/100;
        ($tumor_var_freq > 0)  or die;
        $mer5 =~ tr/agct/AGCT/ ;
        ($chr2 eq $chrom)  or  die;
        ($position == $start+3)  or  die;
        ($position == $end-2)  or  die;

        print fileA "$chr2\t$start\t$end\t$name\t$tumor_var_freq\t$strand\t$position\t$mer5\n";

}





#####
