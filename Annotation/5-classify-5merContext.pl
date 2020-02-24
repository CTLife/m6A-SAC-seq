#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "all.fasta", global variable.
my $output_g = '';  ## such as "All-5mer-Contexts", global variable.
{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Classify A sites into 256 groups based on their 5-mer contexts.

        Usage:
               perl  4-classify-5merContext.pl    [-version]    [-help]     [-in inputDir]    [-out outDir]
        For instance:
               perl  4-classify-5merContext.pl    -in all.fasta   -out All-5mer-Contexts    > 4-classify-5merContext.runLog.txt

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
$input_g  = 'all.fasta';             ## This is only an initialization value or suggesting value, not default value.
$output_g = 'All-5mer-Contexts';     ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  4-classify-5merContext.pl  -help' \n";
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

open(my $FH_input_g, $input_g)  ||  die;
my @lines = <$FH_input_g> ; 


open(fileA1, ">", "$output_g/DRACH.txt")  or  die;
open(fileA2, ">", "$output_g/DRACH.fa")  or  die;
open(fileB1, ">", "$output_g/nonDRACH.txt")  or  die;
open(fileB2, ">", "$output_g/nonDRACH.fa")  or  die;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
for(my $i=0; $i<$#lines; $i=$i+2) {
        my $head = $lines[$i];
        my $sequ = $lines[$i+1];
        $head =~ m/^>/ or die;
        $sequ =~ m/^[AGCTN][AGCTN]A[AGCTN][AGCTN]\n$/i or die "##$sequ##";

        if( $sequ =~ m/^[AGT][AG]AC[ACT]\n$/i ) {
            $sequ =~ m/^([AGCT][AGCT]A[AGCT][AGCT])\n$/i or die;
            my $seq2 = $1;
            $seq2 =~ tr/agct/AGCT/ ;
            print fileA2 "$head$sequ";
            print fileA1 "$head";
        }else{
            $sequ =~ m/^([AGCTN][AGCTN]A[AGCTN][AGCTN])\n$/i or die "##$sequ##";
            my $seq2 = $1;
            $seq2 =~ tr/agctn/AGCTN/ ;
            print fileB2 "$head$sequ";
            print fileB1 "$head";
        }
}





#####
