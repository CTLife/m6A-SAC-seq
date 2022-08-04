#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';      ## such as "1-raw", global variable.
my $background_g = '';  ## such as 0.05 (5%)
my $P_g    = '';
my $diff_g = '';
my $cov_g  = '';

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------

        Usage:
               perl  identify.pl    [-version]    [-help]     [-in inputDir]    [-background 0.05]    [-p 0.1]    [-diff 0.1]    [-cov 10]
        For instance:
               perl  identify.pl    -in 1-raw    -background 0.05    -p 0.1    -diff 0.1    -cov 10    > identify.runLog.txt

        Optional arguments:
        -version        Show version number of this program and exit.
        -help           Show this help message and exit.

        Required arguments:
        -in inputDir            Path of input files.  (no default)
        -background mutation    Maximum background mutation rate.  (no default)
        -p pvalue               P-value.  (no default)
        -diff difference        Mimimum mutation rate (background was removed).  (no default)
        -cov coverage           Reads coverage.  (no default)
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
$input_g  = '1-raw';         
$background_g = 0.05;  ## such as 0.05 (5%)
$P_g    = 0.1;
$diff_g = 0.1;
$cov_g  = 10;   

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  identify.pl  -help' \n";
    exit 0;
}

## Get Arguments
if ( exists $args{'-version'   }   )     { say  "\n$version\n";    exit 0; }
if ( exists $args{'-help'      }   )     { say  "\n$HELP\n";       exit 0; }
if ( exists $args{'-in'        }   )     { $input_g      = $args{'-in'          }; }else{say   "\n -in     is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-background'}   )     { $background_g = $args{'-background'  }; }else{say   "\n -background    is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-p'         }   )     { $P_g          = $args{'-p'           }; }else{say   "\n -p       is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-diff'      }   )     { $diff_g       = $args{'-diff'        }; }else{say   "\n -diff    is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-cov'       }   )     { $cov_g        = $args{'-cov'         }; }else{say   "\n -cov     is required.\n";   say  "\n$HELP\n";    exit 0; }

## Conditions
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$background_g =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$P_g    =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$diff_g =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$cov_g  =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Input   File:  $input_g
                Background  :  $background_g
                P-value     :  $P_g
                Mimimum mutation rate  :  $diff_g
                Coverage    :  $cov_g
        ##########################################################
\n";
}
###################################################################################################################################################################################################



system( " perl  2-sepAGCT_varscan.pl    -in $input_g/minus.snp     -out 2-sepAGCT/minusStrand ");
system( " perl  2-sepAGCT_varscan.pl    -in $input_g/plus.snp      -out 2-sepAGCT/plusStrand ");

system( " perl  3-BEDregion.pl    -in 2-sepAGCT/minusStrand/A.txt    -out 3-BEDregion/minusStrand/A ");
system( " perl  3-BEDregion.pl    -in 2-sepAGCT/minusStrand/T.txt    -out 3-BEDregion/minusStrand/T ");

system( " perl  3-BEDregion.pl    -in 2-sepAGCT/plusStrand/A.txt    -out 3-BEDregion/plusStrand/A ");
system( " perl  3-BEDregion.pl    -in 2-sepAGCT/plusStrand/T.txt    -out 3-BEDregion/plusStrand/T ");

system( " bash  4-ToFasta.sh ");
 
system( " perl  5-classify-5mer.pl    -in 4-fasta/all.fasta   -out 5-classified ");

system( " perl  6-BED.pl    -in 5-classified/DRACH.txt      -out 6-BED/DRACH ");
system( " perl  6-BED.pl    -in 5-classified/nonDRACH.txt   -out 6-BED/nonDRACH ");






