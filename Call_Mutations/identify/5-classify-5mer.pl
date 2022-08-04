#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "4-fasta/all.fasta", global variable.
my $output_g = '';  ## such as "5-classified", global variable.
{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Classify A sites into 256 groups based on their 5-mer contexts.

        Usage:
               perl  5-classify-5mer.pl    [-version]    [-help]     [-in inputFile]    [-out outDir]
        For instance:
               perl  5-classify-5mer.pl    -in 4-fasta/all.fasta   -out 5-classified    > 5-classify-5mer.runLog.txt

        Optional arguments:
        -version        Show version number of this program and exit.
        -help           Show this help message and exit.

        Required arguments:
        -in inputFile        "inputFile" is the name of your input file.  (no default)
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
$input_g  = '4-fasta/all.fasta';             ## This is only an initialization value or suggesting value, not default value.
$output_g = '5-classified';                  ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  5-classify-5mer.pl  -help' \n";
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
my $output2_g = "$output_g/each_5mer";
&myMakeDir($output2_g);
open(my $FH_input_g, $input_g)  ||  die;
my @lines = <$FH_input_g> ; 

open(fileA1, ">", "$output_g/DRACH.txt")  or  die;
open(fileA2, ">", "$output_g/DRACH.fasta")  or  die;
open(fileB1, ">", "$output_g/nonDRACH.txt")  or  die;
open(fileB2, ">", "$output_g/nonDRACH.fasta")  or  die;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
## >chr1-14362-T-G---234-1-0.43%-T---164-1-0.61%-T---Reference-1.0-0.6554511278195432---0-164-0-1---0-234-0-1......0.18......0.6554511278195432......165::chr1:14359-14364(-)  
for(my $i=0; $i<$#lines; $i=$i+2) {
        my $head = $lines[$i];
        my $sequ = $lines[$i+1];
        $head =~ m/^>/ or die;
        $head =~ s/\n$// or die;
        $sequ =~ tr/agct/AGCT/ ;
        $sequ =~ m/^[AGCTN][AGCTN]A[AGCTN][AGCTN]\n$/ or die "##$sequ##";
       
        if( $sequ =~ m/^[AGT][AG]AC[ACT]\n$/ ) {
            $sequ =~ m/^([AGCT][AGCT]A[AGCT][AGCT])\n$/ or die;
            my $seq2 = $1;
            $seq2 =~ tr/agct/AGCT/ ;
            print fileA2 "$head\n$seq2\n";
            print fileA1 "$head......$seq2\n";
        }else{
            $sequ =~ m/^([AGCTN][AGCTN]A[AGCTN][AGCTN])\n$/ or die "##$sequ##";
            my $seq2 = $1;
            $seq2 =~ tr/agctn/AGCTN/ ;
            print fileB2 "$head\n$seq2\n";
            print fileB1 "$head......$seq2\n";
        } 
        $sequ =~ s/\n$// or die;
        open(fileC, ">>", "$output2_g/$sequ.txt")  or  die "\n\n\n$output2_g/$sequ.txt\n\n\n";
        print fileC "$head......$sequ\n";
}





#####
