#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "sets/11_A_B.bed", global variable.
my $output_g = '';  ## such as "11_A_B.bed", global variable.
{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Separate each line based on the 3rd column "ref" in the varscan2 output files.

        Usage:
               perl  select.pl    [-version]    [-help]     [-in inputDir]    [-out outDir]
        For instance:
               perl  select.pl    -in sets/11_A_B.bed   -out 11_A_B.bed    > select.runLog.txt

        Optional arguments:
        -version        Show version number of this program and exit.
        -help           Show this help message and exit.

        Required arguments:
        -in inputDir        "inputDir" is the name of your input file.  (no default)
        -out outDir         "outDir" is the name of output path that contains your running results of this step.  (no default)
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
$input_g  = 'sets/11_A_B.bed';         ## This is only an initialization value or suggesting value, not default value.
$output_g = '11_A_B.bed';         ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  select.pl  -help' \n";
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

 
open(my $FH_input_g, $input_g)  ||  die;
open(fileA, ">", "select1.$output_g")  or  die;
open(fileB, ">", "select2.$output_g")  or  die;

open(fileA1, ">", "Others_select1.$output_g")  or  die;
open(fileB1, ">", "Others_select2.$output_g")  or  die;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
 
while(my $line=<$FH_input_g>) {
        $line =~ m/^(\S+)\t(\d+)\t(\d+)\t(\S+)\t/ or die "##$line##";
        my $name  = $4;
        $name =~ m/\-([\d.]+)%\-\S+\-([\d.]+)%\-/ or die;
        my $normal_var_freq = $1;
        my $tumor_var_freq  = $2;
        if( ($normal_var_freq < $tumor_var_freq)  )  {
             print  fileA   $line ; 
             if($normal_var_freq < 5)  {print  fileB   $line ;}else{print  fileB1   $line ;}
        }else{
             print  fileA1   $line ; 
        }
        
}

 



#####
