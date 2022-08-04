#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "1-finalSites/HEK293.bed", global variable.
my $output_g = '';  ## such as "2-getFraction/HEK293.bed",    global variable.
my $file_g = '';    ## such as "A_C_merge.spikeins.fit.txt",    global variable.

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Usage:
               perl  2-getFraction.pl    [-version]    [-help]  [-f file]   [-in inputFile]    [-out outFile]
        For instance:
               perl  2-getFraction.pl    -f A_C_merge.spikeins.fit.txt   -in 1-finalSites/HEK293.bed   -out 2-getFraction/HEK293.bed    

        Optional arguments:
        -version        Show version number of this program and exit.
        -help           Show this help message and exit.

        Required arguments:
        -f file   Linear regression results.
        -in inputFile        "inputFile" is the name of your input file.  (no default)
        -out outFile         "outFile" is the name of output file.  (no default)
        
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    Version 0.2,  2021-04-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '1-finalSites/HEK293.bed';      ## This is only an initialization value or suggesting value, not default value.
$output_g = '2-getFraction/HEK293.bed';     ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  -f ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  2-getFraction.pl  -help' \n";
    exit 0;
}

## Get Arguments
if ( exists $args{'-version' }   )     { say  "\n$version\n";    exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP\n";       exit 0; }
if ( exists $args{'-f'      }   )     { $file_g   = $args{'-f'      };  }else{say   "\n -f     is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in'      }; }else{say   "\n -in     is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'     }; }else{say   "\n -out    is required.\n";   say  "\n$HELP\n";    exit 0; }

## Conditions
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$output_g =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                File: $file_g
                Input   File:  $input_g
                Output  File:  $output_g
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
system( "rm -rf  $output_g" );

my $file_a = $file_g;
open(fileIn, "<",  $input_g )        or  die;
open(fileAa, "<",  $file_a  )        or  die;
open(fileO1, ">",  $output_g)        or  die;

my @lines1 = <fileIn> ; 
my @lines2 = <fileAa> ; 
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my @motifs        = @lines2;
my @coefficients1 = @lines2;
my @coefficients2 = @lines2;
for(my $i=0; $i<=$#lines2; $i=$i+1) {
        my $temp = $lines2[$i];
        $temp =~ m/^(\S+)\s+([e\-\.\d]+)\s+([e\-\.\d]+)\n/  or  die "##$temp##"; 
        $motifs[$i] = $1;	
        $coefficients1[$i] = $2;		
        $coefficients2[$i] = $3;		        
}


for(my $i=0; $i<=$#lines1; $i=$i+1) {
        my $temp = $lines1[$i];
        $temp =~ m/^(chr\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+([\.\d]+)\s+([-+])\s/  or  die "##$temp##"; 
        my $chrom = $1;	
        my $start = $2;		
        my $end   = $3;		
        my $name  = $4;
        my $mutation = $5;
        my $strand   = $6;
        $mutation = $mutation/100; 
        next unless $chrom =~ m/^chr[\dXYM]+$/;

        $name =~ m/\.{6}([AGCT]{5})/  or  die "##$name\n##";
        my $fiveMer = $1;

        my $bool = 0;
        my $fraction = -1; 

        for(my $j=0; $j<=$#lines2; $j=$j+1) {
            if($motifs[$j] eq $fiveMer) {$fraction = ($mutation-$coefficients2[$j])/$coefficients1[$j];  $bool=$bool+1; }			
        }

        ($bool == 1) or die "##$bool##$fiveMer##\n";
        if($fraction < 0.01) {$fraction = 0.01}; 
        $fraction = $fraction*100; 
        $mutation = $mutation*100; 
        if($fraction > 100) {$fraction = 100; } 

        print  fileO1   "$chrom\t$start\t$end\t$name\t$mutation\t$strand\t$fraction\n";
}

###################################################################################################################################################################################################






say   "Done ......";






#####
