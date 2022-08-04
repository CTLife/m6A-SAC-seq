#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "1_rawFASTQ", global variable.
my $output_g = '';  ## such as "Spikeins", global variable.
{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Usage:
               perl  Spikeins.pl    [-version]    [-help]     [-in inputDir]    [-out outDir]
        For instance:
               nohup   time  perl  Spikeins.pl    -in 1_rawFASTQ   -out Spikeins    > Spikeins.runLog.txt   2>&1  &   

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
my $version = "    Version 1.2,  2022-08-02.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '1_rawFASTQ';       ## This is only an initialization value or suggesting value, not default value.
$output_g = 'Spikeins';         ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  Spikeins.pl  -help' \n";
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
opendir(my $FH_input1_g, $input_g)  ||  die;
my @files = readdir($FH_input1_g);
###################################################################################################################################################################################################



# TATCTGTCTCGACGTNNANNGGCCTTTGCAACTAGAATTACACCATAATTGCT
# TATCTGTCTCGACGTNNANNGGCATTCAAGCCTAGAATTACACCATAATTGCT
# TATCTGTCTCGACGTNNANNGGCGAGGTGATCTAGAATTACACCATAATTGCT
# TATCTGTCTCGACGTNNANNGGCTTCAACAACTAGAATTACACCATAATTGCT
# TATCTGTCTCGACGTNNANNGGCGATGGTTTCTAGAATTACACCATAATTGCT

# TATCTGTCTCGACGT[AGCT]{5}GGCCTTTGCAACTAGAATTACACCATAATTGCT 
# TATCTGTCTCGACGT[AGCT]{5}GGCATTCAAGCCTAGAATTACACCATAATTGCT 
# TATCTGTCTCGACGT[AGCT]{5}GGCGAGGTGATCTAGAATTACACCATAATTGCT 
# TATCTGTCTCGACGT[AGCT]{5}GGCTTCAACAACTAGAATTACACCATAATTGCT 
# TATCTGTCTCGACGT[AGCT]{5}GGCGATGGTTTCTAGAATTACACCATAATTGCT 

#RC: AGCAATTATGGTGTAATTCTAGTTGCAAAGGCC[ACGT]{5}ACGTCGAGACAGATA     
#RC: AGCAATTATGGTGTAATTCTAGGCTTGAATGCC[ACGT]{5}ACGTCGAGACAGATA     
#RC: AGCAATTATGGTGTAATTCTAGATCACCTCGCC[ACGT]{5}ACGTCGAGACAGATA     
#RC: AGCAATTATGGTGTAATTCTAGTTGTTGAAGCC[ACGT]{5}ACGTCGAGACAGATA     
#RC: AGCAATTATGGTGTAATTCTAGAAACCATCGCC[ACGT]{5}ACGTCGAGACAGATA     
       
        
###################################################################################################################################################################################################

for(my $i=0; $i<=$#files; $i=$i+1) {
        my $temp = $files[$i];
        ##next unless $temp =~ m/^CHe/; 
        next unless $temp =~ m/\.fastq/;
        my $input = "$input_g/$temp";
        my $out = "$output_g/$temp/$temp";
        $out =~ s/\.fastq\.gz//g or die;
        #$out =~ s/\.fastq//g or die;
        &myMakeDir($out);
        print "$temp\n";
        system("zcat  $input  |  grep -oP GTCTCGACGT[AGCT]{5}GGCCTTTGCAACTAGAATTACACCAT     > $out.spikein-1.seq");
        system("zcat  $input  |  grep -oP GTCTCGACGT[AGCT]{5}GGCATTCAAGCCTAGAATTACACCAT     > $out.spikein-2.seq");
        system("zcat  $input  |  grep -oP GTCTCGACGT[AGCT]{5}GGCGAGGTGATCTAGAATTACACCAT     > $out.spikein-3.seq");
        system("zcat  $input  |  grep -oP GTCTCGACGT[AGCT]{5}GGCTTCAACAACTAGAATTACACCAT     > $out.spikein-4.seq");
        system("zcat  $input  |  grep -oP GTCTCGACGT[AGCT]{5}GGCGATGGTTTCTAGAATTACACCAT     > $out.spikein-5.seq");

}

###################################################################################################################################################################################################






say   "Done ......";






#####
