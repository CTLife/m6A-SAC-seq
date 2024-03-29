#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################



my $input1="2-getFraction";
my $out="7-getMotifs";


###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Running ......";

sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die;       }
}
&myMakeDir($out);

opendir(my $FH_input1_g, $input1)  ||  die;
my @files = readdir($FH_input1_g);
###################################################################################################################################################################################################



 
        
###################################################################################################################################################################################################

for(my $i1=0; $i1<=$#files; $i1=$i1+1) {
        my $temp1 = $files[$i1];
        next unless $temp1 =~ m/\.bed$/;
        print "$temp1\n";
        $temp1 =~ s/\.bed$// or die;

        open(my $FH1_g,  "$input1/$temp1.bed" )  ||  die;
        my @lines1 = <$FH1_g> ;
        open(fileA,  ">", "$out/$temp1.txt" )  or  die;

	for(my $i=0; $i<=$#lines1; $i=$i+1) {
       	 my $temp = $lines1[$i];
       	 $temp =~ m/\.\.\.\.\.\.([AGCT]{5})\s/i  or  die "##$temp##";                                
       	 my $motif = $1;	   
       	 print   fileA    "$motif\n";       
	}
 
        close(fileA);
        system( "sort    $out/$temp1.txt        |   uniq   -c     | sort -n     -k 1   -r  -o   $out/$temp1.number" );
        system( "Rscript  7-histogram.R  -i $out/$temp1.number  -o $out/$temp1.pdf" );
}

###################################################################################################################################################################################################






say   "Done ......";






#####
