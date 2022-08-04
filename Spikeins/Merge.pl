#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Running ......";

sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die;       }
}
my $output_g = "Merge";
&myMakeDir($output_g);

opendir(my $FH_input1_g, "Spikeins_RC")  ||  die;
my @files1 = readdir($FH_input1_g);
###################################################################################################################################################################################################




 
###################################################################################################################################################################################################

for(my $i=0; $i<=$#files1; $i=$i+1) {
        my $temp1 = $files1[$i];
        next unless $temp1 =~ m/\.R1$/;
        my $temp2 = $temp1;
        $temp2 =~ s/\.R1$/.R2/  or  die;
        
        my $file_A1 = "Spikeins_RC/$temp1/$temp1.spikein-1.seq";
        my $file_A2 = "Spikeins_RC/$temp1/$temp1.spikein-2.seq";
        my $file_A3 = "Spikeins_RC/$temp1/$temp1.spikein-3.seq";
        my $file_A4 = "Spikeins_RC/$temp1/$temp1.spikein-4.seq";
        my $file_A5 = "Spikeins_RC/$temp1/$temp1.spikein-5.seq";
                
        my $file_B1 = "Spikeins/$temp2/$temp2.spikein-1.seq";
        my $file_B2 = "Spikeins/$temp2/$temp2.spikein-2.seq";
        my $file_B3 = "Spikeins/$temp2/$temp2.spikein-3.seq";
        my $file_B4 = "Spikeins/$temp2/$temp2.spikein-4.seq";
        my $file_B5 = "Spikeins/$temp2/$temp2.spikein-5.seq";
        
        open(A1, "<",  $file_A1 ) or die;
        open(A2, "<",  $file_A2 ) or die;
        open(A3, "<",  $file_A3 ) or die;
        open(A4, "<",  $file_A4 ) or die;
        open(A5, "<",  $file_A5 ) or die;       
        open(B1, "<",  $file_B1 ) or die;
        open(B2, "<",  $file_B2 ) or die;
        open(B3, "<",  $file_B3 ) or die;
        open(B4, "<",  $file_B4 ) or die;
        open(B5, "<",  $file_B5 ) or die;
        
        my @lines_A1 = <A1>; 
        my @lines_A2 = <A2>; 
        my @lines_A3 = <A3>; 
        my @lines_A4 = <A4>; 
        my @lines_A5 = <A5>; 
        my @lines_B1 = <B1>; 
        my @lines_B2 = <B2>; 
        my @lines_B3 = <B3>; 
        my @lines_B4 = <B4>; 
        my @lines_B5 = <B5>; 
        
        for(my $i1=0; $i1<=$#lines_A1; $i1=$i1+1) {
            my $L1 = $lines_A1[$i1];
            $L1 =~ s/\n//  ;
            my $revcomp = reverse($L1);
            $revcomp =~ tr/ATGCatgc/TACGtacg/;
            $lines_A1[$i1] = "$revcomp\n";
        }
        
        for(my $i1=0; $i1<=$#lines_A2; $i1=$i1+1) {
            my $L1 = $lines_A2[$i1];
            $L1 =~ s/\n//  ;
            my $revcomp = reverse($L1);
            $revcomp =~ tr/ATGCatgc/TACGtacg/;
            $lines_A2[$i1] = "$revcomp\n";
        }
        
        for(my $i1=0; $i1<=$#lines_A3; $i1=$i1+1) {
            my $L1 = $lines_A3[$i1];
            $L1 =~ s/\n//  ;
            my $revcomp = reverse($L1);
            $revcomp =~ tr/ATGCatgc/TACGtacg/;
            $lines_A3[$i1] = "$revcomp\n";
        }     
        for(my $i1=0; $i1<=$#lines_A4; $i1=$i1+1) {
            my $L1 = $lines_A4[$i1];
            $L1 =~ s/\n//  ;
            my $revcomp = reverse($L1);
            $revcomp =~ tr/ATGCatgc/TACGtacg/;
            $lines_A4[$i1] = "$revcomp\n";
        }         
        for(my $i1=0; $i1<=$#lines_A5; $i1=$i1+1) {
            my $L1 = $lines_A5[$i1];
            $L1 =~ s/\n//  ;
            my $revcomp = reverse($L1);
            $revcomp =~ tr/ATGCatgc/TACGtacg/;
            $lines_A5[$i1] = "$revcomp\n";
        }  
        
        $temp1 =~ s/\.R1$//  or  die;
        &myMakeDir("$output_g/$temp1");
        open(C1, ">",  "$output_g/$temp1/$temp1.spikein-1.seq" ) or die;
        open(C2, ">",  "$output_g/$temp1/$temp1.spikein-2.seq" ) or die;
        open(C3, ">",  "$output_g/$temp1/$temp1.spikein-3.seq" ) or die;
        open(C4, ">",  "$output_g/$temp1/$temp1.spikein-4.seq" ) or die;
        open(C5, ">",  "$output_g/$temp1/$temp1.spikein-5.seq" ) or die;
        
        print C1 @lines_A1;
        print C2 @lines_A2;
        print C3 @lines_A3;
        print C4 @lines_A4;
        print C5 @lines_A5;
        
        print C1 @lines_B1;
        print C2 @lines_B2;
        print C3 @lines_B3;
        print C4 @lines_B4;
        print C5 @lines_B5;
        
}

###################################################################################################################################################################################################






say   "Done ......";






#####
