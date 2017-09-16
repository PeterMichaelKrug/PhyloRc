#!/usr/bin/perl -w
use strict;
use warnings;
my @infile=qw(RAxML_bipartitions.Alignment.phylip.partition_pergene.txt.excludingrps4.txt.bip.tre.list RAxML_bipartitions.Alignment.phylip.partition_pergene.txt.excludingnad2.txt.bip.tre.list);
#my $infile2="";
my $first=1;#first file
my @table;

foreach my $file (@infile){
   my $countr=0;
    open IN,"<",$file;
    while (<IN>){
         
           #first file defines age and leafs, then bips
           next if (/Age,BS,node/gi);
           if (/(\d+.\d+)\,(\d+)\,(.+)/gi){
            
            if($first){
                #my $string="$1,$3,$2";
                my $string="$1,$2";
                push @table, $string;
                
            }else{
                
                $table[$countr]=$table[$countr].",".$2;
                print $table[$countr],"\n";
                
            }
            $countr++;
            print "\nNODE: ",($countr+1),"\n";
           }
    }
    $first=0;
}
#foreach (@table){
#    print$_,"\n";
#}