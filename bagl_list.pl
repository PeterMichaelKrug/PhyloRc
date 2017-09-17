#!/usr/bin/perl -w
use strict;
use Cwd;
#load all lists:
my $file;#="partition_cp_mt_coding_noncoding.txtreplicates/matK_vs_mt_noncoding/RAxML_bipartitions.matK_vs_mt_noncodingbipplot.list"; #make @list
my $out;#="reduced_ages.txt";

print "\n\ntool bagl_list.pl started with \n";
if(@ARGV){
    $file=$ARGV[0];
    $out=$ARGV[1];
}
print "$file\n\n";
open FILE,"<",$file or die $!;
my (@ages,%age_support,$median);
while (<FILE>){
    if(/(\d+\.\d+)\,(\d+)\,(.+)/ig){
        #print $1,"\t\t",$2,"\n";
        push @ages, $1;
        $age_support{$1}=$2;
    }
}
@ages=sort {$b <=> $a} @ages;
#print "Ages ",join "\n",@ages,"\n";
$median=&median(@ages);
print "\nMedian: ",$median;
my($shallow,$deep,$snodes,$dnodes,$nodes,$all);
foreach (@ages){
    $nodes++;
    if ($_<$median){
        $deep+=$age_support{$_} if ($age_support{$_}>50);
        $all+=$age_support{$_} if ($age_support{$_}>50);
        #print $deep,"\t\t",$age_support{$_},"\n";
        $dnodes++;
        
    }else{
        $shallow+=$age_support{$_} if $age_support{$_}>50;
        $all+=$age_support{$_} if ($age_support{$_}>50);
        $snodes++;
        
    }
    #print "all :\t\t:",$all,"\n";
}

my $name=$file;
$name=~s/RAxML_bipartitions(\w+)bipplot.list/$1/g; #funxt nicht
open BAGL,'>>',$out;
print BAGL "marker,gradient: ";
print BAGL $name,",";

print "\ndeep: ",join " ",$deep,"deep nodes: ",$dnodes,"\n";
print "quotient ",$deep/$dnodes,"\n";
print "\nshallow: ",join " ",$shallow,"shallow nodes: ",$snodes,"\n";
print "quotient ",$shallow/$snodes,"\n";
#print "all: ",$all," nodes ",$nodes,"\n";
print BAGL "(($shallow/$snodes)-($deep/$dnodes))/($all/$nodes)";
print BAGL (($shallow/$snodes)-($deep/$dnodes))/($all/$nodes),"\n";
print "cwd: ",cwd;
close BAGL;
print "\nWritten to $out\n";
# subroutines
sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}