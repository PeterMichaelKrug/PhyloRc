#!/usr/bin/perl -w

#todo argument bs-threshold (50/70) see Hillis & Bull (1993)
my $ok=1; my $noRs=0;
my $time=time;
use strict;
use warnings;
use Number::Range;
use Data::Dumper;
use File::Path(qw(make_path remove_tree));
use List::Compare;
use Statistics::Descriptive;
use Bio::Phylo::IO;
use Cwd;
use Bio::TreeIO;
use Bio::Tree::NodeI;
use List::Compare;
##################################################
my ($template);

#++++++++++++++++++
if (@ARGV){
    $template=$ARGV[0];
}
my %options=&get_info($template);
print Dumper(\%options);
#get values
my $alnfile=$options{alignment_file};
my $partfile=$options{partition_file};
my $reftop=$options{reference_topology};
my $bootslarge=$options{nboots_large};
my $bootsmall=$options{nboots_small};
my $replicates=$options{nreplicates};
my $rax=$options{path_to_raxml};
my $threads=$options{num_threads};
my $model=$options{substitution_model};
my $cleanup=$options{cleanup};
my $outgroup=$options{outgroup};
#exit;
#*********************

#other variables
my ($reftop_forest);
#make raxml options 

my $opt="-f a ";
my $xp=" -x 12345 -p 12345 ";
my $raxmlpath="$rax -T $threads -m $model -o $outgroup";
#my $outgroup="-o Sphagnum_capillifolium_,Takakia_lepidozioides_";
#++++++++++++++++++ delete data from prev runs
#goto BAGL_LIST;
system "rm RAxML_*" if ($cleanup);

#parse topology, do some calculations
open TREE,"<",$reftop;
while(<TREE>){$reftop_forest.=$_ };
my $forest = Bio::Phylo::IO->parse(-format => 'newick',-string => $reftop_forest);
my $reftop_tree= $forest -> first;
#calculate a lineage through time plot via Bio::Phylo
my $ltt_ref=$reftop_tree -> calc_ltt;
my @ltt=@{$ltt_ref};
unshift @ltt, $alnfile;
&do_ltt(@ltt);          #prints tab delimited file
print "made a ltt plot";
#**********************
#calculate the pybus gamma statistic via Bio::Phylo
my $pybusgamma=$reftop_tree->calc_gamma();
print "\nThe Pybus_Gamma statistic:",$pybusgamma,"\n";
#**********************

#first make a plain phylogeny with all partitions;
my $raxmltop="$raxmlpath $xp -# $bootslarge $opt -s $alnfile -q $partfile -n $alnfile\.$partfile\.tre";# >/dev/null 2>&1";
print "\n\nCalculating a master phylogeny with all partitions\n";
print "with $raxmltop\n";
system $raxmltop if ($ok);
 #exit;         #produces RAxML_bootstrap\.$alnfile\.$partfile.tre
open ALI,"<", $alnfile or die "$!";
my $first=1; my($taxanum,$length);
my @RefSeq; my @taxlist;
while (<ALI>){
    #get first line: $taxa and $length
    if($first){
        ($taxanum,$length)=split /\s/;
        print "Taxa: ",$taxanum," and length: ",$length,"\n";
        $first=0;
    }elsif(/^(\S+)\s(.+)/){
        #print "A taxon obviously: ",$1,"\n";
        push @taxlist, $1;
        my @seq=split // ,$2;
        #print "\tand its sequence(first 100 only): \n";
        #for (1300..1400)
            #{print $seq[$_]}
        #print "\n";
        push  @RefSeq,\@seq;
    }
}


open PART,"<", $partfile or die "$!";
my (@lines,%partitions, %nucleotides,%positions,%nucperpos);
#evaluate partitions
while(<PART>){
		
    push @lines,$_;
    if(/(BIN),\s*(\S+)\s*\=(.+)/i){
        #print $1,"\t",$2,"\n";         #put values into $5,$6 etc if needed
        my $pname=$2;                   #get name of partition
        print "found ",$pname, " ";
        my $num=$3;
        chomp $num;
        $num=~s/\s//ig;
        $num=~s/-/../ig;
        #print $num,"\n und dann \n";
        my $range=Number::Range->new($num);
        my @numbers=$range->range;
        #calc length of Alignment = positions
        print "\tlength in Characterpositions ",$#numbers+1,", ";      #again @number starts @ 0
        $positions{$pname}=$#numbers+1;
        #count nucs
        my $nucs;
        foreach $num(@numbers){        #partitions begin with 1, array with 0
            my $nucs2;
            #print $num, " has ";
            for (0..$taxanum-1){
                #print $_;       #[0..63]
                #print $RefSeq[$_][$num-1];
                #counting
                if ($RefSeq[$_][$num-1] =~ /[ATGCMRWSYKVHDBN01\?]/i){      #get Ambiguities also, everything except "\-" does not work ok.
                    #print "0";
                    $nucs++;
                    $nucs2++;           #count nucs per site, into hash{position}=[\#nucs2]
                }else{
                    #print;     #[----]
                }
                $nucperpos{$num}=[$nucs2];          #numbers of nucleotides per position
            }
            #print"\n";
        }
        print " and ",$nucs," nucelotides in total\n";
        $nucleotides{$pname}=$nucs;
        $partitions{$pname}=\@numbers;
    }
};

my $allpos=scalar keys %nucperpos;
#print "\nnucperpos is now ",$allpos," in length\n";
#sort partitions by length, descending, printing to .rankedstats.txt
my @rank=sort sorting keys %nucleotides;

open STATS,">","$partfile.\.rankedstats.txt";
print STATS "marker\t#nucleotides\t#positions\t#nuc/#pos\n";
foreach (@rank){
    print STATS $_,"\t", $nucleotides{$_},"\t",$positions{$_},"\t",$nucleotides{$_}/$positions{$_},"\n"; #(get gappiness:$nucleotides{$_}/$positions{$_})/#taxa
}
close STATS;
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#start resampling
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#create folder
my $folder=$partfile."replicates";      #folder to pass to raxml via -w (makeshort!)
rmdir $folder if ($ok);
mkdir $folder if ($ok);
chdir $folder;
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#           start processing partitions:
#   (1)resampling,
#   (2)calc bs_age_gradient, todo: pass to bs_age_gradient.pl ##ok!
#   (3)calc R, todo: pass trees with support to doallcomps.pl (or a new version of that) ##ok!!
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
my @repcoll;            #collect all replicates here to start raxml with
my @pearsons;           #get all
my @perpart;            #list all bootstrap_age results
#outer: @rank contains nucleotides in ASCENDING order wrt #nuc  #reason: .. cannot count backwards
#outer loop:each smaller partition
        my $round=1;
 
for (my $i=0;$i<$#rank;$i++){   #do not need <= because last one does not need no comparison
    
    print "Comparing smaller alignment ",$rank[$i]," with $nucleotides{$rank[$i]} nucleotides and.. \n";
    #start raxml on smaller alignment
    #extract $partitions{$rank[$i]} from Alignment to new file
    print "Calculating plain phylogeny for ",$rank[$i],"\n";
    system "cp ../$alnfile $alnfile";            #as we need the Alignment file here
    my @sm=@{$partitions{$rank[$i]}};
    my $toexclsm=$rank[$i]."\.excl.txt";
    my $uniquestrgsm;
    for(1){     #make sub of this
        open OUT,">", $toexclsm;      
        my @all=1..$allpos;
        my @repsm=sort {$a <=> $b} @sm;           #but we want excludefiles, so..
        my $lc = List::Compare->new('-u', \@all, \@repsm);
        my @uniquesm = $lc->get_unique;           #returns only those values unique to @all, meaning all exluding the reps
        @uniquesm = sort {$a <=> $b} @uniquesm;
        $uniquestrgsm = ranges(@uniquesm);
        print OUT $uniquestrgsm,"\n";               #exclusion file must look like 1-1234 1236-2345 without anything else!!
        close OUT;
    }
    #print "\n\=\^\.\^\=  Excluding $uniquestrgsm \n";
    my $rax_smaller_excl="\.\./".$raxmlpath." -E ".$toexclsm." -n somethingnotneeded -s $alnfile";# >/dev/null 2>&1";          #output will be saved in this level
    print "\n\tExecuting $rax_smaller_excl\n";
    system $rax_smaller_excl if ($ok);
    my $rax_smaller_cmd="\.\./".$raxmlpath.$xp.$opt." -s ".$alnfile."\.".$toexclsm." -n ".$toexclsm."\.tre -\# ".$bootsmall." >/dev/null 2>&1";
    #print "\n\tExecuting $rax_smaller_cmd\n";
    system $rax_smaller_cmd if ($ok);
#now count node frequency of smaller set on fixed topology
    my $fixedsm="\.\./".$raxmlpath." -f b -t \.\./".$reftop." -z RAxML_bootstrap\.$toexclsm\.tre -n ".$toexclsm."\.tre\.fixplot";#  >/dev/null 2>&1";
    print "Did ", $fixedsm;
#call bsagegrad.pl on smaller set
    if ($round eq 1){                       #we need this only on first iteration
        #$perpart[0]="$folder/RAxML_bipartitions\.$toexclsm\.tre\.fixplot";
        my $bsagegrad="perl \.\./bs_age_gradient.pl \.\./$reftop RAxML_bipartitions\.$toexclsm\.tre\.fixplot";
        my $list="$folder\/RAxML_bipartitions\.$toexclsm\.tre\.fixplot.list";
        push @perpart,$list;
        print "HERE: ",$bsagegrad;
        system $bsagegrad;
        print "\n\t\tjust calculated $toexclsm\.fixplot.list\n";
        #exit;
    }
    #print @perpart;
    system $fixedsm ;#if ($ok); #written to RAxML_bipartitions.matK.excl.txt.tre.fixplot
    my $treefile_sm="RAxML_bipartitions\.$toexclsm\.tre";
    my $treefile_smfix="RAxML_bipartitions\.$toexclsm\.tre\.fixplot";
    my @smR=&sub_calcR($treefile_sm);
    my @smRfix=&sub_calcR($treefile_smfix);   
    my ($meanSM, $meanRSM)=(@{$smR[0]},@{$smR[1]});
    print "\n\tmean: ",$meanSM,"\tmeanR: ",$meanRSM,"\n";
    my ($meanSMfix, $meanRSMfix)=(@{$smRfix[0]},@{$smRfix[1]});
    print "\n\tmean: ",$meanSMfix,"\tmeanR: ",$meanRSMfix,"\n";
    #exit;
##inner loop: each bigger partition
    for (my $j=$i+1;$j<=$#rank;$j++){
        #print "to be compared: $rank[$j] with $nucleotides{$rank[$j]} nucleotides\n";
        #make replicates with $nucleotides{$rank[$i]}
        
        my $tsfref = $partitions{$rank[$j]};            #to sample from 
        #~~~~~~~~~~~~~~~~~~~~~~~
        #change directory
        #~~~~~~~~~~~~~~~~~~~~~~~
        my $dir=$rank[$i]."_vs_".$rank[$j];
        if ($ok){
            remove_tree($dir) ;                           #savely delete directory if exsist
            mkdir $dir,0755 or die $!;
        }
        chdir $dir;                                         #go there
        print "here we are: "; cwd;
        system "cp ../../$alnfile $alnfile";
        my (@bootcollect,@bipartcollect,@fixcollect);       #jam memory
        for(my $s=1;$s<=$replicates;$s++){                  #next replicate
            my @tsf=@{$tsfref};                             #@tsf contains positions from bigger marker
            my $size = 0; my $sum=0;
            print "replicate No. $s\n";
            my $partition_rep="part_".$rank[$i]."_vs_".$rank[$j]."_bs_".$s."\.txt";
            #print "\$partition_rep\n";
            my @rep;                                        #put rand set positions in
            BLOCK:for (my$k=0;;$k++){                       #infinite loop!!!end with last BLOCK
                
                my $currRandNr=int(rand(scalar(@tsf)));     #count rand pos per replicate -> decide if
                my $currpos=splice @tsf, $currRandNr, 1;    #splice currpos from bigger alignment
                push @rep,$currpos;
                #calc #nucs in @reps
                $size=@{$nucperpos{$currpos}}[0];
                if($size){
                    $sum+=$size;
                }else{
                    #print "\nempty position in $currpos\n";
                }
                #random machine to stop sampling
                    #my $randstop=$nucleotides{$rank[$i]}-$sum;					
                    #my $randnum=rand(1);
                    #my $tempstop=$randstop/$nucleotides{$rank[$j]};
                #print $size;                       #$size has now the #nucs in currpos
                if ($sum>=$nucleotides{$rank[$i]}){         #quit always if bigger
                    #print "writing to $partition_rep\n";
                    #write partition file with bootstraps
                    open OUT,">", $partition_rep;
                    
                    my @all=1..$allpos;
                    my @rep=sort {$a <=> $b} @rep;           #but we want excludefiles, so..
                    my $lc = List::Compare->new('-u', \@all, \@rep);
                    my @unique = $lc->get_unique;           #returns only those values unique to @all, meaning all exluding the reps
                    @unique = sort {$a <=> $b} @unique;
                    my $uniquestrg = ranges(@unique);
                    print OUT $uniquestrg,"\n";         #exclusion file must look like 1-1234 1236-2345 without anything else!!
                    close OUT;
                    chmod 0666,$partfile;
                    push @repcoll,"$dir\\$partition_rep";
                    last BLOCK;
                    
                }#elsif (zufallsgenerator?!?){
                #    do exactly the same;
                #}
            }
        #kick racksemel to reduce alignment and...
            my $reducer="\.\./\.\./".$raxmlpath." -E ".$partition_rep." -n somethingnotneeded -s $alnfile  >/dev/null 2>&1";
            #print "reducing via \t",$reducer,"\n";
            system "$reducer" if ($ok);
        #...generate bootstrap replicates
            my $replicator="\.\./\.\./".$raxmlpath.$xp." -f a -s ".$alnfile.".".$partition_rep." -n ".$partition_rep."_rep.tre -\# ".$bootslarge." >/dev/null 2>&1";
            #print "replicating now ",$replicator,"\n";
            system $replicator if ($ok);
        #count node frequencies against fixed topology 
            my $fixed="\.\./\.\./".$raxmlpath." -f b -t "."\.\./\.\./".$reftop
                        ." -z RAxML_bootstrap\.".$partition_rep."_rep.tre -n ".$partition_rep."\.fixplot >/dev/null 2>&1";
            system $fixed;
            #print $fixed; 
            push @bootcollect, "RAxML_bootstrap\.".$partition_rep."_rep.tre ";
            push @bipartcollect, "RAxML_bipartitions\.".$partition_rep."_rep.tre ";
            push @fixcollect, "RAxML_bipartitions\.".$partition_rep."\.fixplot";
        }
        #concat treefiles
        #print @bootcollect;
        my $catboot=join ' ', @bootcollect;     #not used at all
        print "\n catboot:",$catboot,"\n";  
        my $catbip=join' ', @bipartcollect;
        my $catfix=join' ',@fixcollect;
    #calculate Rs via sub_calcR:
    #1) classic, each one on its own ml-topology
        my @calcedR=&sub_calcR(@bipartcollect);
            print "\n\tmean: ",@{$calcedR[0]},"\n\tmeanR: ",@{$calcedR[1]},"\n\n";
    #2) new, this time node frequencies plottet to a given topology
        my @calcedRfix=&sub_calcR(@fixcollect);
            print "\n\tmean: ",@{$calcedRfix[0]},"\n\tmeanR: ",@{$calcedRfix[1]},"\n\n";
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #start sub_compareR here for @smR and all the @calcedRs/@caledRfixes
        #1)
        my @Rs_to_comp=@{$smR[0]};
        my @Ms_to_comp=@{$smR[1]};
        push @Rs_to_comp, @{$calcedR[0]};
        push @Ms_to_comp, @{$calcedR[1]};
        my $statsoutfile="\.\./$alnfile\.statistics_out.txt";
        &sub_compRs([@Ms_to_comp],[@Rs_to_comp],$statsoutfile,$dir,$round);       #should it be returning something?
        #2)
        my @Rfixs_to_comp=@{$smRfix[0]};
        my @Mfixs_to_comp=@{$smRfix[1]};
        push @Rfixs_to_comp, @{$calcedRfix[0]};
        push @Mfixs_to_comp, @{$calcedRfix[1]};
        my $statsoutfile_fix="\.\./$alnfile\.statistics_fix_out.txt";
        &sub_compRs([@Mfixs_to_comp],[@Rfixs_to_comp],$statsoutfile_fix,$dir,$round);       #should it be returning something?
        #print "means ",join ", ",@Rs_to_comp;
        #print"\n";
        #print "meanRs ",join ", ",@Ms_to_comp;
       
        open STAT,"<",$statsoutfile or die $!;
        my @rstat=<STAT>;
        #print @rstat;
        #last;
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        #print "cat ".$catboot." > ".$dir."\.allboot";
        system "cat ".$catboot." > ".$dir."\.allboot";
        #concat bipartition trees also
        #print "cat ".$catbip." > ".$dir."\.allbip";
        system "cat ".$catbip." > ".$dir."\.allbip";
        #
        #call raxml to plot node frequencies to topology;
        my $bstree="RAxML_bootstrap\.$alnfile\.$partfile.tre";
        #call Raxml to plot node freqs to $bstree
        my $rax_bip_cmd="\.\./\.\./".$raxmlpath." -f b -t "."\.\./\.\./".$reftop." -z ".$dir."\.allboot -n ".$dir."bipplot";
        print "\n",$rax_bip_cmd,"\n";
        system "$rax_bip_cmd" ;#if $ok;
        print "\n\n Ready plotting bipartition frequencies in $dir\n\n";
        #now repeat for each rep to get sdev and 95%CI from replicates
        my @repstat;
        for (0..$replicates-1){
            #print "doing ",$bootcollect[$_],"\n";
            my $rax_per_rep="\.\./\.\./".$raxmlpath." -f b -t "."\.\./\.\./".$reftop." -z ".$bootcollect[$_]." -n ".$dir."repplot".($_+1)."  >/dev/null 2>&1";
            system $rax_per_rep;
            #now we have files like .$dir."repplot".($_+1) ~ RAxML_bipartitions.matK_vs_cp_CDSrepplot10
            #open each, tableize them
            my $rep="RAxML_bipartitions.".$dir."repplot".($_+1);
            push @repstat, $rep;
        }
        print getcwd(),"\n";
        #print join "\n ",@repstat;
       
            my ($repstattable)=&sub_rep_stat(@repstat);
            
        
        #call Raxml to calc Pearson corr coeff
        #*************************************
            my $pearsoncmd="\.\./\.\./".$raxmlpath." -f m -t ".$dir."\.allboot -z \.\./\.\./RAxML_bootstrap.$alnfile.$partfile.tre -s $alnfile -n $alnfile.$dir.pearson.txt >/dev/null 2>&1";
            print $pearsoncmd,"\n";
            
            system $pearsoncmd;
            my $pear="RAxML_info.$alnfile.$dir.pearson.txt";
            print $pear;
            my ($pearson)=&sub_pearson($pear); 
            print "\n\nPearson: ",$pearson;
            push @pearsons,$pearson;
        #***********************************
        #call bs_age_gradient.pl with ARGS:
                    #ARGV[0]=chronogram [pruned_chronogram_mosses_nwk.tre]
                    #ARGV[1]=tree with plotted node frequencies
        my $bsagegrad="perl \.\./\.\./bs_age_gradient.pl \.\./\.\./$reftop RAxML_bipartitions\.".$dir."bipplot";
        #cycle bs_age_gradient.pl through bipartition files from each replicate
        my $list="$folder\/$dir\/RAxML_bipartitions\.".$dir."bipplot\.list";
        push @perpart,$list if ($round<2);
        print "\n",$bsagegrad,"\n\n";
        system $bsagegrad ;#if $ok;
        #uncomment all system,mkdir,print OUT etc
        #last;                   #rmove this!!
        chdir "..";             #done with files in this folder, now create and change to next folder/marker
        
        #print $tsf[0]," to ",$tsf[-1],"\n";         #hash %partitions contain references to
        # get the duration so far
        my $now=time;
        my $duration=convert_seconds_to_hhmmss($now-$time);
        print "\n*****************************************";
        print "\nduration ",$now-$time," secs. or $duration";
        print "\n*****************************************\n";
        print "now calculating age dependent bootstrap level";
        #put Bootstrap_alignment_partition.pl here
        
    }
    
   $round++;
   #shrink to smaller partition size only
   last if($noRs); #stop here!!! 
    #print "++++ make replicates with $nucleotides{$rank[$i]} nucleotides for $rank[$i]\n";
}
chdir "..";
BAGL_LIST:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       use sub table to generate an output file for all files if $round <1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print"\n\n*********************\nstarting bagl_list.pl with:\n";
#unlink "$folder\.append.txt" or die "$!";
print @perpart,"\n";exit;
foreach(@perpart){
    print "\n $_";
    print "cwd: ",cwd;
    print "\nAnalysing $_\n";
    my $bagl_cmd="perl bagl_list.pl $_ $folder\.append.txt";
    #print $bagl_cmd,"\n";
    system "$bagl_cmd";
    #exit;
}
my $pearstring=join ",",@pearsons;
print "\n\n all pearsons",$pearstring;

my $now=time;
my $duration=convert_seconds_to_hhmmss($now-$time);
print "\n*****************************************";
print "\nall calculations took ",$now-$time," secs. or $duration";
print "\n*****************************************\n\n";
#print "\n\n$catboot";
#********************************************************************************************
#       subroutines
#       1)sorting
#       2)ranges
#       3)sub_pearson
#       4)sub_calcR
#       5)sub_compRs
#       6)do_ltt
#       7)convert_seconds_to_hhmmss
#       8)table
#       9)sub_rep_stat
#       10)get_info
#********************************************************************************************
sub sorting{
    $nucleotides{$a}<=>$nucleotides{$b}   #values in alphabetical order, ASCENDING!!!
    or                          #may happen: same # nucleotides
    $a cmp $b                   #keys in Ascii order
}
sub ranges {
    my @vals = @_;
    my $min = $vals[0];
    my $max;
    my @list;
    for my $i (0 .. (scalar(@vals)-2)) {
        if (($vals[$i+1] - $vals[$i]) != 1) {
            $max = $vals[$i];
            push @list, ($min == $max) ? "$max-$max" : "$min-$max";         #will print single numbers as 123-123 and ranges as 124-127;
            $min = $vals[$i+1];
        }
    }
    $max = $vals[-1];
    push @list, ($min == $max) ? "$max-$max" : "$min-$max";
    return join ' ', @list;
}

sub sub_pearson{            
    my $pearson;
    open PEA,"<",$_[0] or die "$!";
    while(<PEA>){
       if (/^Pearson\:\s(0\.\d+)/gi){
        print$1;
        $pearson=$1;
       }
    }
    $pearson;
}

sub sub_bs_age_grad{
    print @_;
}

sub sub_calcR{
    my (@meanc,@meanRc);
    my @treefiles=@_;
    #print "Treefiles: ",@treefiles;
    foreach my $treefile(@treefiles){
        #print"\n\n Printing subcalcRs\n$treefile\n";
        open( TREE, $treefile ) or die "Cannot open $treefile\n";

        my @treelines = <TREE>;
        my $theTree = "";
        foreach my $line (@treelines) {
            $theTree = $line;   #raxml bipartition trees are plain newick, one line, one tree
            last;
        }
        my @bp = ();
        @bp = $theTree =~ /([0-9]{1,3})\:/g;
        #searches for 5, 50 and 100, ignores branch length 0.\d+
        my $bp_str = join ", ",@bp;
        
        #print "the bootstraps: ",$bp_str,"\n" ;
        
        #*******************************************************************************
        #++++++++++++++++++++++++++ calculate support statistics:+++++++++++++++++++++++
        #*******************************************************************************
        
        #sum all  support values
        my $sum  = 0;
        my $sumR = 0;
        
        #print "values\tsingleRs\n\n";
        foreach my $val (@bp) {
            #if ($isMB) { $val *= 100 } #transform MB PPs to percentages
            $sum += $val;
        
        #	print $val, "\t";
            if ( $val >= 50.0 ) {
                my $singelR = ( $val - 49.0 ) / 51.0;
                $sumR += $singelR;
        
        #		print "$singelR\n";
            }
            else {
        
        #		print "-\n";
            }
        }
        #divide by ntax-2
        #=> average
        my $mean  = $sum /  ( ( scalar @bp ) * 100 );
        my $meanR = $sumR / ( ( scalar @bp ) );
        push @meanc,$mean;
        push @meanRc,$meanR;        
        #open RF, ">>$R_outputfile" or die "Cannot open $R_outputfile\n";
        #print RF "$mean\t$meanR\t$treefile\n";
    }
    my @out=([@meanc], [@meanRc]);
    return @out;
    #return @meanc, @meanRc;
}
sub sub_compRs{
    #called like &sub_compRs(@Ms_to_comp,@Rs_to_comp,$statsoutfile);
    #print"\n ", join ": ",@_; #gives two refs and a filename
    my ($Ms_to_comp,$Rs_to_comp,$stat_outputfile,$dir,$round)=@_;
    my ($smaller_marker,$bigger_marker);
        if ($dir =~ /(\w+)_vs_(\w+)/gi){
        ($smaller_marker,$bigger_marker)=($1,$2);
    }
   
    #my @allRs;
    my $meanSm=${$Ms_to_comp}[0];
    my $RSm=${$Rs_to_comp}[0];
    my $statM = Statistics::Descriptive::Full->new();
    my $statR = Statistics::Descriptive::Full->new();
    
    #my @diffsMean=(); my @diffsR=();
    #print $meanSm, $RSm, $treeFile, "\n";
    for ( my $i = 1 ; $i <= $#{$Ms_to_comp} ; $i++ ) {          #first line is already used
          my $meanLa=${$Ms_to_comp}[$i];
          my $RLa=${$Rs_to_comp}[$i];
    
    #calc diffs between [0] (meanSm) and Rest [1]..[n] (meanLa)
    #=>structure in first gene (smaller dataset) in file name minus structure in second gene (larger dataset) in  file name
        my $diffMean = $meanSm - $meanLa;
        my $diffR    = $RSm - $RLa;
    
        #print "$diffMean\n";
        $statM->add_data($diffMean);
        $statR->add_data($diffR);
    }
    my $varM   = $statM->variance();
    my $qu     = $varM/ $statM->count();
    my $seM    = sqrt($qu);
    my $lowBdM = $statM->mean() - ( 1.96 * $seM );
    my $upBdM  = $statM->mean() + ( 1.96 * $seM );
    
    my $varR = $statR->variance();
    $qu = $varR / $statR->count();
    my $seR    = sqrt($qu);
    my $lowBdR = $statR->mean() - ( 1.96 * $seR );
    my $upBdR  = $statR->mean() + ( 1.96 * $seR );
    
    print "all Rs read & significance tested.\n";
    open ST, ">>$stat_outputfile" or die "Cannot open $stat_outputfile\n"; ###appending 
    print ST "\n$round\t";
    print ST "$dir\t";
    print ST "mean\t", sprintf( "%.4f", $statM->mean() ), "\t",
      sprintf( "%.4f", $seM ), "\t", sprintf( "%.4f", $lowBdM ), "\t",
      sprintf( "%.4f", $upBdM ), "\t";
    
    #case meanSm - meanLa <0 => meanLa greater!
    if ( $lowBdM < 0 && $upBdM < 0 ) {
        print ST "Better in\t$bigger_marker\t";
        print ST "mean\t", sprintf( "%.4f", abs($statM->mean()) ), "\t",
        sprintf( "%.4f", $seM ), "\t", sprintf( "%.4f", abs($upBdM)), "\t",
        sprintf( "%.4f", abs($lowBdM) ), "\t";
    }
    #case meanSm - meanLa >0 => meanSm greater!
    elsif ( $lowBdM > 0 && $upBdM > 0 ) {
        print ST "Better in\t$smaller_marker\t";
        print ST "mean\t", sprintf( "%.4f", $statM->mean() ), "\t",
        sprintf( "%.4f", $seM ), "\t", sprintf( "%.4f", $lowBdM ), "\t",
        sprintf( "%.4f", $upBdM ), "\t";
    }
    #case insignif:
    else {
        print ST "insignificant\n";
    }
    print ST "\n$round\t";
    print ST "$dir\t";
    print ST "R\t", sprintf( "%.4f", $statR->mean() ), "\t",
      sprintf( "%.4f", $seR ), "\t", sprintf( "%.4f", $lowBdR ), "\t",
      sprintf( "%.4f", $upBdR ), "\t";
    #case meanSm - meanLa <0 => meanLa greater!
    if ( $lowBdR < 0 && $upBdR < 0 ) {
        print ST "Better in\t$bigger_marker\t";
        print ST "R\t", sprintf( "%.4f", abs($statR->mean()) ), "\t",
        sprintf( "%.4f", $seR ), "\t", sprintf( "%.4f", abs($upBdR)), "\t",
        sprintf( "%.4f", abs($lowBdR) ), "\t";
    }
    #case meanSm - meanLa >0 => meanSm greater!
    elsif ( $lowBdR > 0 && $upBdR > 0 ) {
        print ST "Better in \t$smaller_marker\t";
        print ST "R\t", sprintf( "%.4f", $statR->mean() ), "\t",
        sprintf( "%.4f", $seR ), "\t", sprintf( "%.4f", $lowBdR ), "\t",
        sprintf( "%.4f", $upBdR ), "\t";
    }
    #case insignif:
    else {
        print ST "insignificant\n";
    }
    print ST "\n";
   # print "output to $stat_outputfile finished.\n\n";
}

sub do_ltt{
    my $name = shift @_;
    #print $name;
    my (@ltt)=@_;
    my $ctr=0;
    open OUT,">",$name."ltt_out.txt";
    foreach (@ltt){
        
        print OUT $ctr,"\t";
        print OUT @{$_}[1],"\t";
        print OUT @{$_}[2],"\n";
        $ctr++;
    }
    close OUT;
}

 sub convert_seconds_to_hhmmss {
  my $hourz=int($_[0]/3600);
  my $leftover=$_[0] % 3600;
  my $minz=int($leftover/60);
  my $secz=int($leftover % 60);
  return sprintf ("%02d:%02d:%02d", $hourz,$minz,$secz)
 }
 
sub table{
    my $first=1;#first file
    my @table;
    my $head="Age";
    print "\n\nlength of \@\_ ", $#_,"\n";
    for(my $i=0;$i<=$#_;$i++){
        my @infile=$_[$i];
        foreach my $file (@infile){
            my $part=$file;
            $part=~s/RAxML\_bipartitions\.Alignment\.phylip\.partition\_pergene\.txt\.excluding//i;
            $part=~s/\.txt\.bip\.tre\.list//;
            my $countr=0;
            open IN,"<",$file;
            while (<IN>){
             
               #first file defines age and leafs, then bips
               if (/Age,BS,node/gi){
                    next ;
               }
               if (/(\d+.\d+)\,(\d+)\,(.+)/gi){
                    if($first){
                        #my $string="$1,$3,$2";
                        my $string="$1,$2";
                        push @table, $string;
                        
                    }else{
                        #$head=$head.$part;
                        $table[$countr]=$table[$countr].",".$2;
                        #print $table[$countr],"\n";
                        
                    }
                    $countr++;
               }
            }
        $first=0;
        $head=$head.",".$part;
        }
    }
    
    unshift @table,$head;
    @table;
    
}

# summarize per replicate bootstrap support, prints a table with terminals, bs[1..N]
sub sub_rep_stat{               #slowest possible approach, but easy to make. OPTIMIZE!! e.g. open all files first
    my @grepd=grep /\S*/,@_;
    #print "\n",$#grepd+1;
    #get bootstraps from each tree, the calc stdev
    my $treefile=$grepd[0];     #expect all trees to have same taxa
    my $tree;
    open FIRST, $treefile;
    while(<FIRST>){$tree.=$_};
    close FIRST;
    my $first_forest=Bio::Phylo::IO->parse(-format=>'newick',-string => $tree);
    my $first_tree=$first_forest->first;
    if($first_tree->is_ultrametric(0.01)){
        print "\t\t$first_tree is ultrametric\n";
    
        }else{
        print "\t\t$first_tree is not ultrametric\n";
    }
    open OUT,">","$treefile.repstatistics";
    foreach my $taxnode (@{$first_tree->get_entities}){
        if($taxnode->is_terminal){
            print "\nnode is leaf: ",$taxnode->get_name;
        }elsif($taxnode->is_root){
            print "\n skip root\n";
        }else{
            print "\nlooking at node\n", $taxnode->get_name;
            my @bpterm=@{$taxnode->get_terminals};          #reference to terminals in first tree
            my @derefterm;
            foreach (@bpterm){
                                push @derefterm, $_->get_name;
                           }  
            #print "\nnode has these terminals",;
            print OUT join " ",@derefterm,",";
            #my @locterm;
            foreach my $reptile(@grepd){
                #print "good";
                    if($reptile){
                    open REP, $reptile or die $!;
                    while(<REP>){
                        
                            my $forest=Bio::Phylo::IO->parse(-format=>'newick',-string => $_);
                            my $rep=$forest->first;
                                print "\.";
                            #print $rep,"\n";
                             foreach my $node(@{$rep->get_entities}){
                                my @bplocal=@{$node->get_terminals};
                                #get node data then.
                                my @derefbp;
                                foreach (@bplocal){
                                     push @derefbp, $_->get_name;
                                }
                                #print @derefbp," comp with ",@derefterm," \n ";
                                my $bs=$node->get_name;
                                my $lc = List::Compare->new(\@derefterm, \@derefbp);      #compare two array irrespective order
                                if ($lc->is_LeqvlntR ){                 #test: both lists are equivalent
                                    print OUT $bs.", ";
                                }
                            }
                           
                        }
                        close REP;
                    }else{
                        print "\nempty line\n";
                        }
                }
            print OUT "\n";
        }
    }
    close OUT;
}
sub get_info{
    my $file=$_[0];
    my %out;
    open TEMPLATE, $file;
    while(<TEMPLATE>){
        chomp;
        next if (/\#/ or m/^\s*$/);
        
        /\t*(\S+)\s*\=\s*\"*(\S+)\b\"*\s*/;
        
        $out{$1}=$2;
        #print $1,"----->",$2,"\n";
    }
    %out;
}
__END__

