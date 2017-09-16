#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::Phylo::IO qw(parse);
use Bio::TreeIO;
use Bio::Tree::NodeI;

use List::Compare;

#open tree with bipartitions from raxml plotted on topology from complete data
#open tree 
my ($bptreein,$toptreein,$bptreestring,$toptreestring,$replicates);
if(@ARGV){
    $toptreein=$ARGV[0];
    $bptreein=$ARGV[1];
   
}
my $list=$bptreein.".list";
print "executing bs_age_gradient with \n\t$toptreein and \t$bptreein\n";
open TOP,"<",$toptreein;
while(<TOP>){$toptreestring.=$_ };
open TREE,"<",$bptreein or die $!;
while(<TREE>){$bptreestring.=$_ };          #concatenates all lines in file!! must be newick
#now TOP handles the topology and ages
#now TREE handles the bootstraps

#get topology with Bio::Phylo(can handle age);
my $toptreeforest=Bio::Phylo::IO->parse(-format=>'newick',-string => $toptreestring);
my $toptree=$toptreeforest -> first;

if ($toptree->is_ultrametric(0.01)){
    print "\t\t$toptreein is ultrametric\n" ;
    print "and has ",$toptree->calc_number_of_internals," nodes.";
}else{
    print "\t\t$toptreein is not ultrametric\n" ;
}
$toptree->calc_node_ages;
#get bootstraps from raxml 
my $bptreeforest=Bio::Phylo::IO->parse(-format=>'newick',-string => $bptreestring);
my $bptree=$bptreeforest -> first;
my @nodes=$bptree->get_nodes();
my @taxa=("Sphagnum_capillifolium","Takakia_lepidozioides");
my @BIN;
foreach (@nodes){
    my ($out);
    if($_->is_terminal){
        
            $out=$_->get_name;
           
            if ($out eq $taxa[0]){
                print "\nfound $out";
                push @BIN,$_;
            }elsif($out eq $taxa[1]){
                print "\nfound $out";
                push @BIN,$_;
            }
    }
}
print @BIN;
my $outnode=$bptree->get_mrca(\@BIN) ;
print "\nThe_outnode:$outnode";
#reroot
$bptree->reroot($outnode);


#$bptree->set_root_node($node);
if($bptree->is_ultrametric(0.01)){
    print "\t\t$bptreein is ultrametric\n";
    
}else{
    print "\t\t$bptreein is not ultrametric\n";
}
open OUT,">",$list;         #needs something unique from alignment.part.tre
print OUT "Age,BS,node(Terminals)\n";
my $n;
my @biglist;
# cycle through all nodes from topo
my $countr; my $countr2;
 foreach my $topnode ( @{ $toptree->get_entities } ) {

    if( $topnode->is_terminal ){
            #print "\nNode is leaf: ",$topnode->get_name,"\n";
    }elsif($topnode->is_root){
            print "Rootage: ",$topnode->get_generic('age'),"\n";
    }else{
        $countr++;
            my @term=@{$topnode->get_terminals};
            my $age=$topnode->get_generic('age');
            my @topterm;
            foreach(@term){
                push @topterm, $_->get_name;
            }
            #print "\n\nTopterm: ",@topterm;
        #cycle through second treee with bootstrapvalues
        #put here: out loop cycles
            foreach my $bpnode(@{$bptree->get_entities}){
                if ($bpnode->is_terminal){
                   # print "\nNode is leaf: ",$bpnode->get_name,"\n";  
                }else{
                my @bpterm=@{$bpnode->get_terminals};
                my @locterm;
                foreach(@bpterm){
                    #print $_->get_name;
                    push @locterm,$_->get_name;             #get name of the terminal node
                    }
                #get bootstrap
                    my $bs = $bpnode->get_name;             #'name' stores bootstrap here (only raxml??)
                #identify nodes by taxa
                    my $lc = List::Compare->new(\@locterm, \@topterm);      #compare two array irrespective order
                    if ($lc->is_LeqvlntR ){                 #test: both lists are equivalent
                        #print "\nEqual:\nLocterm ",@locterm,"\nTopterm ",@topterm ;
                        #print "\nThe Support: ",$bs," and the age: ",$age;
                        print OUT "$age,$bs,",join(" ",@locterm),"\n";
                    $countr2++;
                    }
                }
                
            }   
        }
 }
 print "\nNodes in $toptree: $countr\n ";
 print "found in $bptree: $countr2\n";
 