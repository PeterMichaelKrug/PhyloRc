# PhyloRc

PhyloRc (or whatever it will be named after publication) is a couple of perl scripts from my third chapter of my thesis. It calculates three variants of phylogenetic structure (sensu Müller et al. 2006 (MPE)):
	
	Rc,free: the phylogenetic structure of a given partition (of an alignment) without constraining the topology
	
	Rc,fix: the phylogenetic structure of a given partition with a constraint topology (i.e. a chronogram)
	
	Rc,fix(S): the phylogenetic structure of a given partition under a constraint topology for a given subset of nodes (i.e. older or younger then the median age of the tree)

This tools uses RAxML v8.x ff (Stammatakis et al. 2014 (Bioinformatics)) to compute bootstrap samples and to calculate the frequency of nodes either on the partition's ML-tree or the constraint topology.
Since partitions (e.g. genes, introns, spacers) most likely differ in total nucleotide count, a resampling has to be done first to enable unbiased comparisons between each partitions (as much as possible). (See Krug et al. 2017 for details)

Thus, the results of Rc are presented as a list of differences between partitions.

How To:
The Main script is configured via template.txt. Here, you specify name and location of the alignment (in relaxed phylip format), the tree file and the partition file, furthermore raxml options, number of resampling replicates and bootstrap replicates.

After the template file has been configured, invoke the scripts via "perl Mainresampling.pl template.txt".
