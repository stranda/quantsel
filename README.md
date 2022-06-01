# quantsel

This package implements 2d continuous space genetic simulations.  It
is roughly based on the old package kernelPop, but has so many
innovations, that it has been renamed.  Here is a summary of features:

  * continuous space
  * discrete-time
  * 1000s of Mendelian loci possible (but not really intended as a genomics simulator)
  * haploid as well as diploid inheritance
  * additive polygenic phenotypes (specified by user upon Mendelian loci above)
  * plasticity
  * mapping phenotypes to both dispersal characteristics and fitness components
	* implement selection on phenotypes
	* allows dispersal to evolve
  * user-specified arbitrary life history complexity
  * density dependent population regulation (tolerance to crowding is a fitness component)
  * dispersal is broken into two components, short and long, with a mixture parameter
  * grid-based habitat quality layer, each habitat can specify different:
	* carrying capacity (survival decreases as N increases based on logistic function)
	* plasticity effect on each phenotype
	* selection gradients for each phenotype

