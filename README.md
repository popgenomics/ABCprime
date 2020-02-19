DILSprime  
=================
   * [Model](#model)
   * [Usage](#usage)

# Model  
**DILSprime** is made to investigate various models of speciation by using msprime as simulator.  
Currently, only models with two populations are handled by DILSprime : **SI**, **AM**, **SC** and **IM**.  

# Usage  
priorgen_2pop.py SC_2M_1N 50 config.yaml | python simMSprime.py SC 250000 | mscalc_2pop_SFS.py 0  
  
priorgen generates prior distribution for the Hudson's *ms* coalescent simulator, in coalescent units used by *ms*. The first argument is a model specifier (here: **S**econdary **C**ontact, with genomic heterogeneity for migration and homogeneity for *Ne*). The second argument is the number of multilocus simulations. And the last argument contains the boundary of the prior distribution. 
simMSprime simulate samples by using msprime. The first argument is the demographic model (in [SI, AM, SC, IM]) and the second argument is the effective population size of the (virtual) reference population used by priorgen. This last argument is used to convert the coalescent units (produced by priorgen) in demographic units (used by simMSprime).  
mscalc produces summary statistics (ABCstat.txt) and join site-frequency-spectrum (ABCjsfs.txt) from the ms's and msprime's output.  


