#!/usr/bin/python
# before executing on IFB cluster: module load msprime/0.7.3

import msprime
import math
import numpy as np
import sys

help = "\t\033[1;31;40mExecute msprime for different demographic models by reading from the stdin produced by priorgen. The output respect the Hudson's ms output format\033[0m\n\n\t"
help += "\033[1;31;40mTakes one model specifier and the effective population size of the reference population used by priorgen (expressed in number of diploid individuals):\033[0m\n\t\t"
help += "\n\t\t".join(['SI', 'AM', 'SC', 'IM'])
help += "\n\n"
help += "\t\033[1;32;40mExample: ./simMSprime.py SC 250000\033[0m\n"

if len(sys.argv) != 3:
	print(help)
	sys.exit()

model = sys.argv[1]
Nref = int(sys.argv[2])
plot = False # if True: plot the different trees along a gene

infile = open('bpfile', 'r')
tmp = infile.readline()
L = [ float(i) for i in infile.readline().strip().split('\t') ]
nA = [ int(i) for i in infile.readline().strip().split('\t') ]
nB = [ int(i) for i in infile.readline().strip().split('\t') ]
theta = [ float(i) for i in infile.readline().strip().split('\t') ]
rho = [ float(i) for i in infile.readline().strip().split('\t') ]

nLoci = len(L)

mu_rate = [ theta[i]/(4.0 * Nref * L[i]) for i in range(nLoci) ]
r_rate = [ rho[i]/(4.0 * Nref * L[i]) for i in range(nLoci) ]

def msOutput(ts):
	sim = ts['msprime']
	L = sim.get_sequence_length()
	segsites = sim.get_num_mutations()
	nInd = sim.get_sample_size()
		
	position = []
	genotypes = []
	for i in sim.variants():
		position.append(i.site.position / L)
		genotypes.append(i.genotypes)

	res = '// {0}\n'.format('\t'.join([ str(i) for i in ts['parameters'] ]))
	res += 'segsites: {0}\n'.format(segsites)
	if segsites > 0:
		res += 'positions: {0}\n'.format('\t'.join( [ str(i) for i in position ] ))
		tmp = ''
		for ind in range(nInd):
			for pos in range(segsites):
				tmp += '{0}'.format(genotypes[pos][ind])
			tmp += '\n'
		res += tmp 
	else:
		res += ''
	return(res)

def SI( nsamA=5, nsamB=5, length=1000, mutation_rate=0.000000003, recombination_rate=0.00000003, N1=1000, N2=1000, Na=2000, Tsplit=5000 ):
	# Tsplit : time of divergence in generations
	# N1; N2 : effective (diploid) population size
	res = {}
	res['parameters'] = [ length, mutation_rate, recombination_rate, nsamA, nsamB, N1, N2, Na, Tsplit ]

	# Migration
	M = np.array([ [0, 0], [0, 0] ]) # using ms annotation : [ [m 1 1, m 1 2], [m 2 1, m 2 2] ]
	
	# 2 populations
	pop0 = msprime.PopulationConfiguration(sample_size = nsamA, initial_size = N1, growth_rate = 0.00)
	pop1 = msprime.PopulationConfiguration(sample_size = nsamB, initial_size = N2, growth_rate = 0.00)

	# divergence
	divergence_event = msprime.MassMigration(time = Tsplit, source = 1, dest = 0, proportion = 1) # mass migration from 1 to zero; 'source' and 'dest' have a backward definition; here, all lineages from population 1 move backward to population 0.
	
	# no migration prior the Time of split
	divergence_rate_change = msprime.MigrationRateChange(time = Tsplit, rate = 0, matrix_index=None)
	
	# set the ancestral population size	
	ancestral_pop = msprime.PopulationParametersChange(time = Tsplit, initial_size = Na, population = 0 )
	
	ts = msprime.simulate(population_configurations = [pop0, pop1], length=length, mutation_rate=mutation_rate, recombination_rate=recombination_rate, migration_matrix = M, 
		demographic_events = [divergence_event, divergence_rate_change, ancestral_pop ])
	res['msprime'] = ts

	if plot == True:
		tree = ts.first()
		print( tree.draw( format = 'unicode' ) )
		for tree in ts.trees():
			print('interval:', tree.interval)
			print(tree.draw(format='unicode'))
	return(res)

def IM( nsamA=5, nsamB=5, length=1000, mutation_rate=0.000000003, recombination_rate=0.000000003, N1=1000, N2=1000, Na=2000, Tsplit=5000, m_1_2=0.01, m_2_1=0.05 ):
	# Tsplit : time of divergence in generations
	# N1; N2 : effective (diploid) population size
	res = {}
	res['parameters'] = [ length, mutation_rate, recombination_rate, nsamA, nsamB, N1, N2, Na, Tsplit, m_1_2, m_2_1 ]

	# Migration
	M = np.array([ [0, m_1_2], [m_2_1, 0] ]) # using ms annotation : [ [m 1 1, m 1 2], [m 2 1, m 2 2] ]
	
	# 2 populations
	pop0 = msprime.PopulationConfiguration(sample_size = nsamA, initial_size = N1, growth_rate = 0.00)
	pop1 = msprime.PopulationConfiguration(sample_size = nsamB, initial_size = N2, growth_rate = 0.00)

	# divergence
	divergence_event = msprime.MassMigration(time = Tsplit, source = 1, dest = 0, proportion = 1) # mass migration from 1 to zero; 'source' and 'dest' have a backward definition; here, all lineages from population 1 move backward to population 0.
	
	# no migration prior the Time of split
	divergence_rate_change = msprime.MigrationRateChange(time = Tsplit, rate = 0, matrix_index=None)
	
	# set the ancestral population size	
	ancestral_pop = msprime.PopulationParametersChange(time = Tsplit, initial_size = Na, population = 0 )
	
	ts = msprime.simulate(population_configurations = [pop0, pop1], length=length, mutation_rate=mutation_rate, recombination_rate=recombination_rate, migration_matrix = M, 
		demographic_events = [divergence_event, divergence_rate_change, ancestral_pop ])
	res['msprime'] = ts
	if plot == True:
		tree = ts.first()
		print( tree.draw( format = 'unicode' ) )
		for tree in ts.trees():
			print('interval:', tree.interval)
			print(tree.draw(format='unicode'))
	return(res)

def AM( nsamA=5, nsamB=5, length=1000, mutation_rate=0.000000003, recombination_rate=0.000000003, N1=1000, N2=1000, Na=2000, Tsplit=5000, Tam=2000, m_1_2=0.01, m_2_1=0.05 ):
	# Tsplit : time of divergence in generations
	# N1; N2 : effective (diploid) population size
	# m_1_2=x : x% of population 1 made by migrants from population 2 per generation
	res = {}
	res['parameters'] = [ length, mutation_rate, recombination_rate, nsamA, nsamB, N1, N2, Na, Tsplit, Tam, m_1_2, m_2_1 ]

	# Migration
#	M = np.array([ [0, m_1_2], [m_2_1, 0] ]) # using ms annotation : [ [m 1 1, m 1 2], [m 2 1, m 2 2] ]
	M = np.array([ [0, 0], [0, 0] ])
	
	# 2 populations
	pop0 = msprime.PopulationConfiguration(sample_size = nsamA, initial_size = N1, growth_rate = 0.00)
	pop1 = msprime.PopulationConfiguration(sample_size = nsamB, initial_size = N2, growth_rate = 0.00)

	# ancestral migration
	ancestral_migration_M_1_2 = msprime.MigrationRateChange(time = Tam, rate = m_1_2, matrix_index=(0,1))
	ancestral_migration_M_2_1 = msprime.MigrationRateChange(time = Tam, rate = m_2_1, matrix_index=(1,0))
	
	# divergence
	divergence_event = msprime.MassMigration(time = Tsplit, source = 1, dest = 0, proportion = 1) # mass migration from 1 to zero; 'source' and 'dest' have a backward definition; here, all lineages from population 1 move backward to population 0.
	
	# no migration prior the Time of split
	divergence_rate_change = msprime.MigrationRateChange(time = Tsplit, rate = 0, matrix_index=None)
	
	# set the ancestral population size	
	ancestral_pop = msprime.PopulationParametersChange(time = Tsplit, initial_size = Na, population = 0 )
	
	ts = msprime.simulate(population_configurations = [pop0, pop1], length=length, mutation_rate=mutation_rate, recombination_rate=recombination_rate, migration_matrix = M, demographic_events = [ ancestral_migration_M_1_2, ancestral_migration_M_2_1, divergence_event, divergence_rate_change, ancestral_pop ])
	
	tree = ts.first()
	print( tree.draw( format = 'unicode' ) )
	res['msprime'] = ts
	return(res)

def SC( nsamA=5, nsamB=5, length=1000, mutation_rate=0.000000003, recombination_rate=0.000000003, N1=1000, N2=1000, Na=2000, Tsplit=5000, Tsc=2000, m_1_2=0.01, m_2_1=0.05 ):
	# Tsplit : time of divergence in generations
	# N1; N2 : effective (diploid) population size
	# m_1_2=x : x% of population 1 made by migrants from population 2 per generation
	res = {}
	res['parameters'] = [ length, mutation_rate, recombination_rate, nsamA, nsamB, N1, N2, Na, Tsplit, Tsc, m_1_2, m_2_1 ]

	# Migration
	M = np.array([ [0, m_1_2], [m_2_1, 0] ]) # using ms annotation : [ [m 1 1, m 1 2], [m 2 1, m 2 2] ]
#	M = np.array([ [0, 0], [0, 0] ])
	
	# 2 populations
	pop0 = msprime.PopulationConfiguration(sample_size = nsamA, initial_size = N1, growth_rate = 0.00)
	pop1 = msprime.PopulationConfiguration(sample_size = nsamB, initial_size = N2, growth_rate = 0.00)

	# ancestral migration
	ancestral_migration_M_1_2 = msprime.MigrationRateChange(time = Tsc, rate = 0, matrix_index=(0,1))
	ancestral_migration_M_2_1 = msprime.MigrationRateChange(time = Tsc, rate = 0, matrix_index=(1,0))
	
	# divergence
	divergence_event = msprime.MassMigration(time = Tsplit, source = 1, dest = 0, proportion = 1) # mass migration from 1 to zero; 'source' and 'dest' have a backward definition; here, all lineages from population 1 move backward to population 0.
	
	# no migration prior the Time of split
	divergence_rate_change = msprime.MigrationRateChange(time = Tsplit, rate = 0, matrix_index=None)
	
	# set the ancestral population size	
	ancestral_pop = msprime.PopulationParametersChange(time = Tsplit, initial_size = Na, population = 0 )
	
	ts = msprime.simulate(population_configurations = [pop0, pop1], length=length, mutation_rate=mutation_rate, recombination_rate=recombination_rate, migration_matrix = M, demographic_events = [ ancestral_migration_M_1_2, ancestral_migration_M_2_1, divergence_event, divergence_rate_change, ancestral_pop ])
	
	res['msprime'] = ts
	
	plot = False
	if plot == True:
		tree = ts.first()
		print( tree.draw( format = 'unicode' ) )
		for tree in ts.trees():
			print('interval:', tree.interval)
			print(tree.draw(format='unicode'))
	return(res)

cnt = -1
for line in sys.stdin:
	cnt += 1
	line = line.strip().split('\t')
	if model == 'SI':
		#	      0		   1	  2   3	       4   5	      6	       7       8	   9  10
		# SI : msnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN tbs tbs
		nsamA = int( line[4] )
		nsamB = int( line[5] )
		mutation_rate = mu_rate[ cnt%nLoci ]
		recombination_rate = r_rate[ cnt%nLoci ]
		length = L[ cnt%nLoci ]
		N1 = int(round(float(line[6])* Nref, 0))
		N2 = int(round(float(line[7])* Nref, 0))
		Na = int(round(float(line[10])* Nref, 0))
		Tsplit = int(round(float(line[8]) * 4 * Nref, 0))

		# check the parameters
		if N1 <= 0:
			N1 = 1
		if N2 <= 0:
			N2 = 1
		if Na <= 0:
			Na = 1
		if Tsplit <= 0:
			Tsplit = 1
	
		x = SI( nsamA=nsamA, nsamB=nsamB, length=length, mutation_rate=mutation_rate, recombination_rate=recombination_rate, N1=N1, N2=N2, Na=Na, Tsplit=Tsplit )
		print(msOutput(x))

	if model == 'AM':
		#	      0		   1	  2   3	       4   5	      6	       7        8       9   10	      11          12  13
		# AM : msnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs
		nsamA = int( line[4] )
		nsamB = int( line[5] )
		mutation_rate = mu_rate[ cnt%nLoci ]
		recombination_rate = r_rate[ cnt%nLoci ]
		length = L[ cnt%nLoci ]
		N1 = int(round(float(line[6])* Nref, 0))
		N2 = int(round(float(line[7])* Nref, 0))
	
		Tsplit = int(round(float(line[11]) * 4 * Nref, 0))
		Na = int(round(float(line[13])* Nref, 0))
		
		Tam = int(round(float(line[8]) * 4 * Nref, 0))

		# check the parameters
		if N1 <= 0:
			N1 = 1
		if N2 <= 0:
			N2 = 1
		if Na <= 0:
			Na = 1
		if Tsplit <= 0:
			Tsplit = 1
		if Tam <= 0:
			Tam = 1
	
		m_1_2 =  float(line[9]) / 4.0 / N1
		m_2_1 =  float(line[10]) / 4.0 / N2
		
		x = AM( nsamA=nsamA, nsamB=nsamB, length=length, mutation_rate=mutation_rate, recombination_rate=recombination_rate, N1=N1, N2=N2, Na=Na, Tsplit=Tsplit, Tam=Tam, m_1_2=m_1_2, m_2_1=m_2_1 )
		print(msOutput(x))

	if model == 'SC':
		#	      0		   1	  2   3	       4   5		6	   7	    8	     9	     10	       11	   12  13
		# SC : msnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs
		nsamA = int( line[4] )
		nsamB = int( line[5] )
		mutation_rate = mu_rate[ cnt%nLoci ]
		recombination_rate = r_rate[ cnt%nLoci ]
		length = L[ cnt%nLoci ]
		N1 = int(round(float(line[8])* Nref, 0))
		N2 = int(round(float(line[9])* Nref, 0))
	
		Tsplit = int(round(float(line[11]) * 4 * Nref, 0))
		Na = int(round(float(line[13])* Nref, 0))
		
		Tsc = int(round(float(line[10]) * 4 * Nref, 0))

		# check the parameters
		if N1 <= 0:
			N1 = 1
		if N2 <= 0:
			N2 = 1
		if Na <= 0:
			Na = 1
		if Tsplit <= 0:
			Tsplit = 1
		if Tsc <= 0:
			Tsc = 1
	
		m_1_2 =  float(line[6]) / 4.0 / N1
		m_2_1 =  float(line[7]) / 4.0 / N2
		
		x = SC( nsamA=nsamA, nsamB=nsamB, length=length, mutation_rate=mutation_rate, recombination_rate=recombination_rate, N1=N1, N2=N2, Na=Na, Tsplit=Tsplit, Tsc=Tsc, m_1_2=m_1_2, m_2_1=m_2_1 )
		print(msOutput(x))

	if model == 'IM':
		#	      0		   1      2   3	       4   5	      6		7	  8	     9	     10		 11  12
		# IM : msnsam tbs 10000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -n 1 tbs -n 2 tbs -m 1 2 tbs -m 2 1 tbs -ej tbs 2 1 -eN tbs tbs
		nsamA = int( line[4] )
		nsamB = int( line[5] )
		mutation_rate = mu_rate[ cnt%nLoci ]
		recombination_rate = r_rate[ cnt%nLoci ]
		length = L[ cnt%nLoci ]
		N1 = int(round(float(line[6])* Nref, 0))
		N2 = int(round(float(line[7])* Nref, 0))
	
		Tsplit = int(round(float(line[10]) * 4 * Nref, 0))
		Na = int(round(float(line[12])* Nref, 0))

		# check the parameters
		if N1 <= 0:
			N1 = 1
		if N2 <= 0:
			N2 = 1
		if Na <= 0:
			Na = 1
		if Tsplit <= 0:
			Tsplit = 1
		
		m_1_2 =  float(line[8]) / 4.0 / N1
		m_2_1 =  float(line[9]) / 4.0 / N2
		
		x = IM( nsamA=nsamA, nsamB=nsamB, length=length, mutation_rate=mutation_rate, recombination_rate=recombination_rate, N1=N1, N2=N2, Na=Na, Tsplit=Tsplit, m_1_2=m_1_2, m_2_1=m_2_1 )
		print(msOutput(x))

