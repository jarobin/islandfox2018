# Script to simulate San Nicolas island fox genomes under neutrality in  msprime (v0.5.0) 
# for analysis of peaks of heterozygosity.
#
# Performs 1000 simulations, each generating 220 10 Mb "chromosomes" in four individuals.
# Each simulation generates and writes to its own output folder.
#
# Requires two input files:
#    - SNI_ABC_posteriors_100.txt, containing the 100 best parameter sets for San Nicolas
#       demographic models inferred in Robinson et al. 2016, DOI:10.1016/j.cub.2016.02.062
#    - Wong2010_dogRecRates_perbp.txt, containing per-site recombination rates 
#       corresponding to 1 Mb chunks across the dog genome by Wong et al. 2010, 
#       DOI:10.1534/genetics.109.106831
#
# Note: 12 years added to all times in demographic models from 2016 study, because here
# T=0 corresponds to the year 2000, and in the 2016 study T=0 corresponded to the year
# 1988.
#
# Usage (msprime must be installed): python ./sim_SNIgenomes_msprime.py
# Implemented with python2.7


import sys
import os
import random
import msprime


# Define simulator function, which outputs a VCF file for one 10 Mb "chromosome"
def runSimulator(simNum, chrNum, seedNum, Ne4, T, Ne3, Ne1, recRate):
    # San Nicolas demographic model, moving backwards in time from the year 2000
    demographic_events=[
        msprime.PopulationParametersChange(time=40, initial_size=10, population_id=0),
        msprime.PopulationParametersChange(time=42, initial_size=Ne3, population_id=0),
        msprime.PopulationParametersChange(time=T, initial_size=Ne4, population_id=0),
        msprime.PopulationParametersChange(time=8012, initial_size=20000, population_id=0)
    ]
    # Sample one individual (two haplotypes) in 1929, two individuals in 1988, and
    # one individual in 2000
    samples=[
        msprime.Sample(population=0, time=71), 
        msprime.Sample(population=0, time=71), 
        msprime.Sample(population=0, time=12), 
        msprime.Sample(population=0, time=12), 
        msprime.Sample(population=0, time=12), 
        msprime.Sample(population=0, time=12),
        msprime.Sample(population=0, time=0),
        msprime.Sample(population=0, time=0)
    ]
    # Define parameters for the simulation
    tree_sequence=msprime.simulate(
        demographic_events=demographic_events, 
        samples=samples, 
        Ne=Ne1, 
        length=1e7, 
        recombination_rate=recRate, 
        mutation_rate=2e-8, 
        random_seed=seedNum
    )
    # Output VCF file
    with open('peak_sim'+str(seedNum)+'.vcf', 'w') as vcf_file:
        tree_sequence.write_vcf(vcf_file, 2, str(chrNum))


# Function to run the simulator, which outputs a log file and executes the simulator 
# 220 times

def executeJobs(simNum):
    # Create and move into a directory for each simulation
    if not os.path.exists('sim'+'{:04}'.format(simNum)):
        os.makedirs('sim'+'{:04}'.format(simNum))
    os.chdir('sim'+'{:04}'.format(simNum)) 
 
    # Set demographic parameters by randomly selecting a set from the input parameter file
    secure_random=random.SystemRandom()
    myParams=secure_random.choice(allParams).split(' ')
    Ne4=int(myParams[0])
    Ne3=int(myParams[2])
    Ne1=int(myParams[3])
    # Note: 12 years added to T
    T=int(myParams[1])+12

    # Write the selected parameters to a log file
    logfile=open('peak_sim'+str(simNum)+'.log', 'w')
    logfile.write("# Simulation %s\n" % simNum)
    logfile.write("# Ne4=%s\n" % Ne4)
    logfile.write("# T=%s\n" % T)
    logfile.write("# Ne3=%s\n" % Ne3)
    logfile.write("# Ne1=%s\n" % Ne1)
    logfile.write("#\n# Chromosome\tSeed\tRecombinationRate\n")

    # Simulate 220 10Mb regions ("chromosomes")
    chrNum=1
    while chrNum<221:
        # Seed number is the sim number and the chromosome number for reproducibility
        seedNum=int(str(simNum)+'{:03}'.format(chrNum))
        # Randomly select a recombination rate
        secure_random=random.SystemRandom()
        recRate=secure_random.choice(allRecRates)
        # Write seed and rec rate to log file
        logfile.write("%s\t%s\t%s\n" % (chrNum, seedNum, recRate))
        # Convert rec rate to float
        recRate=float(recRate)
        # Simulate!
        runSimulator(simNum, chrNum, seedNum, Ne4, T, Ne3, Ne1, recRate)
        chrNum+=1
    logfile.close()
    os.chdir('..')


# Read in the parameter sets file
with open('SNI_ABC_posteriors_100.txt') as f:
    allParams=f.read().splitlines()

f.close()


# Read in the recombination rates file
with open('Wong2010_dogRecRates_perbp.txt') as f:
    allRecRates=f.read().splitlines()

f.close()


# Start simNum counter and set the number of simulations to run (in this case, 10)
simNum=1
endSim=simNum+10


# Execute simulations and print progress to terminal
while simNum<endSim:
    executeJobs(simNum)
    print("Simulation %s complete" % simNum)
    simNum+=1


exit()
