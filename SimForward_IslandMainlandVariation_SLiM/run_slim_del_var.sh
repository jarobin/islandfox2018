# Submission script for forward simulations in SLiM to examine levels of neutral and 
# deleterious mutations in a small island relative to a large mainland population.
# Implemented with SLiM version 2.4.2
# 
# This script writes a SLiM job script and then executes it.
#
# SLiM must be installed. Change paths as necessary.
#
# Scripts to make jobs specific to each model are required (e.g. make_slim_split.job.sh)
#
# Usage: ./run_slim_del_var.sh

SLIMDIR=/opt/SLiM/bin

# Set number of genes
# Note: the number of genes per "chromosome" is hard-coded in the SLiM job scripts,
# so if a number other than 2000 is used here, the numbers of genes per chromosome
# in make_slim_${m}.job.sh must also be changed
g=2000

# Set number of burn-in generations to build up variation in the ancestral population
t=100000

# Set per-site mutation rate
u=1e-8

# Set long-term population size on the island (e.g. 50, 100, 200, 500, 1000)
n=1000

# Set dominance coefficient for deleterious mutations 
# (e.g. h=0.0 for fully recessive, h=0.5 for additive)
h=0.5

# Model type (e.g. split, bigisland, bottancient, bottrecent, bottserial, inbred)
m=bigisland

# Make the SLiM script
./make_slim_${m}.job.sh ${g} ${t} ${u} ${n} ${h}

# Execute the simulation
${SLIMDIR}/slim slim_${m}_${g}genes_${t}burn_n${n}_h${h}.job > \
slim_${m}_${g}genes_${t}burn_n${n}_h${h}.job_out.txt
